#!/usr/bin/env python3
"""Encode raw CABBAGE JSONL (from scripts/fetch-cabbage.py) into compact
per-view snapshots for the browser: databases/cabbage/<view>.cbg.gz plus a
manifest.json.

CBG1 format (little-endian, then gzipped):
  magic "CBG1", u8 version, release str, u8 view_id, view name str,
  u32 rowCount, u16 columnCount, then per column in order:
    id str, u8 kind
    kind 0 (dictionary string): u32 nDict, entries freq-desc (u16 len + utf8),
      rowCount LEB128 varints into the dictionary (0 = null)
    kind 1 (numeric): rowCount float64 (NaN = null)
js/cabbage.js decodeSnapshot() is the reference reader — keep in sync.

The JSONL is streamed twice: once to count values (dictionaries, kinds),
once to emit all columns' row data together. Cells are normalised through a
per-column memo so repeated values are parsed once.
"""
import argparse
import gzip
import hashlib
import json
import re
import os
import struct
import sys
import time
from collections import Counter

RAW_DIR = os.path.join(os.path.dirname(__file__), "..", "databases", "cabbage", "raw")
OUT_DIR = os.path.join(os.path.dirname(__file__), "..", "databases", "cabbage")

# Columns stored as float64 (sortable, range-filterable); everything else is
# dictionary-encoded text.
NUMERIC_COLS = {
    "collection_year", "host_age", "isolation_latitude", "isolation_longitude",
    "region_start", "region_end", "reference_sequence_coverage",
    "reference_sequence_identity", "taxon_id",
}
# JSON-array cells joined into one display string.
LIST_COLS = {"AMR_associated_publications"}


def put_str(buf, s):
    b = s.encode("utf-8")
    buf += struct.pack("<H", len(b))
    buf += b


def put_varint(buf, v):
    while True:
        b = v & 0x7F
        v >>= 7
        if v:
            buf.append(b | 0x80)
        else:
            buf.append(b)
            return


def num_str(x):
    """Canonical display form of a float cell ('2010', not '2010.0')."""
    if x == int(x) and abs(x) < 1e15:
        return str(int(x))
    return repr(x)


def make_norm(col_suffix):
    """f(raw) -> ('s'|'n', value) | None for one column."""
    numeric = col_suffix in NUMERIC_COLS
    is_list = col_suffix in LIST_COLS

    def norm(raw):
        if raw is None:
            return None
        if is_list and isinstance(raw, list):
            return ("s", "; ".join(str(x) for x in raw))
        if isinstance(raw, bool):
            return ("s", str(raw))
        if numeric:
            if isinstance(raw, (int, float)):
                return ("n", float(raw))
            s = str(raw).strip()
            if not s:
                return None
            try:
                return ("n", float(s))
            except ValueError:
                return ("s", s)  # rare stray text in a numeric column
        if isinstance(raw, float):
            return ("s", num_str(raw))
        return ("s", str(raw)) if str(raw) else None
    return norm


def build(name):
    meta = json.load(open(os.path.join(RAW_DIR, f"{name}.jsonl.meta.json")))
    col_ids = meta["columns"]
    suffixes = [c.split("-", 1)[1] for c in col_ids]
    norms = [make_norm(s) for s in suffixes]
    path = os.path.join(RAW_DIR, f"{name}.jsonl")

    # Pass 1: row count, column kinds, string-value counts (under the
    # canonical string form, so mixed columns still get one dictionary).
    kinds = [None] * len(col_ids)
    counts = [Counter() for _ in col_ids]
    memos = [{} for _ in col_ids]  # raw -> normalised, shared with pass 2
    rows = 0
    with open(path) as f:
        for line in f:
            rows += 1
            rec = json.loads(line)
            for i, cid in enumerate(col_ids):
                raw = rec.get(cid)
                if isinstance(raw, list):
                    v = norms[i](raw)
                elif raw in memos[i]:
                    v = memos[i][raw]
                else:
                    v = norms[i](raw)
                    try:
                        memos[i][raw] = v
                    except TypeError:
                        pass  # unhashable (never happens for non-list JSON)
                if v is None:
                    continue
                if kinds[i] is None:
                    kinds[i] = v[0]
                elif kinds[i] != v[0]:
                    kinds[i] = "s"  # mixed column -> text dictionary
                counts[i][v[1] if v[0] == "s" else num_str(v[1])] += 1
    kinds = [k or "s" for k in kinds]  # all-null column -> empty text dict

    # Dictionary index assignment by descending frequency; 0 = null.
    entries_out, index = [], []
    for i in range(len(col_ids)):
        if kinds[i] != "s":
            entries_out.append(None)
            index.append(None)
            continue
        entries = [v for v, _ in counts[i].most_common()]
        entries_out.append(entries)
        index.append({v: j + 1 for j, v in enumerate(entries)})

    # Pass 2: stream the rows once, appending to per-column buffers. The
    # final layout must match decodeSnapshot(): header per column (id, kind,
    # dict entries for strings OR the f64 block for numerics) followed by
    # every string column's varint block in column order.
    head = bytearray(b"CBG1")
    head.append(1)
    put_str(head, meta["release"])
    head.append(meta["view_id"])
    put_str(head, name)
    head += struct.pack("<I", rows)
    head += struct.pack("<H", len(col_ids))

    hdrs = []      # id + kind (+ dict entries) per column
    numbufs = []   # f64 bytes per numeric column
    varbufs = []   # varint bytes per string column
    for i, cid in enumerate(col_ids):
        hb = bytearray()
        put_str(hb, cid)
        if kinds[i] == "s":
            hb.append(0)
            hb += struct.pack("<I", len(entries_out[i]))
            for e in entries_out[i]:
                put_str(hb, e)
            hdrs.append(hb)
            numbufs.append(None)
            varbufs.append(bytearray())
        else:
            hb.append(1)
            hdrs.append(hb)
            numbufs.append(bytearray())
            varbufs.append(None)

    with open(path) as f:
        for line in f:
            rec = json.loads(line)
            for i, cid in enumerate(col_ids):
                raw = rec.get(cid)
                if isinstance(raw, list):
                    v = norms[i](raw)
                else:
                    v = memos[i].get(raw)
                if kinds[i] == "s":
                    if v is None:
                        put_varint(varbufs[i], 0)
                    else:
                        sval = v[1] if v[0] == "s" else num_str(v[1])
                        put_varint(varbufs[i], index[i][sval])
                else:
                    numbufs[i] += struct.pack("<d", v[1] if (v and v[0] == "n") else float("nan"))

    out = head
    for hb, nb in zip(hdrs, numbufs):
        out += hb
        if nb is not None:
            out += nb
    for vb in varbufs:
        if vb is not None:
            out += vb

    gz = gzip.compress(bytes(out), compresslevel=9, mtime=0)
    gpath = os.path.join(OUT_DIR, f"{name}.cbg.gz")
    with open(gpath, "wb") as f:
        f.write(gz)
    md5 = hashlib.md5(gz).hexdigest()
    print(f"{name}: {rows:,} rows x {len(col_ids)} cols -> {gpath} "
          f"({len(out) / 1e6:.1f} MB raw, {len(gz) / 1e6:.1f} MB gz)", flush=True)
    return {
        "file": f"{name}.cbg.gz",
        "view_id": meta["view_id"],
        "release": meta["release"],
        "rows": rows,
        "columns": col_ids,
        "raw_bytes": len(out),
        "bytes": len(gz),
        "md5": md5,
    }


def build_associations(min_n=20):
    """Join the genotypes and phenotypes views by BioSample and precompute
    the table the results card queries: for every (species, gene,
    antibiotic), how many isolates carrying the gene tested resistant /
    intermediate / susceptible.

    gene is the AMRFinderPlus element symbol (amr_element_symbol), acquired
    genes only (element_subtype AMR): point mutations cannot be called from
    reads in openpathogen and are left out. The annotation's gene_symbol is
    not used: it is null for many AMR calls (mph(C), tet(K), fusB, dfrG ...)
    and spells one gene several ways (mecA_2, ermC', aacA-aphD).

    Per gene the table also keeps the AMRFinderPlus class and the antibiotic
    names the portal links the gene to (genotype-antibiotic_name), so the
    card can separate a mechanistic association from co-occurrence. Two
    backgrounds are kept: all AST isolates (the species antibiogram) and AST
    isolates that also have a genome (the cohort gene carriers come from).

    Phenotype per (isolate, antibiotic): Updated_phenotype_CLSI, else
    Updated_phenotype_EUCAST (both 2025 breakpoints, per the CABBAGE paper),
    else the submitted category; the first record wins. Rows with fewer than
    min_n isolates are dropped (the card never shows them).
    """
    def effective_phenotype(r, prefix):
        for col in ('Updated_phenotype_CLSI', 'Updated_phenotype_EUCAST', 'resistance_phenotype'):
            v = r.get(prefix + col)
            if v:
                return v
        return None

    # Canonical gene symbol: quotes, parentheses and _N row suffixes go, so
    # erm(C) / ermC' / ermC_1 all count as one gene. geneCanon() in
    # js/cabbage.js mirrors this.
    def canon_gene(g):
        g = re.sub(r"['\u2019()]", "", g)
        while True:
            m = re.match(r"^(.*)_(\d+)$", g)
            if not m:
                break
            g = m.group(1)
        return g

    def bucket(ph):
        if ph is None:
            return None
        if ph in ('resistant', 'non-susceptible'):
            return 'r'
        if 'intermediate' in ph:
            return 'i'
        if ph.startswith('susceptible'):
            return 's'
        return None

    # Pass 1: genotype side. BioSample -> (species, set of canonical symbols);
    # per symbol: spelling, AMRFinderPlus class and linked antibiotic names.
    print("  association: reading genotypes...", flush=True)
    isolate_genes = {}
    gene_info = {}
    with open(os.path.join(RAW_DIR, "genotypes.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            bs = r.get("genotype-BioSample_ID")
            if not bs or r.get("genotype-element_type") != "AMR" or r.get("genotype-element_subtype") != "AMR":
                continue
            sym = r.get("genotype-amr_element_symbol") or r.get("genotype-gene_symbol")
            if not sym:
                continue
            c = canon_gene(sym)
            e = isolate_genes.get(bs)
            if e is None:
                e = isolate_genes[bs] = (r.get("genotype-species"), set())
            e[1].add(c)
            gi = gene_info.get(c)
            if gi is None:
                gi = gene_info[c] = {"symbols": Counter(), "classes": Counter(), "links": set()}
            gi["symbols"][sym] += 1
            if r.get("genotype-class"):
                gi["classes"][r["genotype-class"]] += 1
            ab = r.get("genotype-antibiotic_name")
            if ab:
                gi["links"].add(ab)
    print(f"  association: {len(isolate_genes):,} genotyped isolates, {len(gene_info):,} acquired genes", flush=True)

    # Pass 2: phenotype side.
    print("  association: joining phenotypes...", flush=True)
    assoc = {}           # (species, gene, antibiotic) -> [r, i, s]
    background = {}      # (species, antibiotic) -> [r, i, s]   all AST isolates
    background_seq = {}  # (species, antibiotic) -> [r, i, s]   AST isolates with a genome
    species_counts = {}  # species -> genotyped isolates with AST
    seen_ab = {}         # BioSample -> antibiotics already counted
    with open(os.path.join(RAW_DIR, "phenotypes.jsonl")) as f:
        for line in f:
            r = json.loads(line)
            bs = r.get("phenotype-BioSample_ID")
            ab = r.get("phenotype-antibiotic_name")
            b = bucket(effective_phenotype(r, "phenotype-"))
            if not bs or not ab or not b:
                continue
            sabs = seen_ab.get(bs)
            new_isolate = sabs is None
            if new_isolate:
                sabs = seen_ab[bs] = set()
            if ab in sabs:
                continue
            sabs.add(ab)

            g = isolate_genes.get(bs)
            species = r.get("phenotype-species") or (g[0] if g else None)
            if not species:
                continue
            if new_isolate and g:
                species_counts[species] = species_counts.get(species, 0) + 1
            idx = 'ris'.index(b)
            e = background.get((species, ab))
            if e is None:
                e = background[(species, ab)] = [0, 0, 0]
            e[idx] += 1
            if g:
                e = background_seq.get((species, ab))
                if e is None:
                    e = background_seq[(species, ab)] = [0, 0, 0]
                e[idx] += 1
                for gene in g[1]:
                    key = (species, gene, ab)
                    ae = assoc.get(key)
                    if ae is None:
                        ae = assoc[key] = [0, 0, 0]
                    ae[idx] += 1

    kept = {k: e for k, e in assoc.items() if sum(e) >= min_n}
    used_genes = sorted({g for (_, g, _) in kept})
    print(f"  association: {len(assoc):,} (species, gene, antibiotic) rows, {len(kept):,} with n >= {min_n}", flush=True)

    out = {
        "version": 2,
        "min_n": min_n,
        "rows": [
            {"species": sp, "gene": g, "antibiotic": ab, "r": e[0], "i": e[1], "s": e[2]}
            for (sp, g, ab), e in sorted(kept.items())
        ],
        "background": [
            {"species": sp, "antibiotic": ab, "r": e[0], "i": e[1], "s": e[2]}
            for (sp, ab), e in sorted(background.items())
        ],
        "background_seq": [
            {"species": sp, "antibiotic": ab, "r": e[0], "i": e[1], "s": e[2]}
            for (sp, ab), e in sorted(background_seq.items())
        ],
        "genes": [
            {
                "gene": g,
                "symbol": gene_info[g]["symbols"].most_common(1)[0][0],
                "class": (gene_info[g]["classes"].most_common(1) or [("", 0)])[0][0],
                "links": sorted(gene_info[g]["links"]),
            }
            for g in used_genes
        ],
        "species": sorted(species_counts, key=lambda s: -species_counts[s]),
    }
    payload = json.dumps(out, separators=(",", ":")).encode()
    gz = gzip.compress(payload, compresslevel=9, mtime=0)
    path = os.path.join(OUT_DIR, "associations.json.gz")
    with open(path, "wb") as f:
        f.write(gz)
    md5 = hashlib.md5(gz).hexdigest()
    print(f"associations: {len(out['rows']):,} rows, {len(out['background']):,} background, "
          f"{len(out['genes']):,} genes, {len(out['species'])} species -> {path} ({len(gz) / 1e6:.2f} MB gz)", flush=True)
    return {
        "file": "associations.json.gz",
        "version": 2,
        "min_n": min_n,
        "rows": len(out["rows"]),
        "background_rows": len(out["background"]),
        "genes": len(out["genes"]),
        "species": len(out["species"]),
        "bytes": len(gz),
        "md5": md5,
    }

def main():
    ap = argparse.ArgumentParser()
    ap.add_argument("--views", default="phenotypes,genotypes,combined")
    ap.add_argument("--associations-only", action="store_true",
                    help="rebuild associations.json.gz and its manifest entry; keep the views")
    args = ap.parse_args()

    os.makedirs(OUT_DIR, exist_ok=True)
    if args.associations_only:
        man_path = os.path.join(OUT_DIR, "manifest.json")
        if not os.path.exists(man_path):
            sys.exit("manifest.json missing: build the views first")
        manifest = json.load(open(man_path))
        manifest["built"] = time.strftime("%Y-%m-%d")
        manifest["associations"] = build_associations()
        with open(man_path, "w") as f:
            json.dump(manifest, f, indent=1)
        print("manifest updated")
        return
    # Merge with any existing manifest so partial rebuilds keep other views.
    man_path = os.path.join(OUT_DIR, "manifest.json")
    manifest = {"format": "CBG1", "built": time.strftime("%Y-%m-%d"), "views": {}}
    if os.path.exists(man_path):
        try:
            old = json.load(open(man_path))
            if old.get("format") == "CBG1":
                manifest["views"].update(old.get("views", {}))
        except (ValueError, OSError):
            pass
    for name in args.views.split(","):
        manifest["views"][name] = build(name)
        manifest["release"] = manifest["views"][name]["release"]
    if os.path.exists(os.path.join(RAW_DIR, "genotypes.jsonl")):
        manifest["associations"] = build_associations()
    with open(man_path, "w") as f:
        json.dump(manifest, f, indent=1)
    print(f"manifest written ({sum(v['bytes'] for v in manifest['views'].values()) / 1e6:.1f} MB total)")


if __name__ == "__main__":
    sys.exit(main())
