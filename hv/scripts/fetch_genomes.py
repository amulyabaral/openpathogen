#!/usr/bin/env python3
"""Fetch one representative RefSeq complete genome per organism (NCBI eutils).

Writes genomes into genomes_ref/ and the build-db manifest (manifest_ref.tsv).
Genus for rank-collapse is parsed from the organism name.
"""
import json
import sys
import time
import urllib.request
import urllib.parse
from pathlib import Path

EUTILS = "https://eutils.ncbi.nlm.nih.gov/entrez/eutils"

# WHO priority pathogens + common commensals + a control
ORGANISMS = [
    "Staphylococcus aureus",
    "Staphylococcus epidermidis",
    "Enterococcus faecalis",
    "Enterococcus faecium",
    "Escherichia coli",
    "Klebsiella pneumoniae",
    "Pseudomonas aeruginosa",
    "Acinetobacter baumannii",
    "Salmonella enterica",
    "Streptococcus pyogenes",
    "Streptococcus pneumoniae",
    "Streptococcus agalactiae",
    "Haemophilus influenzae",
    "Neisseria gonorrhoeae",
    "Neisseria meningitidis",
    "Campylobacter jejuni",
    "Clostridioides difficile",
    "Bacteroides fragilis",
    "Fusobacterium nucleatum",
    "Corynebacterium diphtheriae",
    "Listeria monocytogenes",
    "Bacillus subtilis",
    "Moraxella catarrhalis",
    "Helicobacter pylori",
    "Mycobacterium tuberculosis",
]


def get(url: str) -> bytes:
    for attempt in range(4):
        try:
            req = urllib.request.Request(url, headers={"User-Agent": "openpathogen-hv/0.1"})
            with urllib.request.urlopen(req, timeout=60) as r:
                return r.read()
        except Exception as e:  # noqa: BLE001
            print(f"  retry {attempt + 1}: {e}", file=sys.stderr)
            time.sleep(2 + 2 * attempt)
    raise RuntimeError(f"failed: {url}")


def best_refseq(organism: str) -> dict | None:
    """Plain-text esearch (field qualifiers misbehave on db=assembly),
    then pick the first summary that is a complete RefSeq genome,
    preferring reference/representative category."""
    q = urllib.parse.quote(organism)
    data = json.loads(get(f"{EUTILS}/esearch.fcgi?db=assembly&retmode=json&retmax=25&term={q}"))
    uids = data["esearchresult"]["idlist"]
    if not uids:
        return None
    summ = json.loads(get(f"{EUTILS}/esummary.fcgi?db=assembly&retmode=json&id={','.join(uids)}"))
    res = summ["result"]
    candidates = []
    for uid in uids:
        r = res[uid]
        acc = r.get("assemblyaccession", "")
        if not acc.startswith("GCF"):
            continue
        if r.get("assemblystatus") != "Complete Genome":
            continue
        if not r.get("ftppath_refseq"):
            continue
        cat = r.get("refseq_category", "")
        prio = 0 if "reference" in cat.lower() or "representative" in cat.lower() else 1
        candidates.append((prio, r))
    candidates.sort(key=lambda x: x[0])
    return candidates[0][1] if candidates else None


def main() -> None:
    out_dir = Path(__file__).resolve().parent.parent / "genomes_ref"
    out_dir.mkdir(exist_ok=True)
    manifest = out_dir.parent / "manifest_ref.tsv"
    rows = []
    for org in ORGANISMS:
        print(f"[fetch] {org}", flush=True)
        r = best_refseq(org)
        if r is None:
            print(f"  !! no complete RefSeq found, skipping", file=sys.stderr)
            continue
        acc = r["assemblyaccession"]
        asmname = r["assemblyname"]
        org_name = r["organism"]
        ftp = r["ftppath_refseq"]
        fname = f"{acc}_{asmname}_genomic.fna.gz"
        dest = out_dir / fname
        if dest.exists() and dest.stat().st_size > 0:
            print(f"  cached {acc}", flush=True)
        else:
            url = f"{ftp}/{fname}"
            data = get(url)
            dest.write_bytes(data)
            print(f"  {acc} ({r.get('refseq_category', '?')}): {len(data) / 1e6:.1f} MB", flush=True)
        genus = org_name.split(" ")[0]
        short = "".join(w[:3] for w in org_name.split()[0:2]).lower()
        rows.append((str(dest), short, org_name, genus))
        time.sleep(0.4)

    with manifest.open("w") as fh:
        fh.write("#fasta\tid\tname\tgenus\tflags\n")
        for path, sid, name, genus in rows:
            fh.write(f"{path}\t{sid}\t{name}\t{genus}\t0\n")
    print(f"[fetch] manifest → {manifest} ({len(rows)} genomes)")


if __name__ == "__main__":
    main()
