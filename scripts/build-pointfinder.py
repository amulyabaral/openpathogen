#!/usr/bin/env python3
"""build-pointfinder.py - openpathogen's PointFinder files, from a checkout
of the PointFinder database (https://bitbucket.org/genomicepidemiology/pointfinder_db).

  python3 scripts/build-pointfinder.py <pointfinder_db dir> --amrfinder <AMRProt-mutation.tsv>
                                       [--out databases/pointfinder]

--amrfinder is the point mutation table of the AMRFinderPlus database
(AMRProt-mutation.tsv in https://ftp.ncbi.nlm.nih.gov/pathogen/Antimicrobial_resistance/AMRFinderPlus/database/,
public domain; the files here were built with 2026-08-07.1). A row whose
AMRFinderPlus symbol is in it for the species is marked "amrfinder": NCBI
curates the same mutation. PointFinder marks every Klebsiella row as under
development; js/pointfinder.js reports such a row as found only when
AMRFinderPlus curates it too.

For every species folder in the database's config (P. falciparum is left out:
openpathogen screens bacteria), databases/pointfinder/<folder>/ gets
  <folder>.comp.b.gz .length.b .name .seq.b
      a KMA 1.5.1 index (kma/kma_index) of <folder>.fsa, the folder's combined
      FASTA. INSTALL.py indexes every *.fsa of the folder, which puts each
      template in twice, and a few per-gene copies differ from the combined
      one (E. coli parC, N. gonorrhoeae pilQ; M. tuberculosis's
      ethA-promoter file holds embA), so only the combined file is indexed.
      The hash table (.comp.b) has at least 2^20 slots whatever the database
      size, about 4.3 MB that is nearly all empty, so it ships gzipped (170
      to 410 KB) and the browser unpacks it; gzip runs without a timestamp,
      so the files are the same on every build.
  catalogue.json.gz
      the templates and the active mutation rows of phenotypes.txt, normalised
      for js/pointfinder.js: what kind of change each row is, where it sits on
      its template, the antibiotics, the required mutations. Every row is
      checked against its template; a row whose reference does not match is
      left out with the reason in "dropped", unless OVERRIDES fixes it.
      Templates that are not of their species are replaced (REPLACE), and
      their rows checked against the replacement.
and databases/pointfinder/manifest.json names the database version and commit
and lists the species with their files.

Row positions (phenotypes.txt, as ResFinder 4.7.2 reads them):
  codon n of a coding template      template bases P+3n-2 .. P+3n
  nucleotide p of an RNA gene, or   template base P+p
    of the gene part of a promoter
  promoter position -k              template base P+1-k
where P is the promoter length in genes.txt (ampC-promoter-size-53bp), else 0.
"""
import argparse, gzip, json, os, re, shutil, subprocess, sys, tempfile

HERE = os.path.dirname(os.path.abspath(__file__))
REPO = os.path.normpath(os.path.join(HERE, '..'))
KMA_INDEX = os.path.join(REPO, 'kma', 'kma_index')

SKIP_FOLDERS = {'plasmodium_falciparum'}
# AMRFinderPlus taxgroups (AMRProt-mutation.tsv) of each folder.
TAXGROUPS = {
    'campylobacter': ['Campylobacter'], 'enterococcus_faecalis': ['Enterococcus_faecalis'],
    'enterococcus_faecium': ['Enterococcus_faecium'], 'escherichia_coli': ['Escherichia'],
    'helicobacter_pylori': ['Helicobacter_pylori'], 'klebsiella': ['Klebsiella_pneumoniae'],
    'mycobacterium_tuberculosis': [], 'neisseria_gonorrhoeae': ['Neisseria_gonorrhoeae'],
    'salmonella': ['Salmonella'], 'staphylococcus_aureus': ['Staphylococcus_aureus'],
}
CAUTION_RE = re.compile(r'under development|used? with caution', re.I)
LABELS = {
    'campylobacter': 'Campylobacter (C. jejuni, C. coli)',
    'enterococcus_faecalis': 'Enterococcus faecalis',
    'enterococcus_faecium': 'Enterococcus faecium',
    'escherichia_coli': 'Escherichia coli',
    'helicobacter_pylori': 'Helicobacter pylori',
    'klebsiella': 'Klebsiella pneumoniae complex',
    'mycobacterium_tuberculosis': 'Mycobacterium tuberculosis',
    'neisseria_gonorrhoeae': 'Neisseria gonorrhoeae',
    'salmonella': 'Salmonella',
    'staphylococcus_aureus': 'Staphylococcus aureus',
}

# Fixes to database rows, by (folder, line number in phenotypes.txt):
#   template   a typo in the template name
#   pos        the row is written in E. coli numbering against the species'
#              own template; the local position of the same site
OVERRIDES = {
    ('campylobacter', 12): {'template': 'gyrA_2_CP035927.1', 'why': "typo gyrA_2-CP035927.1"},
    ('staphylococcus_aureus', 106): {'pos': 2527, 'why': 'T2500A is E. coli numbering; the site is 2527 on this S. aureus 23S template'},
    ('staphylococcus_aureus', 107): {'pos': 2603, 'why': 'G2576T is E. coli numbering; the site is 2603 on this S. aureus 23S template'},
    ('salmonella', 58): {'pos': 1191, 'why': 'C1192T is E. coli numbering; the site is 1191 on this Salmonella 16S template'},
}

# Templates replaced by the species' own gene, from
# scripts/pointfinder-references.fsa (each header gives its source).
# PointFinder's Klebsiella gyrA and parC carry K. pneumoniae accessions
# (AF052258.1 and AF303641.1, short QRDR sequences), but the full-length
# templates are E. coli genes: 97.9% and 98.8% identical to E. coli and
# about 85% to K. pneumoniae over their whole length, the QRDR included.
# KMA puts K. pneumoniae reads on them at about 65% template identity,
# below the 80% a gene needs, so gyrA S83/D87 and parC S80/E84 were never
# read in K. pneumoniae. The type strain's genes (ATCC 13883, 99.6% and
# 99.8% identical to K. pneumoniae genomes) have the reference residue of
# every row, at the same codon numbers.
REPLACE = {
    ('klebsiella', 'gyrA_1_AF052258.1'): 'gyrA_2_NZ_KN046818.1',
    ('klebsiella', 'parC_1_AF303641.1'): 'parC_2_NZ_KN046818.1',
}
REFERENCES = os.path.join(HERE, 'pointfinder-references.fsa')

# AMRFinderPlus gene symbols where PointFinder's differ, for the CABBAGE
# look-up of a mutation (gyrA_S83L); a coding change on a promoter template
# is named for the gene (mtrR_G45D). PointFinder numbers residues as
# AMRFinderPlus does: at all 148 rows whose gene and position also have a
# CABBAGE symbol in the species (release 2026-07), the reference residue
# is the same, the enterococcal gyrA and parC included (their templates
# start at codon 2 and are numbered as E. coli's, as the literature's S83
# and E87 are). Klebsiella rpsL and Campylobacter porA get no symbols:
# their rows do not follow the protein's numbering.
AMR_GENE = {
    ('staphylococcus_aureus', 'grlA'): 'parC',
    ('staphylococcus_aureus', 'grlB'): 'parE',
    ('neisseria_gonorrhoeae', 'porB'): 'porB1b',
}
NO_SYMBOL = {('klebsiella', 'rpsL'), ('campylobacter', 'porA')}


def symbol_gene(folder, gene):
    if (folder, gene) in NO_SYMBOL:
        return None
    return AMR_GENE.get((folder, gene), re.sub(r'-promoter$', '', gene))

PROMOTER_RE = re.compile(r'(.+promoter)-size-(\d+)(?:bp)')
CODE = {}
_aa = 'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG'
_i = 0
for x in 'TCAG':
    for y in 'TCAG':
        for z in 'TCAG':
            CODE[x + y + z] = _aa[_i]
            _i += 1


def read_lines(path):
    if not os.path.exists(path):
        return []
    with open(path, newline='') as f:
        return [l.rstrip('\r\n') for l in f]


def read_fasta(path):
    seqs, name, buf = [], None, []
    for line in read_lines(path):
        line = line.strip()
        if not line:
            continue
        if line.startswith('>'):
            if name:
                seqs.append((name, ''.join(buf).upper()))
            name, buf = line[1:].split()[0], []
        else:
            buf.append(line)
    if name:
        seqs.append((name, ''.join(buf).upper()))
    return seqs


def tpos(t, pos):
    """0-based template index of a row position (nucleotide coordinates)."""
    p = t['promoter']
    return p + pos - 1 if pos > 0 else p + pos


def codon_at(t, n):
    i = t['promoter'] + 3 * (n - 1)
    return t['seq'][i:i + 3] if i >= 0 and i + 3 <= t['length'] else ''


def parse_requires(text, templates_by_name, templates_by_gene):
    """'gyrA_1_CP073768.1_S83L.W.A.V;pmrA_..._S39I,pmrA_..._R81S' ->
    [[{template, ref, pos, alts}], ...] (OR of AND-groups), and the parts
    that cannot be resolved to a template."""
    groups, unresolved = [], []
    for g in [x.strip() for x in text.split(';') if x.strip()]:
        members = []
        for m in [x.strip() for x in g.split(',') if x.strip()]:
            mm = re.match(r'^(.*)_([A-Za-z*]+)(-?\d+)([A-Za-z*]+(?:\.[A-Za-z*]+)*)$', m)
            if not mm:
                unresolved.append(m)
                members = None
                break
            name, ref, pos, alts = mm.group(1), mm.group(2).upper(), int(mm.group(3)), mm.group(4).upper().split('.')
            t = templates_by_name.get(name)
            if t is None and '_' not in name:
                # old style without the accession (gyrA_S83Y.F.A): the
                # species' template of that gene, when there is one
                cands = templates_by_gene.get(name, [])
                t = cands[0] if len(cands) == 1 else None
            if t is None:
                unresolved.append(m)
                members = None
                break
            members.append({'template': t['name'], 'ref': ref, 'pos': pos, 'alts': alts})
        if members:
            groups.append(members)
    return groups, unresolved


def read_amrfinder(path):
    """taxgroup -> the standard mutation symbols AMRFinderPlus curates."""
    out = {}
    for line in read_lines(path):
        if not line or line.startswith('#'):
            continue
        f = line.split('\t')
        if len(f) > 3:
            out.setdefault(f[0], set()).add(f[3])
    return out


def build_folder(db, folder, out, amrfinder):
    path = os.path.join(db, folder)
    genes = [g.strip() for g in read_lines(os.path.join(path, 'genes.txt')) if g.strip()]
    rna = {g.strip() for g in read_lines(os.path.join(path, 'RNA_genes.txt')) if g.strip()}
    promoters = {}
    screened = set()
    for g in genes:
        m = PROMOTER_RE.match(g)
        if m:
            promoters[m.group(1)] = int(m.group(2))
            screened.add(m.group(1))
        else:
            screened.add(g)

    replacements = dict(read_fasta(REFERENCES))
    templates, by_name, by_gene = [], {}, {}
    for name, seq in read_fasta(os.path.join(path, folder + '.fsa')):
        replaced = None
        if (folder, name) in REPLACE:
            replaced, name = name, REPLACE[(folder, name)]
            seq = replacements[name]
        parts = name.split('_')
        gene = parts[0]
        locus = gene.split('-')[0] if gene in rna else gene
        t = {
            'name': name, 'gene': gene, 'locus': locus,
            'accession': '_'.join(parts[2:]) if len(parts) > 2 else (parts[1] if len(parts) == 2 else ''),
            'length': len(seq), 'promoter': promoters.get(gene, 0), 'rna': gene in rna,
            'screened': gene in screened, 'seq': seq,
        }
        if not t['rna']:
            t['symbolGene'] = symbol_gene(folder, gene)
        if replaced:
            t['replaces'] = replaced
            by_name[replaced] = t  # the rows and requirements name the old template
        templates.append(t)
        by_name[name] = t
        by_gene.setdefault(gene, []).append(t)

    rows, dropped = [], []
    section = 'sub'
    for lineno, line in enumerate(read_lines(os.path.join(path, 'phenotypes.txt')), 1):
        if not line.strip():
            continue
        if line.startswith('#'):
            low = line.lower()
            if 'indel' in low:
                section = 'indel'
            elif low.startswith('#stop codons'):
                section = 'stop'
            elif section == 'stop':
                section = 'sub'
            continue
        f = line.split('\t') + [''] * 14
        ov = OVERRIDES.get((folder, lineno), {})
        raw = ov.get('template') or re.sub(r'-size-\d+bp', '', f[0].strip())
        t = by_name.get(raw)
        rec = {'line': lineno, 'template': raw, 'text': '\t'.join(x for x in f[:9] if x)[:200]}
        if t is None:
            dropped.append({**rec, 'reason': 'no such template in the combined FASTA'})
            continue
        try:
            pos = int(str(ov.get('pos', f[3].strip())))
        except ValueError:
            dropped.append({**rec, 'reason': f'position {f[3]!r} is not a number'})
            continue
        type_gene = f[1].strip().upper()
        ref_nuc, ref_codon, res = f[4].strip().upper(), f[5].strip(), f[6].strip()
        region = 'rna' if t['rna'] else ('promoter' if t['promoter'] and pos < 0 else 'coding')
        alts = [a.strip().upper() for a in res.split(',') if a.strip()]
        if section == 'stop':
            kind = 'stop_any'
        elif section == 'indel':
            rc = ref_codon.lower()
            nuc = region != 'coding' or type_gene != 'AA'
            if rc == 'ins':
                kind = 'ins' if nuc else 'aa_ins'
            elif rc == 'del':
                kind = 'del' if nuc else 'aa_del'
            else:
                dropped.append({**rec, 'reason': f'indel row with Ref_codon {ref_codon!r}'})
                continue
        else:
            kind = 'aa_sub' if region == 'coding' and type_gene == 'AA' else 'nuc_sub'

        # Check the row against its template.
        problem = None
        ref = ref_codon.upper()
        if kind == 'aa_sub':
            codon = codon_at(t, pos)
            aa = CODE.get(codon, '?')
            if not codon:
                problem = f'codon {pos} is outside the template'
            elif aa != ref:
                problem = f'row says {ref}{pos}, template codon {codon} is {aa}'
            if alts == [ref]:
                problem = problem or f'the alternative equals the reference ({ref}{pos}{ref})'
            elif ref in alts and not problem:
                # M. tuberculosis rpoB L440 lists M, L or P: L is no change
                alts = [a for a in alts if a != ref]
                ov = {**ov, 'why': ov.get('why') or f'{ref}{pos}{ref} left out: the alternative equals the reference'}
        elif kind == 'nuc_sub':
            i = tpos(t, pos)
            base = t['seq'][i] if 0 <= i < t['length'] else ''
            ref = ref.upper() or ref_nuc
            if not base:
                problem = f'position {pos} is outside the template'
            elif base != ref:
                problem = f'row says {ref} at {pos}, template has {base}'
        elif kind == 'del':
            i = tpos(t, pos)
            ok = any(t['seq'][i + d:i + d + len(a)] == a for a in alts for d in (-1, 0, 1)) if alts else False
            if not ok:
                problem = f'deleted sequence {",".join(alts)} is not at position {pos}'
        elif kind == 'ins':
            i = tpos(t, pos)
            if not 0 <= i < t['length']:
                problem = f'position {pos} is outside the template'
        elif kind in ('aa_ins', 'aa_del'):
            if not codon_at(t, pos):
                problem = f'codon {pos} is outside the template'
        if problem:
            dropped.append({**rec, 'reason': problem})
            continue

        pmids = [[p.strip() for p in a.split(';') if p.strip() and p.strip() != '-'] for a in f[9].split(',')]
        abx = [a.strip() for a in re.split(r',', f[8]) if a.strip() and a.strip().lower() not in ('unknown', 'none')]
        requires, unresolved = parse_requires(f[12].strip(), by_name, by_gene) if f[12].strip() else ([], [])
        row = {
            'line': lineno, 'template': t['name'], 'gene': t['gene'], 'kind': kind, 'region': region, 'pos': pos,
            'ref': ref if kind in ('aa_sub', 'nuc_sub') else '', 'alts': alts,
            'pmids': pmids, 'antibiotics': abx, 'class': f[7].strip(), 'mechanism': f[10].strip(),
            'notes': f[11].strip(), 'requires': requires, 'requiresText': f[12].strip(),
        }
        if unresolved:
            row['unresolved'] = unresolved
        if ov.get('why'):
            row['fixed'] = ov['why']
        if t.get('replaces'):
            row['fixed'] = row.get('fixed') or f"read on {t['name']} in place of {t['replaces']}, an E. coli sequence"

        if kind == 'aa_sub' and t.get('symbolGene'):
            row['symbols'] = [f"{t['symbolGene']}_{ref}{pos}{'Ter' if a == '*' else a}" for a in alts]
        if CAUTION_RE.search(row['notes']):
            row['caution'] = True
        curated = set().union(*[amrfinder.get(g, set()) for g in TAXGROUPS.get(folder, [])])
        if any(x in curated for x in row.get('symbols', [])):
            row['amrfinder'] = True
        rows.append(row)

    os.makedirs(out, exist_ok=True)
    # Index the combined FASTA (line endings normalised; kma_index would
    # strip the CRs anyway).
    with tempfile.TemporaryDirectory() as tmp:
        fa = os.path.join(tmp, folder + '.fsa')
        with open(fa, 'w') as f:
            for t in templates:
                f.write(f">{t['name']}\n")
                for i in range(0, len(t['seq']), 70):
                    f.write(t['seq'][i:i + 70] + '\n')
        subprocess.run([KMA_INDEX, '-i', fa, '-o', os.path.join(out, folder)], check=True,
                       stdout=subprocess.DEVNULL, stderr=subprocess.DEVNULL)
    comp = os.path.join(out, folder + '.comp.b')
    with open(comp, 'rb') as f:
        packed = gzip.compress(f.read(), compresslevel=9, mtime=0)
    with open(comp + '.gz', 'wb') as f:
        f.write(packed)
    os.remove(comp)
    cat = {
        'species': folder, 'label': LABELS.get(folder, folder),
        'templates': templates, 'rows': rows, 'dropped': dropped,
    }
    with open(os.path.join(out, 'catalogue.json.gz'), 'wb') as f:
        f.write(gzip.compress(json.dumps(cat, separators=(',', ':')).encode(), compresslevel=9, mtime=0))
    return cat


def main():
    ap = argparse.ArgumentParser()
    ap.add_argument('db')
    ap.add_argument('--amrfinder', required=True, help="AMRFinderPlus's AMRProt-mutation.tsv")
    ap.add_argument('--out', default=os.path.join(REPO, 'databases', 'pointfinder'))
    a = ap.parse_args()
    amrfinder = read_amrfinder(a.amrfinder)
    afp_version = read_lines(os.path.join(os.path.dirname(os.path.abspath(a.amrfinder)), 'version.txt'))
    afp_version = afp_version[0].strip() if afp_version else ''
    version = read_lines(os.path.join(a.db, 'VERSION'))[0].strip()
    try:
        commit = subprocess.run(['git', '-C', a.db, 'rev-parse', 'HEAD'], capture_output=True, text=True, check=True).stdout.strip()
        date = subprocess.run(['git', '-C', a.db, 'log', '-1', '--format=%cs'], capture_output=True, text=True, check=True).stdout.strip()
    except (subprocess.CalledProcessError, FileNotFoundError):
        commit, date = '', ''
    folders = []
    for line in read_lines(os.path.join(a.db, 'config')):
        if line.startswith('#') or not line.strip():
            continue
        folder = line.split('\t')[0].strip()
        if folder in SKIP_FOLDERS:
            continue
        folders.append(folder)
    if os.path.isdir(a.out):
        shutil.rmtree(a.out)
    manifest = {'database': 'PointFinder', 'version': version, 'commit': commit, 'date': date,
                'amrfinder_database': afp_version, 'species': []}
    for folder in folders:
        out = os.path.join(a.out, folder)
        cat = build_folder(a.db, folder, out, amrfinder)
        files = sorted(os.listdir(out))
        size = sum(os.path.getsize(os.path.join(out, x)) for x in files)
        manifest['species'].append({'folder': folder, 'label': cat['label'], 'files': files, 'bytes': size,
                                    'templates': len(cat['templates']), 'rows': len(cat['rows']), 'dropped': len(cat['dropped'])})
        caution = sum(1 for r in cat['rows'] if r.get('caution') and not r.get('amrfinder'))
        print(f"{folder}: {len(cat['templates'])} templates, {len(cat['rows'])} rows ({caution} under development, not in AMRFinderPlus), "
              f"{len(cat['dropped'])} dropped, {size:,} bytes")
        for d in cat['dropped']:
            print(f"    dropped line {d['line']} ({d['template']}): {d['reason']}")
    with open(os.path.join(a.out, 'manifest.json'), 'w') as f:
        json.dump(manifest, f, indent=1)


if __name__ == '__main__':
    sys.exit(main())
