/* hamronize.js — the run's gene hits in the PHA4GE hAMRonization format
 * (hAMRonization 1.3.1, https://github.com/pha4ge/hAMRonization), so they
 * can be merged with other tools' results by `hamronize summarize` or
 * loaded into anything else that reads the format.
 *
 * One row per called hit and database (the hits that pass both sliders;
 * partial genes and hits below the thresholds stay in the .res files),
 * exactly as KMA reported it, with no merging of ResFinder and CARD calls
 * as in the report. Columns, their order, the
 * empty cells and number formatting follow the Python package's TSV and
 * JSON writers. Where the spec leaves a choice, the package's own parsers
 * set the precedent:
 *   - KMA .res columns as in its kmerresistance parser: Template_Coverage →
 *     coverage_percentage, Query_Identity → sequence_identity (identity over
 *     the aligned part, the same number the gene tables and the sliders
 *     use), Depth → coverage_depth, Template_length → reference_gene_length.
 *   - Drug classes, antibiotics and phenotype as in its ResFinder parser,
 *     here from ResFinder's phenotype table (the same annotation as the
 *     report): classes and antibiotics lower-cased, sorted, ", "-joined;
 *     notes and PMIDs condensed into predicted_phenotype_confidence_level.
 */

import { describeHit } from './genes.js';
import { parseResFile } from './resfinder-db.js';
import { classifyResult } from './hits.js';

// hAMRonizedResult's fields in declaration order: the TSV header.
export const HAMR_FIELDS = [
  // mandatory
  'input_file_name', 'gene_symbol', 'gene_name', 'reference_database_name', 'reference_database_version',
  'reference_accession', 'analysis_software_name', 'analysis_software_version', 'genetic_variation_type',
  // optional
  'antimicrobial_agent', 'coverage_percentage', 'coverage_depth', 'coverage_ratio', 'drug_class',
  'input_gene_length', 'input_gene_start', 'input_gene_stop', 'input_protein_length', 'input_protein_start',
  'input_protein_stop', 'input_sequence_id', 'nucleotide_mutation', 'nucleotide_mutation_interpretation',
  'predicted_phenotype', 'predicted_phenotype_confidence_level', 'amino_acid_mutation',
  'amino_acid_mutation_interpretation', 'reference_gene_length', 'reference_gene_start', 'reference_gene_stop',
  'reference_protein_length', 'reference_protein_start', 'reference_protein_stop', 'resistance_mechanism',
  'strand_orientation', 'sequence_identity',
];

// hAMRonization.constants.GENE_PRESENCE: KMA reports genes, not mutations.
const GENE_PRESENCE = 'gene_presence_detected';

// Name and version of each index as the spec wants them. VFDB has no
// version number; the spec then asks for the release date, taken from the
// gzip header of VFDB_setA_nt.fas.gz in the index record (Zenodo 22102687).
const DATABASES = {
  resfinder: { name: 'ResFinder', version: '2.6.0' },
  card_homolog: { name: 'CARD', version: '4.0.1' },
  vfdb_core: { name: 'VFDB', version: '2026-05-01' },
};

// The accession inside a template name:
//   ResFinder  mecA_1_NC_002745                      → NC_002745
//   CARD       gb|KC243783.1|+|0-2007|ARO:3000617|…  → ARO:3000617
//   VFDB       VFG001276(gb|WP_024937002) (lukF-PV)… → VFG001276
// ResFinder names are gene_variant_accession, and genes and accessions can
// both contain underscores (VanHAX_PT_1_DQ018710); the few with no
// accession ("…_NoAcc") keep the whole name.
function accession(template, db, hit) {
  if (db === 'card_homolog') return hit.aro ? `ARO:${hit.aro}` : template;
  if (db === 'vfdb_core') return hit.vf?.id || template.split(/[\s(]/)[0];
  const m = String(template).match(/^.+?_\d+_(.+)$/);
  return m && m[1] !== 'NoAcc' ? m[1] : template;
}

// ResFinderIO._condense_notes: the notes, then "PMIDs: …", joined by ". ".
// The phenotype table keeps a gene's notes in one ";"-separated cell.
function condenseNotes(entry) {
  const lines = String(entry.notes || '').split(';').map(s => s.trim()).filter(Boolean);
  const pmids = [...new Set(String(entry.pmid || '').split(/[,;\s]+/).filter(Boolean))].sort();
  if (pmids.length) lines.push('PMIDs: ' + pmids.join(', '));
  return lines.length ? lines.join('. ') : null;
}

const joinSorted = (items) => {
  const out = [...new Set(items.map(s => String(s).trim().toLowerCase()).filter(Boolean))].sort();
  return out.length ? out.join(', ') : null;
};

const toFloat = (v) => { const x = parseFloat(v); return Number.isFinite(x) ? x : null; };
const toInt = (v) => { const x = parseInt(v, 10); return Number.isFinite(x) ? x : null; };

// hAMRonizedResult.__post_init__: the basename, without these suffixes
// removed in this order (so "x.fa.gz" becomes "x").
function inputFileName(name) {
  let n = String(name).split(/[\\/]/).pop();
  for (const suffix of ['.gz', '.fna', '.fasta', '.fsa', '.faa', '.fa']) {
    if (n.endsWith(suffix)) n = n.slice(0, -suffix.length);
  }
  return n;
}

// One hAMRonized record per called hit, for every database that ran.
//   results     runComprehensive() results (KMA output per database)
//   sample      the sample name (input_file_name)
//   version     openpathogen's version, for analysis_software_version
//   thresholds  the sliders; without them every .res row is exported
export function hamronizeResults(results, { sample, version, thresholds }) {
  const records = [];
  sample = inputFileName(sample);
  for (const r of results) {
    const db = DATABASES[r.database];
    if (!db || r.exitCode !== 0 || !r.resTable) continue;
    const rows = thresholds ? classifyResult(r, thresholds).called.map(h => h.row) : parseResFile(r.resTable).rows;
    for (const row of rows) {
      const hit = describeHit(row.Template, r.database, row);
      const amr = r.database !== 'vfdb_core';
      const entry = amr ? hit.entry : null;
      records.push({
        input_file_name: sample,
        gene_symbol: hit.name,
        // VFDB has a product name; ResFinder and CARD names are the symbol.
        gene_name: hit.vf?.product || hit.name,
        reference_database_name: db.name,
        reference_database_version: db.version,
        reference_accession: accession(row.Template, r.database, hit),
        analysis_software_name: 'openpathogen',
        analysis_software_version: version,
        genetic_variation_type: GENE_PRESENCE,
        antimicrobial_agent: amr ? joinSorted(hit.antibiotics) : null,
        coverage_percentage: toFloat(row.Template_Coverage),
        coverage_depth: toFloat(row.Depth),
        drug_class: amr ? joinSorted(hit.classes) : null,
        predicted_phenotype: amr && hit.classes.length ? 'antimicrobial resistance' : null,
        predicted_phenotype_confidence_level: entry ? condenseNotes(entry) : null,
        reference_gene_length: toInt(row.Template_length),
        resistance_mechanism: entry?.mechanism?.trim() || null,
        sequence_identity: toFloat(row.Query_Identity),
      });
    }
  }
  return records;
}

// Integer columns print as integers, as the Python dataclass casts them.
const INT_FIELDS = new Set(['input_gene_length', 'input_gene_start', 'input_gene_stop', 'input_protein_length',
  'input_protein_start', 'input_protein_stop', 'reference_gene_length', 'reference_gene_start', 'reference_gene_stop',
  'reference_protein_length', 'reference_protein_start', 'reference_protein_stop']);

// A value as Python's str() prints it. Floats keep ".0" when integral
// (100.0); otherwise both languages print the shortest round-trip form.
function value(record, field) {
  const v = record[field];
  if (v == null) return '';
  if (typeof v !== 'number' || INT_FIELDS.has(field)) return String(v);
  return Number.isInteger(v) ? v.toFixed(1) : String(v);
}

// csv.DictWriter(delimiter='\t'), QUOTE_MINIMAL: quote a value holding a
// tab, a double quote or a line break.
function tsvCell(s) {
  return /[\t"\r\n]/.test(s) ? `"${s.replace(/"/g, '""')}"` : s;
}

// The TSV the Python writer produces: a header row, then one row per hit.
// With no hits the writer emits nothing at all; the header alone is kept
// here so the file still says what it is (summarize treats it as empty).
export function toTsv(records) {
  const lines = [HAMR_FIELDS.join('\t')];
  for (const rec of records) lines.push(HAMR_FIELDS.map(f => tsvCell(value(rec, f))).join('\t'));
  return lines.join('\n') + '\n';
}

// json.dumps of each record (ensure_ascii, ", " and ": " separators), every
// value a string, in a list: the Python JSON writer. It writes "" for any
// falsy value, so a 0 becomes "" there too.
export function toJson(records) {
  const dumps = (s) => JSON.stringify(s).replace(/[\u0080-￿]/g, (c) => '\\u' + c.charCodeAt(0).toString(16).padStart(4, '0'));
  const entries = records.map(rec =>
    '{' + HAMR_FIELDS.map(f => `${dumps(f)}: ${dumps(rec[f] === 0 ? '' : value(rec, f))}`).join(', ') + '}');
  return '[' + entries.join(', ') + ']\n';
}
