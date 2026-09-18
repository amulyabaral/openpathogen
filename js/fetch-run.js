/* fetch-run.js — fetch a public sequencing run by accession from ENA
 * (which mirrors all NCBI SRA runs), and detect what it is before any
 * bytes move.
 *
 * The ENA Portal API supplies organism, instrument platform, library
 * layout and per-file sizes, so the app can set the read type itself and
 * show the user what was found. The FASTQs then download over ENA's
 * CORS-enabled HTTPS endpoints — the same channel the "Full dataset"
 * example already uses.
 *
 * Files above MAX_FILE_BYTES are refused up front: they are staged fully in
 * browser memory (and optionally cached in IndexedDB), which is impractical
 * beyond ~1 GB per file.
 */

const API = 'https://www.ebi.ac.uk/ena/portal/api/filereport';
const RUN_ACCESSION = /^[SED]RR\d+$/i;

export const MAX_FILE_BYTES = 1024 * 1024 * 1024; // 1 GB per file
export const MAX_FILE_LABEL = '1 GB';

export function isRunAccession(value) {
  return RUN_ACCESSION.test(String(value || '').trim());
}

export async function lookupRun(accession) {
  const acc = String(accession || '').trim();
  if (!isRunAccession(acc)) {
    throw new Error(`"${acc}" is not a run accession. Enter a single run accession (SRR…, ERR… or DRR…).`);
  }
  const fields = 'run_accession,fastq_ftp,fastq_bytes,library_layout,instrument_platform,instrument_model,scientific_name';
  const url = `${API}?accession=${encodeURIComponent(acc)}&result=read_run&fields=${fields}&format=json`;
  const resp = await fetch(url);
  if (!resp.ok) throw new Error(`ENA lookup failed (HTTP ${resp.status}).`);
  const rows = await resp.json().catch(() => null);
  if (!Array.isArray(rows) || rows.length === 0) {
    throw new Error(`No public run found for ${acc} in ENA. Runs submitted only to NCBI SRA may not be available.`);
  }
  if (rows.length > 1) {
    throw new Error(`${acc} resolves to ${rows.length} runs. Enter a single run accession (SRR…, ERR… or DRR…).`);
  }
  const r = rows[0];
  if (!r.fastq_ftp) {
    throw new Error(`No FASTQ files are available for ${acc} in ENA. This run may be stored only in another format.`);
  }

  const urls = String(r.fastq_ftp).split(';').map(p => 'https://' + p.replace(/^ftp:/, ''));
  const sizes = String(r.fastq_bytes || '').split(';').map(Number);
  const files = urls.map((url, i) => ({ url, name: url.split('/').pop(), bytes: sizes[i] || 0 }));

  const tooBig = files.filter(f => f.bytes > MAX_FILE_BYTES);
  if (tooBig.length) {
    throw new Error(`File${tooBig.length > 1 ? 's' : ''} in this run exceed the ${MAX_FILE_LABEL}-per-file limit (${tooBig.map(f => f.name).join(', ')}). Download the run yourself and upload the files, or choose a smaller run.`);
  }

  const platform = String(r.instrument_platform || '').toUpperCase();
  const layout = String(r.library_layout || '').toUpperCase() === 'PAIRED' ? 'paired' : 'single';
  return {
    accession: r.run_accession || acc,
    organism: r.scientific_name || 'unknown organism',
    platform,
    instrument: r.instrument_model || '',
    layout,
    readType: /NANOPORE/i.test(platform) ? 'nanopore' : layout,
    files,
    totalBytes: files.reduce((n, f) => n + (f.bytes || 0), 0),
  };
}
