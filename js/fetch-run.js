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
 *
 * ENA can list three files for a paired run: the pair as _1 and _2, and a
 * file of reads without a mate (DRR315548.fastq.gz next to DRR315548_1 and
 * _2). The pair is matched by name, never by position, and the odd file is
 * reported and left out.
 */

const API = 'https://www.ebi.ac.uk/ena/portal/api/filereport';
const RUN_ACCESSION = /^[SED]RR\d+$/i;

export const MAX_FILE_BYTES = 1024 * 1024 * 1024; // 1 GB per file
export const MAX_FILE_LABEL = '1 GB';

export function isRunAccession(value) {
  return RUN_ACCESSION.test(String(value || '').trim());
}

// 1 or 2 for an ENA file named as one half of a pair (…_1.fastq.gz), else 0.
export function mateOf(name) {
  const m = /_([12])\.(?:fastq|fq)(?:\.gz)?$/i.exec(String(name || ''));
  return m ? Number(m[1]) : 0;
}

// Which of a run's listed files to use, and in which slot.
//   files   [{name, …}] as ENA lists them
//   layout  'paired' when the run is paired-end and both mates are wanted
// Returns { files: [{…, slot}], skipped: [{…}], layout } where layout is
// what will actually run: 'paired' only when an R1/R2 pair was found.
export function chooseRunFiles(files, layout) {
  const r1 = files.find(f => mateOf(f.name) === 1);
  const r2 = files.find(f => mateOf(f.name) === 2);
  if (layout === 'paired' && r1 && r2) {
    return {
      files: [{ ...r1, slot: 'r1' }, { ...r2, slot: 'r2' }],
      skipped: files.filter(f => f !== r1 && f !== r2),
      layout: 'paired',
    };
  }
  const first = files.find(f => mateOf(f.name) !== 2) || files[0];
  return { files: first ? [{ ...first, slot: 'r1' }] : [], skipped: files.filter(f => f !== first), layout: 'single' };
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
  const listed = urls.map((url, i) => ({ url, name: url.split('/').pop(), bytes: sizes[i] || 0 }));

  const platform = String(r.instrument_platform || '').toUpperCase();
  const nanopore = /NANOPORE/i.test(platform);
  const listedLayout = String(r.library_layout || '').toUpperCase() === 'PAIRED' ? 'paired' : 'single';
  // Nanopore reads have no mates whatever the record says.
  const { files, skipped, layout } = chooseRunFiles(listed, nanopore ? 'single' : listedLayout);

  const tooBig = files.filter(f => f.bytes > MAX_FILE_BYTES);
  if (tooBig.length) {
    throw new Error(`File${tooBig.length > 1 ? 's' : ''} in this run exceed the ${MAX_FILE_LABEL}-per-file limit (${tooBig.map(f => f.name).join(', ')}). Download the run yourself and upload the files, or choose a smaller run.`);
  }

  return {
    accession: r.run_accession || acc,
    organism: r.scientific_name || 'unknown organism',
    platform,
    instrument: r.instrument_model || '',
    layout,
    listedLayout,
    readType: nanopore ? 'nanopore' : layout,
    files,
    skipped,
    totalBytes: files.reduce((n, f) => n + (f.bytes || 0), 0),
  };
}
