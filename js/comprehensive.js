/* comprehensive.js — configurable pipeline: optional fastp QC/trim →
 * optional human read removal (Deacon) → KMA against any subset of
 * ResFinder + CARD + VFDB, and PointFinder's references for the species.
 *
 * Reuses the production runAnalysis() path verbatim for each database, so
 * the pipeline can never drift from what a single-database run does.
 */

import { runAnalysis } from './wasm-runtime.js';
import { ENGINE } from './engine.js';
import { removeHumanReads, HOST_INDEX, HOST_THRESHOLDS } from './host.js';
import { normaliseThresholds } from './thresholds.js';
import {
  pointfinderSpecies, loadPointFinderManifest, loadCatalogue, pointfinderDatabase, evaluatePointFinder,
} from './pointfinder.js';

export const COMPREHENSIVE_DBS = [
  { key: 'resfinder', label: 'Resistance genes (ResFinder)', kind: 'amr' },
  { key: 'card_homolog', label: 'Resistance homologs (CARD)', kind: 'amr' },
  { key: 'vfdb_core', label: 'Virulence factors (VFDB)', kind: 'vf' },
  // One index per species, chosen at run time from the species.
  { key: 'pointfinder', label: 'Point mutations (PointFinder)', kind: 'point' },
];

export function runFastp(files, { paired, nanopore, options } = {}) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL('./fastp-runner.worker.js?v=wasm32a', import.meta.url));
    worker.onmessage = (e) => {
      const m = e.data;
      if (m.type === 'log') {
        fastpLog?.(m.text);
      } else if (m.type === 'done') {
        worker.terminate();
        resolve(m);
      } else if (m.type === 'error') {
        worker.terminate();
        reject(new Error(m.message));
      }
    };
    worker.onerror = (err) => {
      worker.terminate();
      reject(new Error(err.message || 'fastp worker crashed'));
    };
    worker.postMessage({ type: 'run', files, paired, nanopore, options, engine: ENGINE });
  });
}

let fastpLog = null;
export function setFastpLog(fn) { fastpLog = fn; }

export async function runComprehensive(files, readType, {
  onStep, onProgress, log, runQc = true, removeHost = false, sampleName = 'sample', dbKeys, thresholds, fastpOptions,
  species = '',
} = {}) {
  const paired = readType === 'paired';
  const nanopore = readType === 'nanopore';
  // PointFinder runs only for a species it covers.
  const pointFolder = pointfinderSpecies(species);
  const dbs = (dbKeys?.length
    ? COMPREHENSIVE_DBS.filter(d => dbKeys.includes(d.key))
    : COMPREHENSIVE_DBS.filter(d => d.key !== 'pointfinder'))
    .filter(d => d.key !== 'pointfinder' || pointFolder);
  if (dbKeys?.includes('pointfinder') && !pointFolder) {
    log?.(species
      ? `PointFinder has no references for ${species}; point mutations were not checked.`
      : 'No species was given, so point mutations were not checked (PointFinder needs the species).', 'info');
  }
  const total = dbs.length + (runQc ? 1 : 0) + (removeHost ? 1 : 0);
  let stepNo = 0;
  const step = (text) => onStep?.(++stepNo, total, text);

  // ── Optional QC + trim ──
  let qc = null;
  let qcHtml = null;
  let qcReads = null;
  let analysisFiles = files;
  if (runQc) {
    step('Quality control (fastp)…');
    try {
      const out = await runFastp(files, { paired, nanopore, options: fastpOptions });
      qc = out.report;
      qcHtml = out.html || null;
      qcReads = out.reads || [];
      if (qcReads.length) analysisFiles = qcReads;
    } catch (err) {
      qc = null; // QC is best-effort: proceed with original reads
      fastpLog?.(`fastp failed (${err.message}); continuing with unfiltered reads`, 'warn');
    }
  }

  // ── Optional human read removal, on the trimmed reads when fastp ran ──
  // Best-effort like QC: a failure is reported in the results and the
  // databases run on the reads as they were.
  let host = null;
  if (removeHost) {
    const label = 'Removing human reads (Deacon)…';
    step(label);
    const inputs = analysisFiles;
    const names = inputs.length === 2
      ? [`${sampleName}_nohuman_1.fastq.gz`, `${sampleName}_nohuman_2.fastq.gz`]
      : [`${sampleName}_nohuman.fastq.gz`];
    // The equivalent native command, so the step can be reproduced.
    const outArgs = names.length === 2 ? `-o ${names[0]} -O ${names[1]}` : `-o ${names[0]}`;
    log?.(`$ deacon filter --deplete --abs-threshold ${HOST_THRESHOLDS.absThreshold} --rel-threshold ${HOST_THRESHOLDS.relThreshold} `
      + `${HOST_INDEX.path.split('/').pop()} ${inputs.map(f => f.name).join(' ')} ${outArgs}`, 'cmd');
    try {
      const out = await removeHumanReads(inputs, {
        names,
        onProgress: (done, all) => onProgress?.(`${label} ${Math.min(100, Math.round(100 * done / all))}%`),
        log,
      });
      const { readsIn, readsOut } = out.stats;
      log?.(`Removed ${(readsIn - readsOut).toLocaleString('en-US')} of ${readsIn.toLocaleString('en-US')} reads as human `
        + `(${fmtPct(readsIn ? (readsIn - readsOut) / readsIn : 0)}) in ${out.seconds.toFixed(1)}s. Index ${out.info}.`, 'ok');
      host = { ...out, inputNames: inputs.map(f => f.name) };
      analysisFiles = out.reads;
    } catch (err) {
      host = { error: err.message };
      log?.(`Human read removal failed (${err.message}). The databases run on all reads.`, 'warn');
    }
  }

  // ── Gene databases on the (optionally filtered) reads ──
  const config = {
    thresholds: normaliseThresholds(thresholds),
    nanopore,
  };
  const results = [];
  for (const db of dbs) {
    step(db.label + '…');
    if (db.key === 'pointfinder') {
      results.push(await runPointFinder(analysisFiles, pointFolder, config, log));
    } else {
      results.push(await runAnalysis(analysisFiles, db.key, config));
    }
  }
  onStep?.(total, total, 'Done');
  return { qc, qcHtml, qcReads, host, usedFiltered: analysisFiles !== files, results };
}

// KMA against the species' PointFinder references, then every curated
// position read from the result (pointfinder.js). A failure, a download
// that does not come (offline) included, fails this step only: the result
// carries the error, and the report and tables say so.
async function runPointFinder(files, folder, config, log) {
  let manifest, catalogue, r;
  try {
    [manifest, catalogue] = await Promise.all([loadPointFinderManifest(), loadCatalogue(folder)]);
    r = await runAnalysis(files, pointfinderDatabase(folder, manifest), config);
  } catch (err) {
    log?.(`PointFinder could not run: ${err.message}`, 'warn');
    return { database: 'pointfinder', dbLabel: 'PointFinder', exitCode: -1, files: {}, resTable: null, pointfinderError: err.message, pointfinderVersion: manifest?.version || '' };
  }
  r.pointfinderVersion = manifest.version;
  r.pointfinderCommit = manifest.commit;
  if (r.exitCode === 0) {
    try {
      r.pointfinder = evaluatePointFinder(r, catalogue, { nanopore: config.nanopore });
      const ev = r.pointfinder;
      log?.(`PointFinder (${ev.label}): ${ev.found.length} known mutation${ev.found.length === 1 ? '' : 's'} found; `
        + `${ev.read} of ${ev.screened} curated positions read.`, 'ok');
    } catch (err) {
      r.pointfinderError = err.message;
      log?.(`Reading the PointFinder positions failed: ${err.message}`, 'warn');
    }
  }
  return r;
}

// deacon_summary.json in the ZIP: the fields of `deacon filter --summary`
// that apply to this run, from the browser engine's counts.
export function hostSummary(host) {
  if (!host?.stats) return null;
  const { readsIn, readsOut, basesIn, basesOut } = host.stats;
  const [input, input2 = null] = host.inputNames;
  const [output, output2 = null] = host.reads.map(f => f.name);
  const ratio = (a, b) => (b ? a / b : 0);
  const kw = /k=(\d+), w=(\d+)/.exec(host.info || ''); // "k=31, w=61 (111,300,492 minimizers)"
  return {
    version: 'deacon 0.18.0 (deacon-wasm)',
    index: HOST_INDEX.path.split('/').pop(),
    input, input2, output, output2,
    k: kw ? Number(kw[1]) : null,
    w: kw ? Number(kw[2]) : null,
    abs_threshold: HOST_THRESHOLDS.absThreshold,
    rel_threshold: HOST_THRESHOLDS.relThreshold,
    prefix_length: 0,
    deplete: true,
    rename: false,
    ordered: true,
    seqs_in: readsIn,
    seqs_out: readsOut,
    seqs_out_proportion: ratio(readsOut, readsIn),
    seqs_removed: readsIn - readsOut,
    seqs_removed_proportion: ratio(readsIn - readsOut, readsIn),
    bp_in: basesIn,
    bp_out: basesOut,
    bp_out_proportion: ratio(basesOut, basesIn),
    bp_removed: basesIn - basesOut,
    bp_removed_proportion: ratio(basesIn - basesOut, basesIn),
    time: host.seconds,
  };
}

// ── QC summarisation ──

export function summariseQc(report) {
  if (!report || !report.summary) return null;
  const b = report.summary.before_filtering;
  const a = report.summary.after_filtering;
  const dup = report.duplication?.rate ?? null;
  const totalBefore = b.total_reads;
  const totalAfter = a.total_reads;
  const retained = totalBefore ? totalAfter / totalBefore : 0;
  return {
    rawReads: totalBefore,
    rawBases: b.total_bases,
    q30Before: b.q30_rate,
    q30After: a.q30_rate,
    gcBefore: b.gc_content,
    retained,
    adapterBases: report.adapter_cutting?.adapter_trimmed_reads != null
      ? report.adapter_cutting.adapter_trimmed_reads
      : null,
    duplication: dup,
  };
}

export function qcVerdict(m) {
  if (!m) return { tone: 'warn', text: 'Quality control could not run; results are from unfiltered reads.' };
  const bits = [];
  let tone = 'ok';
  if (m.q30After >= 0.85) bits.push(`${(100 * m.q30After).toFixed(1)}% of bases are Q30 after trimming (good)`);
  else if (m.q30After >= 0.7) { bits.push(`${(100 * m.q30After).toFixed(1)}% Q30 after trimming (moderate)`); tone = 'warn'; }
  else { bits.push(`only ${(100 * m.q30After).toFixed(1)}% Q30 after trimming (poor; interpret hits with caution)`); tone = 'bad'; }
  if (m.retained < 0.5) { bits.push(`only ${(100 * m.retained).toFixed(1)}% of reads were retained; check sample quality`); tone = 'bad'; }
  else if (m.adapterBases != null && m.rawReads && m.adapterBases / m.rawReads > 0.2) {
    bits.push(`${(100 * m.adapterBases / m.rawReads).toFixed(1)}% of reads had adapter sequence trimmed`);
  }
  if (m.duplication != null) bits.push(`duplication ${fmtPct(m.duplication)}`);
  return { tone, text: bits.join(' · ') };
}

// Percentage from a fraction; tiny non-zero rates show as "<0.1%" rather
// than a misleading "0.0%".
export function fmtPct(frac, digits = 1) {
  const p = 100 * (frac || 0);
  if (p > 0 && p < 0.1) return '<0.1%';
  return p.toFixed(digits) + '%';
}
