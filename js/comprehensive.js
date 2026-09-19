/* comprehensive.js — configurable pipeline: optional fastp QC/trim → KMA
 * against any subset of ResFinder + CARD + VFDB.
 *
 * Reuses the production runAnalysis() path verbatim for each database, so
 * the pipeline can never drift from what a single-database run does.
 */

import { runAnalysis } from './wasm-runtime.js';
import { ENGINE } from './engine.js';

export const COMPREHENSIVE_DBS = [
  { key: 'resfinder', label: 'Resistance genes (ResFinder)', kind: 'amr' },
  { key: 'card_homolog', label: 'Resistance homologs (CARD)', kind: 'amr' },
  { key: 'vfdb_core', label: 'Virulence factors (VFDB)', kind: 'vf' },
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

export async function runComprehensive(files, readType, { onStep, runQc = true, dbKeys, thresholds, fastpOptions } = {}) {
  const paired = readType === 'paired';
  const nanopore = readType === 'nanopore';
  const dbs = dbKeys?.length
    ? COMPREHENSIVE_DBS.filter(d => dbKeys.includes(d.key))
    : COMPREHENSIVE_DBS;
  const total = dbs.length + (runQc ? 1 : 0);
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

  // ── Gene databases on the (optionally filtered) reads ──
  const config = {
    id_threshold: thresholds?.id_threshold ?? 0.90,
    mrc: thresholds?.mrc ?? 0.60,
    nanopore,
  };
  const results = [];
  for (const db of dbs) {
    step(db.label + '…');
    results.push(await runAnalysis(analysisFiles, db.key, config));
  }
  onStep?.(total, total, 'Done');
  return { qc, qcHtml, qcReads, usedFiltered: analysisFiles !== files, results };
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
