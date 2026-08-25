/* comprehensive.js — one-click pipeline: fastp QC/trim → KMA against
 * ResFinder + CARD + VFDB.
 *
 * Reuses the production runAnalysis() path verbatim for each database, so
 * the simple mode can never drift from what advanced mode runs.
 */

import { runAnalysis } from './wasm-runtime.js';

export const COMPREHENSIVE_DBS = [
  { key: 'resfinder', label: 'Resistance genes (ResFinder)', kind: 'amr' },
  { key: 'card_homolog', label: 'Resistance homologs (CARD)', kind: 'amr' },
  { key: 'vfdb_core', label: 'Virulence factors (VFDB)', kind: 'vf' },
];

export function runFastp(files, { paired, nanopore } = {}) {
  return new Promise((resolve, reject) => {
    const worker = new Worker(new URL('./fastp-runner.worker.js', import.meta.url));
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
    worker.postMessage({ type: 'run', files, paired, nanopore });
  });
}

let fastpLog = null;
export function setFastpLog(fn) { fastpLog = fn; }

export async function runComprehensive(files, readType, { onStep } = {}) {
  const paired = readType === 'paired';
  const nanopore = readType === 'nanopore';
  const step = (i, n, text) => onStep?.(i, n, text);

  // ── Step 1: QC + trim ──
  step(1, 4, 'Quality control (fastp)…');
  let qc = null;
  let analysisFiles = files;
  try {
    const out = await runFastp(files, { paired, nanopore });
    qc = out.report;
    if (out.reads && out.reads.length) analysisFiles = out.reads;
  } catch (err) {
    qc = null; // QC is best-effort: proceed with original reads
    fastpLog?.(`fastp failed (${err.message}); continuing with unfiltered reads`, 'warn');
  }

  // ── Steps 2–4: gene databases on the filtered reads ──
  const config = { id_threshold: 0.90, mrc: 0.60, nanopore };
  const results = [];
  for (let i = 0; i < COMPREHENSIVE_DBS.length; i++) {
    const db = COMPREHENSIVE_DBS[i];
    step(2 + i, 4, db.label + '…');
    results.push(await runAnalysis(analysisFiles, db.key, config));
  }
  step(4, 4, 'Done');
  return { qc, results, usedFiltered: analysisFiles !== files };
}

// ── QC summarisation ──

export function summariseQc(report, paired) {
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
    paired,
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
  if (m.duplication != null) bits.push(`duplication about ${(100 * m.duplication).toFixed(1)}%`);
  return { tone, text: bits.join(' · ') };
}
