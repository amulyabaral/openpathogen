/* education.js — the Learn page.
 *
 * A guided case-study lesson that runs the REAL production pipeline
 * (fastp → KMA against ResFinder + CARD + VFDB, via comprehensive.js)
 * inside a step-by-step narrative. Checkpoint quizzes are graded two ways:
 * a fixed correct answer (the biology of this dataset) plus a live
 * verification line computed from the student's own run output, so every
 * answer is backed by the results on their screen.
 *
 * Instructor view (toggle in the toolbar, or education.html?instructor=1)
 * reveals answer keys, expected results and per-step teaching notes.
 */

import { initWasm, setLogCallback, formatBytes } from './wasm-runtime.js';
import { loadPhenotypesDB, parseResFile } from './resfinder-db.js';
import { cacheDBFile, getCachedDBFile } from './db.js';
import { fetchAssetWithProgress } from './assets.js';
import { runComprehensive, summariseQc, qcVerdict, setFastpLog } from './comprehensive.js';

// ── State ──

const state = {
  files: { r1: null, r2: null },
  dataLoaded: false,
  loadingData: false,
  running: false,
  ran: false,
  runError: null,
  qc: null,
  rows: { resfinder: [], card_homolog: [], vfdb_core: [] },
  runFiles: [],           // [{ database, files }] for the ZIP download
  answers: {},            // quiz id -> { wrongPicks: [i], correct: bool, done: bool }
  instructor: false,
};

const $ = (sel) => document.querySelector(sel);

// ── Terminal (same pattern as the main app) ──

const term = {
  el: null,
  push(msg, level) {
    if (!this.el) this.el = document.getElementById('terminal');
    const line = document.createElement('div');
    line.className = 'term-line term-' + (level || 'info');
    line.textContent = msg;
    this.el.appendChild(line);
    this.el.scrollTop = this.el.scrollHeight;
  },
  text() { return [...this.el.querySelectorAll('.term-line')].map(l => l.textContent).join('\n') + '\n'; },
};

// ── The dataset (same example run the main app loads, straight from ENA) ──

const DATASET = {
  name: 'S. aureus USA300_TCH1516',
  accession: 'SRR10341524',
  files: [
    {
      slot: 'r1', name: 'SRR10341524_1.fastq.gz', bytes: 22326358,
      url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/024/SRR10341524/SRR10341524_1.fastq.gz',
    },
    {
      slot: 'r2', name: 'SRR10341524_2.fastq.gz', bytes: 26024204,
      url: 'https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR103/024/SRR10341524/SRR10341524_2.fastq.gz',
    },
  ],
};

// Genes this lesson highlights in result tables (also used for the live checks).
const KEY_GENES = [
  { re: /^mecA(_|\b)/, label: 'mecA' },
  { re: /^blaZ(_|\b)/, label: 'blaZ' },
  { re: /\(lukS-PV\)/, label: 'lukS-PV' },
  { re: /\(lukF-PV\)/, label: 'lukF-PV' },
  { re: /\(hla\)/, label: 'hla' },
];

// ── Live-verification helpers (run against the student's own results) ──

function findRows(dbKey, re) {
  return (state.rows[dbKey] || []).filter(r => re.test(r.Template || ''));
}

function verifyLine(rows) {
  return rows.map(r =>
    `${r.Template} — ${r.Template_Identity}% identity · ${r.Template_Coverage}% coverage · ${r.Depth}× depth`
  ).join('\n');
}

// ── Lesson content ──

const STEPS = [
  {
    id: 'case', type: 'read', title: 'A clinical case',
    html: `
      <p>A patient in the United States has a severe skin and soft-tissue infection. A
      <em>Staphylococcus aureus</em> isolate from the wound has been sequenced.</p>
      <p><em>S. aureus</em> is a common bacterium. Some strains carry genes that make them
      resistant to certain antibiotics, or more likely to cause severe disease. Your task is to
      examine the sequencing data and answer one question: <strong>is this isolate
      methicillin-resistant <em>S. aureus</em> (MRSA)?</strong></p>
      <div class="edu-defs">
        <div class="edu-def"><b>Reads (FASTQ)</b><span>Short DNA fragments from the sample. Each fragment carries a quality score.</span></div>
        <div class="edu-def"><b>Resistance gene</b><span>A gene that makes a bacterium resistant to an antibiotic, for example by destroying the drug or pumping it out of the cell.</span></div>
        <div class="edu-def"><b>MRSA</b><span><em>S. aureus</em> that carries <code>mecA</code> (or <code>mecC</code>) and is therefore resistant to nearly all beta-lactam antibiotics.</span></div>
        <div class="edu-def"><b>Virulence factor</b><span>A gene product that helps a bacterium cause disease, such as a toxin or a protein that evades the immune system.</span></div>
      </div>`,
    note: `Timing: about five minutes as a walkthrough. A common misconception to address early:
      resistance genes are carried by the bacterium (often on mobile plasmids or transposons),
      not by the patient. The isolate is USA300_TCH1516, the reference strain of the community-
      associated MRSA clone USA300 (sequence type 8), sequenced on a MiSeq (ENA SRR10341524).
      Students analyse real public data.`,
  },
  {
    id: 'how', type: 'read', title: 'How the analysis works',
    html: `
      <p>Laboratories screen sequencing data for resistance and virulence genes in the following steps:</p>
      <div class="edu-flow">
        <span class="edu-flow-node">Sequencing reads</span>
        <span class="edu-flow-arrow">→</span>
        <span class="edu-flow-node">Quality control<small>adapter and quality trimming</small></span>
        <span class="edu-flow-arrow">→</span>
        <span class="edu-flow-node">Alignment<small>reads compared with gene databases</small></span>
        <span class="edu-flow-arrow">→</span>
        <span class="edu-flow-node">Gene list<small>resistance and virulence calls</small></span>
      </div>
      <p>Each read is compared with three curated reference databases:</p>
      <div class="edu-defs">
        <div class="edu-def"><b>ResFinder</b><span>Acquired resistance genes. The classic database used for clinical screening.</span></div>
        <div class="edu-def"><b>CARD</b><span>A comprehensive resistance database. When two databases report the same gene, the finding is more reliable.</span></div>
        <div class="edu-def"><b>VFDB Set A</b><span>Virulence factors that have been verified experimentally.</span></div>
      </div>
      <p>The analysis runs on your own computer. The sequencing files are not uploaded to any
      server.</p>`,
    note: `If students ask why two resistance databases are used: agreement between databases is
      standard practice. Benchmarks show that tools and databases can disagree substantially on
      the genes they report, so finding a gene in both ResFinder and CARD increases confidence.
      This is also the argument the openpathogen paper makes for running several databases in one
      pipeline.`,
  },
  {
    id: 'data', type: 'load', title: 'Get the sequencing data',
    html: `
      <p>Download the reads for this lesson. The download is about 48&nbsp;MB, straight from
      ENA. The data come from a published <em>S. aureus</em> isolate (accession
      <code>${DATASET.accession}</code>).</p>`,
    note: `Classroom bandwidth: 48 MB per machine, downloaded once from ENA — the files are cached
      in the browser afterwards. If your network is slow, open the lesson once on each machine
      before class. Any other public run can be fetched the same way on the main page.`,
  },
  {
    id: 'run', type: 'run', title: 'Run the analysis',
    html: `
      <p>The analysis first checks read quality, then screens the reads against the three
      databases.</p>
      <p>The analysis may take a few minutes. The first run also downloads the reference
      databases. The Logs panel at the bottom of the page shows the output of each tool.</p>`,
    note: `On a typical laptop the example run takes one to two minutes once the databases are
      cached. The run needs about 1–2 GB of free browser memory; closing heavy tabs helps if it fails.
      The Logs panel shows the exact command line, which is useful for teaching reproducibility:
      the same command and databases always give the same answer.`,
  },
  {
    id: 'read-table', type: 'read', title: 'Reading the results',
    html: `
      <p>Each detected gene is reported with four values:</p>
      <div class="edu-defs">
        <div class="edu-def"><b>Identity</b><span>How similar your reads are to the reference gene. A value of 100% means an exact match.</span></div>
        <div class="edu-def"><b>Coverage</b><span>The share of the reference gene that your reads recover. A value of 100% means the whole gene was found.</span></div>
        <div class="edu-def"><b>Depth</b><span>How many reads align to the gene on average. High depth means strong evidence.</span></div>
        <div class="edu-def"><b>p-value</b><span>The probability that the match arose by chance. A very small value means the match is real.</span></div>
      </div>
      <p>A gene call is reliable when identity, coverage and depth are all high. Genes referred
      to in the checkpoints below are highlighted in the result tables.</p>`,
    note: `Expected results for this run: ResFinder reports six genes — mecA and blaZ (both
      beta-lactam resistance, ~100% identity and coverage), the aminoglycoside genes ant(6)-Ia and
      aph(3')-III, and the macrolide genes msr(A) and mph(C). CARD concords on mecA, blaZ, msrA,
      mphC and APH(3')-IIIa (plus efflux and regulatory homologs, depending on thresholds).
      VFDB reports dozens of genes, including lukS-PV and lukF-PV (PVL), hla, hlg and the ica
      operon. Exact depths vary slightly between machines; this is expected and worth discussing.`,
  },
  {
    id: 'q1', type: 'quiz', title: 'Checkpoint: the resistance gene',
    question: 'Treatment with a beta-lactam antibiotic failed. Look at the "Resistance genes (ResFinder)" table. Which gene makes this isolate MRSA?',
    options: [
      { label: 'mecA', why: 'Correct. mecA encodes an altered penicillin-binding protein (PBP2a) that binds beta-lactam antibiotics poorly. It is the mechanism that defines MRSA and is carried on a mobile genetic element called SCCmec.' },
      { label: 'blaZ', why: 'blaZ encodes a beta-lactamase. It breaks down penicillin, but not beta-lactamase-stable antibiotics such as flucloxacillin. It explains penicillin resistance, not methicillin resistance.' },
      { label: 'vanA', why: 'vanA causes resistance to vancomycin. It is not present in this isolate.' },
      { label: 'tetM', why: 'tetM causes resistance to tetracycline. It is not present in this isolate.' },
    ],
    correct: 0,
    verify: () => {
      const rf = verifyLine(findRows('resfinder', /^mecA(_|\b)/));
      const card = verifyLine(findRows('card_homolog', /^mecA(_|\b)|\(mecA\)/));
      const parts = [];
      if (rf) parts.push('ResFinder: ' + rf);
      if (card) parts.push('CARD: ' + card);
      return parts.length ? parts.join('\n') : null;
    },
    note: `Distractor logic: blaZ is the strongest distractor because it does appear in the table.
      Students must distinguish the penicillinase (blaZ) from the MRSA determinant (mecA). A useful
      question to ask the class: which beta-lactam antibiotics does each gene defeat?`,
  },
  {
    id: 'q2', type: 'quiz', title: 'Checkpoint: the penicillinase',
    question: 'ResFinder also detected a second resistance gene: the staphylococcal penicillinase. Which gene is it?',
    options: [
      { label: 'blaZ', why: 'Correct. blaZ encodes a beta-lactamase that breaks down penicillin. It is present in most S. aureus lineages, which is one reason penicillin is generally not effective against S. aureus infections.' },
      { label: 'blaKPC-3', why: 'blaKPC-3 is a carbapenemase found mainly in Gram-negative bacteria such as Klebsiella pneumoniae.' },
      { label: 'blaTEM-1', why: 'blaTEM-1 is a beta-lactamase common in Enterobacterales. It is not the gene reported here.' },
      { label: 'ermC', why: 'ermC causes resistance to macrolide antibiotics, which belong to a different drug class.' },
    ],
    correct: 0,
    verify: () => {
      const rows = findRows('resfinder', /^blaZ(_|\b)/);
      return rows.length ? verifyLine(rows) : null;
    },
    note: `Teaching point: two beta-lactam resistance genes with two different mechanisms, drug
      destruction versus target alteration. This maps onto the classic therapeutic staircase of
      penicillin, then flucloxacillin, then failure because of mecA.`,
  },
  {
    id: 'q3', type: 'quiz', title: 'Checkpoint: the virulence factor',
    question: 'Look at the "Virulence factors (VFDB)" table. Does this isolate carry the genes for Panton–Valentine leukocidin (PVL)?',
    options: [
      { label: 'Yes, both lukS-PV and lukF-PV are present', why: 'Correct. PVL is a two-component toxin formed by LukS-PV and LukF-PV. It damages white blood cells and is associated with severe skin infections and pneumonia. This isolate belongs to a community-associated MRSA lineage, which typically carries PVL.' },
      { label: 'No, neither gene is detected', why: 'Look at the VFDB table again and search for lukS-PV and lukF-PV.' },
      { label: 'Only lukS-PV', why: 'Look closely: both components of the toxin are reported.' },
      { label: 'Only lukF-PV', why: 'Look closely: both components of the toxin are reported.' },
    ],
    correct: 0,
    verify: () => {
      const s = findRows('vfdb_core', /\(lukS-PV\)/);
      const f = findRows('vfdb_core', /\(lukF-PV\)/);
      const parts = [];
      if (s.length) parts.push(verifyLine(s));
      if (f.length) parts.push(verifyLine(f));
      return parts.length ? parts.join('\n') : null;
    },
    note: `PVL is what makes this a community-associated MRSA story rather than a hospital MRSA
      story. Discussion prompt: the virulence genes explain why this strain causes severe skin
      disease, while mecA explains why first-line therapy failed. The clinical picture needs both.`,
  },
  {
    id: 'q4', type: 'quiz', title: 'Checkpoint: interpretation',
    question: 'The treating team asks how the result should be reported. Which statement is correct?',
    options: [
      { label: 'The genotype strongly supports MRSA, but phenotypic susceptibility testing should confirm the result before treatment decisions', why: 'Correct. Detecting mecA at high identity and depth is strong evidence. However, a gene found in the DNA is not always expressed, and databases may not contain new determinants. Clinical laboratories therefore confirm genetic findings with phenotypic tests. This tool is intended for education and research, not for clinical reporting.' },
      { label: 'A detected gene always means that the resistance is expressed', why: 'Detection shows that the gene is present in the genome. It does not show that the gene is active or functional.' },
      { label: 'Every VFDB gene should be reported as a cause of the patient\u2019s disease', why: 'Virulence genes indicate potential. Whether disease develops depends on the bacterium, the host and many other factors.' },
      { label: 'The results are unusable because the reads were trimmed', why: 'Quality trimming removes unreliable data and improves the results. The summary shows how many reads were retained.' },
    ],
    correct: 0,
    verify: null,
    note: `This is the competency most worth testing: interpretive humility. Reinforce the three
      caveats — genotype is not phenotype; only acquired genes were screened (point mutations in
      genes such as rpoB or gyrA require different methods); and the tool is for education and
      research, not diagnosis.`,
  },
  {
    id: 'debrief', type: 'final', title: 'Summary',
    html: `
      <p>You screened sequencing reads from a real isolate and found <code>mecA</code>
      (methicillin resistance, reported by both ResFinder and CARD), <code>blaZ</code>
      (penicillin resistance) and the PVL toxin genes <code>lukS-PV</code> and
      <code>lukF-PV</code>. Together, these findings identify a community-associated MRSA.
      The run also carried aminoglycoside (<code>ant(6)-Ia</code>, <code>aph(3')-III</code>)
      and macrolide (<code>msr(A)</code>, <code>mph(C)</code>) resistance genes — USA300's
      typical multiresistance profile. Agreement between two independent databases
      strengthens the key result.</p>
      <p>Some limits apply. You screened for acquired genes; resistance caused by point mutations
      requires different methods. A gene list is not a phenotype, and this analysis is intended
      for education, not clinical reporting.</p>`,
    note: `Closing prompt for discussion: what would you run next? Good answers include phenotypic
      susceptibility testing, a point-mutation screen, molecular typing, and — for outbreak
      questions — comparing this genome against other isolates.`,
  },
];

// ── Rendering ──

function esc(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}

function renderLesson() {
  const area = $('#lesson');
  area.innerHTML = '';
  STEPS.forEach((step, i) => {
    const el = document.createElement('section');
    el.className = 'edu-step';
    el.id = 'step-' + i;
    el.innerHTML = stepHead(i, step) + stepBody(step) + lockNote(step);
    area.appendChild(el);
    wireStep(el, step, i);
  });
  updateAll();
}

function stepHead(i, step) {
  return `
    <div class="edu-step-head">
      <span class="edu-step-num">${i + 1}</span>
      <span class="edu-step-title">${esc(step.title)}</span>
    </div>`;
}

function stepBody(step) {
  let html = `<div class="edu-step-body">${step.html || ''}`;
  if (step.note) html += `<div class="edu-note"><b>Instructor note.</b> ${step.note}</div>`;
  html += '</div>';
  return html;
}

function lockNote(step) {
  const reason = lockReason(step);
  return reason
    ? `<div class="edu-locknote">${esc(reason)}</div>`
    : '';
}

function lockReason(step) {
  if (step.type === 'quiz' && !state.ran) return 'Complete the "Run the analysis" step to unlock this checkpoint.';
  return null;
}

function wireStep(el, step, i) {
  if (step.type === 'load') wireLoadStep(el);
  else if (step.type === 'run') wireRunStep(el);
  else if (step.type === 'quiz') wireQuiz(el, step);
  else if (step.type === 'final') wireFinal(el);
  // read steps need no wiring
}

function updateAll() {
  updateStepStates();
  updateDots();
  renderFinal();
}

function updateStepStates() {
  STEPS.forEach((step, i) => {
    const el = $('#step-' + i);
    if (!el) return;
    el.classList.toggle('locked', !!lockReason(step));
    el.classList.toggle('done', stepDone(step));
  });
}

function stepDone(step) {
  switch (step.type) {
    case 'read': return true;
    case 'load': return state.dataLoaded;
    case 'run': return state.ran;
    case 'quiz': return !!state.answers[step.id]?.done;
    case 'final': return false;
    default: return false;
  }
}

function updateDots() {
  const dots = $('#edu-dots');
  dots.innerHTML = '';
  STEPS.forEach((step, i) => {
    const done = stepDone(step);
    const dot = document.createElement('span');
    dot.className = 'edu-dot' + (done ? ' done' : '');
    dot.textContent = i + 1;
    dot.title = step.title;
    dots.appendChild(dot);
  });
  const firstOpen = STEPS.findIndex(s => !stepDone(s));
  if (firstOpen >= 0) dots.children[firstOpen]?.classList.add('now');
}

// ── Step: load the dataset ──

function wireLoadStep(el) {
  const body = el.querySelector('.edu-step-body');
  const actions = document.createElement('div');
  actions.className = 'edu-actions';
  actions.innerHTML = `
    <button class="btn btn-primary" id="edu-btn-data" type="button">Download reads (51 MB)</button>
    <span class="edu-filechips" id="edu-chips"></span>`;
  body.appendChild(actions);
  const btn = actions.querySelector('#edu-btn-data');
  btn.addEventListener('click', () => loadDataset(btn));
  if (state.dataLoaded) showChips();
}

async function loadDataset(btn) {
  if (state.loadingData) return;
  state.loadingData = true;
  btn.disabled = true;
  const original = btn.textContent;
  try {
    term.push(`Loading example isolate: ${DATASET.name} · ${DATASET.accession}`, 'info');
    for (const spec of DATASET.files) {
      const cacheKey = 'edu:' + spec.name;
      let data = await getCachedDBFile(cacheKey);
      if (data && data.length) {
        term.push(`${spec.name} (cached, ${formatBytes(data.byteLength)})`, 'info');
      } else {
        term.push(`Downloading ${spec.name}…`, 'progress');
        data = await fetchAssetWithProgress(spec.url, spec.bytes, (got, total) => {
          const pct = total ? Math.min(100, Math.round(100 * got / total)) : 0;
          btn.textContent = `${spec.name}  ${formatBytes(got)} / ${formatBytes(total)}  (${pct}%)`;
        });
        try { await cacheDBFile(cacheKey, data); } catch (_) { /* classroom reload will re-download */ }
        term.push(`${spec.name} (${formatBytes(data.byteLength)})`, 'ok');
      }
      state.files[spec.slot] = new File([data], spec.name, { type: 'application/gzip' });
    }
    state.dataLoaded = true;
    btn.textContent = 'Reads loaded';
    showChips();
    refreshRunBtn();
    updateAll();
    term.push('Reads ready. Continue to step 2 and run the analysis.', 'ok');
  } catch (err) {
    term.push('Failed to load data: ' + err.message, 'error');
    btn.disabled = false;
    btn.textContent = original;
  } finally {
    state.loadingData = false;
  }
}

function showChips() {
  const chips = $('#edu-chips');
  if (!chips) return;
  chips.innerHTML = '';
  for (const spec of DATASET.files) {
    const f = state.files[spec.slot];
    if (!f) continue;
    const chip = document.createElement('span');
    chip.className = 'edu-filechip ok';
    chip.textContent = `${f.name} · ${formatBytes(f.size)}`;
    chips.appendChild(chip);
  }
}

// ── Step: run the pipeline ──

function refreshRunBtn() {
  const btn = $('#edu-btn-run');
  if (btn) btn.disabled = !state.dataLoaded || state.running || state.ran;
}

function wireRunStep(el) {
  const body = el.querySelector('.edu-step-body');
  const holder = document.createElement('div');
  holder.innerHTML = `
    <div class="edu-actions">
      <button class="btn btn-primary" id="edu-btn-run" type="button" disabled>Run the analysis</button>
    </div>
    <div class="run-status" id="edu-run-status" hidden>
      <div class="run-progress"><div class="run-progress-bar"></div></div>
      <span class="run-status-text">Starting…</span>
    </div>
    <div id="edu-results"></div>`;
  body.appendChild(holder);

  refreshRunBtn();
  $('#edu-btn-run').addEventListener('click', doRun);
  if (state.ran) renderRunOutputs();
}

async function doRun() {
  if (state.running || !state.dataLoaded || !state.files.r1 || !state.files.r2) return;
  state.running = true;
  const btn = $('#edu-btn-run');
  const status = $('#edu-run-status');
  const statusText = status.querySelector('.run-status-text');
  btn.disabled = true;
  btn.innerHTML = '<span class="spinner"></span>Running…';
  status.hidden = false;
  document.getElementById('terminal-card').classList.remove('collapsed');

  const onStep = (i, n, text) => {
    statusText.textContent = `Step ${i} of ${n}: ${text}`;
    term.push(`[${i}/${n}] ${text}`, 'progress');
  };
  setFastpLog(t => term.push(t, 'info'));

  try {
    const files = [state.files.r1, state.files.r2];
    const { qc, results } = await runComprehensive(files, 'paired', { onStep });
    state.qc = qc;
    state.runFiles = results;
    for (const r of results) {
      state.rows[r.database] = (r.exitCode === 0 && r.resTable) ? parseResFile(r.resTable).rows : [];
    }
    const failed = results.filter(r => r.exitCode !== 0);
    if (failed.length === results.length) throw new Error('all database steps failed — see Logs');
    state.ran = true;
    statusText.textContent = 'Analysis complete.';
    renderRunOutputs();
    term.push('Analysis complete. Continue with the checkpoints below.', 'ok');
  } catch (err) {
    state.runError = err.message;
    statusText.textContent = 'Analysis failed — open the Logs panel for details.';
    term.push('Error: ' + err.message, 'error');
    btn.disabled = false;
    btn.textContent = 'Try again';
  } finally {
    state.running = false;
    refreshRunBtn();
    updateAll();
  }
}

function renderRunOutputs() {
  const area = $('#edu-results');
  if (!area) return;
  const btn = $('#edu-btn-run');
  if (btn) { btn.disabled = true; btn.textContent = 'Analysis complete ✓'; }
  if (state.runError && !state.runFiles?.some(r => r.exitCode === 0)) {
    area.innerHTML = `<div class="edu-fb bad">The analysis failed: ${esc(state.runError)}.
      Check the Logs panel — then reload the page and try again.</div>`;
    return;
  }

  let html = '<div id="edu-out"></div>';
  area.innerHTML = html;
  const out = area.querySelector('#edu-out');

  // QC card (same summary as the main app's simple mode)
  const metrics = summariseQc(state.qc);
  const verdict = qcVerdict(metrics);
  const qcCard = document.createElement('section');
  qcCard.className = 'card';
  qcCard.style.margin = '16px 0';
  qcCard.innerHTML = `
    <div class="card-title comp-title"><span>Quality control</span><span class="comp-src">fastp</span></div>
    <div class="card-body">
      ${metrics ? `
      <div class="qc-grid">
        <div class="qc-cell"><span class="qc-val">${fmtInt(metrics.rawReads)}</span><span class="qc-key">reads in</span></div>
        <div class="qc-cell"><span class="qc-val">${(100 * metrics.retained).toFixed(1)}%</span><span class="qc-key">retained after trimming</span></div>
        <div class="qc-cell"><span class="qc-val">${(100 * (metrics.q30After ?? 0)).toFixed(1)}%</span><span class="qc-key">Q30 (after)</span></div>
        <div class="qc-cell"><span class="qc-val">${(100 * (metrics.gcBefore ?? 0)).toFixed(1)}%</span><span class="qc-key">GC content</span></div>
      </div>
      <p class="qc-verdict qc-${verdict.tone}">${esc(verdict.text)}</p>` : `
      <p class="qc-verdict qc-warn">${esc(verdict.text)}</p>`}
    </div>`;
  out.appendChild(qcCard);

  // One compact table per database
  const dbs = [
    { key: 'resfinder', label: 'Resistance genes (ResFinder)' },
    { key: 'card_homolog', label: 'Resistance homologs (CARD)' },
    { key: 'vfdb_core', label: 'Virulence factors (VFDB)' },
  ];
  for (const db of dbs) {
    const rows = state.rows[db.key] || [];
    const section = document.createElement('div');
    if (!rows.length) {
      section.innerHTML = `
        <div class="edu-res-head"><span>${esc(db.label)}</span><span class="edu-count">no genes detected</span></div>`;
    } else {
      section.innerHTML = `
        <div class="edu-res-head"><span>${esc(db.label)}</span><span class="edu-count">${rows.length} gene${rows.length === 1 ? '' : 's'} detected</span></div>
        <div class="edu-res-wrap"><table class="edu-table">
          <thead><tr><th>Gene</th><th>Identity</th><th>Coverage</th><th>Depth</th><th>p-value</th></tr></thead>
          <tbody></tbody>
        </table></div>`;
      const tbody = section.querySelector('tbody');
      const PREVIEW = 8;
      rows.forEach((row, idx) => {
        const key = KEY_GENES.find(k => k.re.test(row.Template || ''));
        const tr = document.createElement('tr');
        if (key) tr.className = 'edu-hl';
        const [sym, desc] = prettyGene(row.Template, db.key);
        tr.innerHTML = `
          <td><span class="edu-gene">${esc(sym)}</span>${key ? '<span class="edu-hlchip">key gene</span>' : ''}
            ${desc ? `<span class="edu-gene-desc">${esc(desc)}</span>` : ''}</td>
          <td>${esc(row.Template_Identity || '')}%</td>
          <td>${esc(row.Template_Coverage || '')}%</td>
          <td>${esc(row.Depth || '')}</td>
          <td>${esc(row.p_value || '')}</td>`;
        if (idx >= PREVIEW) tr.classList.add('edu-row-limited');
        tbody.appendChild(tr);
      });
      if (rows.length > PREVIEW) {
        const more = document.createElement('button');
        more.className = 'edu-more';
        more.type = 'button';
        more.textContent = `Show all ${rows.length} genes`;
        more.addEventListener('click', () => {
          section.querySelectorAll('.edu-row-limited').forEach(r => r.classList.remove('edu-row-limited'));
          more.remove();
        });
        section.appendChild(more);
      }
    }
    out.appendChild(section);
  }

  // Downloads: one ZIP with the raw result tables + QC + log (for instructors/lab notebooks)
  const dl = document.createElement('div');
  dl.className = 'edu-actions';
  dl.innerHTML = '<button class="btn" id="edu-btn-zip" type="button">Download results (ZIP)</button>';
  out.appendChild(dl);
  dl.querySelector('#edu-btn-zip').addEventListener('click', () => {
    const enc = new TextEncoder();
    const zip = {};
    if (state.qc) zip['qc_report.json'] = enc.encode(JSON.stringify(state.qc, null, 2));
    for (const r of state.runFiles || []) {
      if (r.resTable) zip[r.database + '.res'] = enc.encode(r.resTable);
      for (const [ext, info] of Object.entries(r.files || {})) {
        if (ext === '.res' || ext === '.log') continue;
        zip[r.database + ext] = info.binary
          ? (info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data))
          : enc.encode(info.data);
      }
    }
    zip['lesson.log'] = enc.encode(term.text());
    const a = document.createElement('a');
    a.href = URL.createObjectURL(new Blob([window.fflate.zipSync(zip)], { type: 'application/zip' }));
    a.download = 'openpathogen_lesson1_results.zip';
    a.click();
    URL.revokeObjectURL(a.href);
  });
}

// ResFinder templates look like "mecA_6_BX571856"; VFDB templates are long
// descriptions whose second parenthesised token is the gene symbol.
function prettyGene(template, dbKey) {
  const t = template || '';
  if (dbKey === 'vfdb_core') {
    const parens = t.match(/\(([^)]+)\)/g) || [];
    const sym = parens.length >= 2 ? parens[1].replace(/[()]/g, '') : t.split(/\s+/)[0];
    const rest = t.replace(/^VFG\d+\([^)]*\)\s*/, '');
    return [sym, rest.length > 96 ? rest.slice(0, 96) + '…' : rest];
  }
  if (dbKey === 'card_homolog') {
    const m = t.match(/\(([^)]+)\)/);
    if (m) return [t.split(/\s+/)[0], m[1]];
    return [t.split('_')[0], t];
  }
  return [t.split('_')[0], ''];
}

function fmtInt(n) {
  if (n == null) return '—';
  if (n >= 1e6) return (n / 1e6).toFixed(1) + 'M';
  if (n >= 1e3) return (n / 1e3).toFixed(0) + 'k';
  return String(n);
}

// ── Steps: quizzes ──

function wireQuiz(el, step) {
  const body = el.querySelector('.edu-step-body');
  body.querySelector('.edu-quiz')?.remove(); // instructor-toggle rewiring replaces, not appends
  const quiz = document.createElement('div');
  quiz.className = 'edu-quiz';
  quiz.innerHTML = `
    <p class="edu-q">${step.question}</p>
    <div class="edu-opts"></div>
    <div class="edu-fb-slot"></div>`;
  body.appendChild(quiz);

  const optsEl = quiz.querySelector('.edu-opts');
  const fbSlot = quiz.querySelector('.edu-fb-slot');

  const render = () => {
    const rec = state.answers[step.id] || { wrongPicks: [], done: false };
    optsEl.innerHTML = '';
    step.options.forEach((opt, i) => {
      const btn = document.createElement('button');
      btn.type = 'button';
      btn.className = 'edu-opt';
      if (rec.done && i === step.correct) btn.classList.add('correct');
      if (!rec.done && rec.wrongPicks.includes(i)) btn.classList.add('wrong');
      if (state.instructor && i === step.correct) btn.classList.add('is-answer');
      const key = String.fromCharCode(65 + i);
      btn.innerHTML = `<span class="opt-key">${key}</span><span>${opt.label}</span>`;
      btn.disabled = rec.done || rec.wrongPicks.includes(i);
      btn.addEventListener('click', () => pick(i));
      optsEl.appendChild(btn);
    });
    // Instructor answer-key chip
    optsEl.querySelector('.is-answer')?.insertAdjacentHTML('beforeend',
      '<span class="edu-hlchip" style="margin-left:auto">answer key</span>');
  };

  const pick = (i) => {
    let rec = state.answers[step.id];
    if (!rec) { rec = state.answers[step.id] = { wrongPicks: [], firstTryCorrect: false, done: false }; }
    if (rec.done || rec.wrongPicks.includes(i)) return;
    const correct = i === step.correct;
    if (correct) {
      rec.done = true;
      rec.firstTryCorrect = rec.wrongPicks.length === 0;
    } else {
      rec.wrongPicks.push(i);
    }
    fbSlot.innerHTML = feedbackHtml(step, i, correct);
    render();
    updateAll();
  };

  const feedbackHtml = (step, i, correct) => {
    const chosen = step.options[i];
    let html = `<div class="edu-fb ${correct ? 'ok' : 'bad'}">
      <strong>${correct ? 'Correct.' : 'Not quite. Try again.'}</strong> ${chosen.why}`;
    if (correct && step.verify) {
      const v = step.verify();
      if (v) {
        html += `<div class="edu-verify"><div class="vu-head">In your results</div>${esc(v).replace(/\n/g, '<br>')}</div>`;
      }
    }
    return html + '</div>';
  };

  render();
}

function quizScore() {
  const quizzes = STEPS.filter(s => s.type === 'quiz');
  const answered = quizzes.filter(s => state.answers[s.id]?.done).length;
  const firstTry = quizzes.filter(s => state.answers[s.id]?.firstTryCorrect).length;
  return { answered, firstTry, total: quizzes.length };
}

// ── Step: final card ──

function wireFinal(el) {
  const body = el.querySelector('.edu-step-body');
  const holder = document.createElement('div');
  holder.id = 'edu-final';
  body.appendChild(holder);
  renderFinal();
}

function renderFinal() {
  const holder = $('#edu-final');
  if (!holder) return;
  const { answered, firstTry, total } = quizScore();
  holder.innerHTML = `
    <div class="edu-final-score">${answered < total ? `${answered} of ${total} checkpoints completed` : `${firstTry} of ${total} correct on the first attempt`}</div>
    <p class="edu-final-sub">${answered < total
      ? 'Complete the checkpoints above to finish the lesson.'
      : 'Lesson complete. Each answer was checked against the results of your own analysis.'}</p>
    <div class="edu-next">
      <a class="btn" href="index.html">Analyze your own sample</a>
      <button class="btn btn-ghost" id="edu-btn-restart" type="button">Restart lesson</button>
    </div>
    <p class="edu-cite">Data: S. aureus USA300_TCH1516, ENA <a href="https://www.ebi.ac.uk/ena/browser/view/SRR10341524" target="_blank" rel="noopener">SRR10341524</a> (study PRJNA579343).
    Methods: <a href="https://doi.org/10.1093/bioinformatics/bty560" target="_blank" rel="noopener">fastp</a>,
    <a href="https://doi.org/10.1186/s12859-018-2336-6" target="_blank" rel="noopener">KMA</a>,
    ResFinder, CARD, VFDB. Educational analysis — not a diagnostic result.</p>`;
  holder.querySelector('#edu-btn-restart').addEventListener('click', () => location.reload());
}

// ── Instructor view ──

function setInstructor(on) {
  state.instructor = on;
  localStorage.setItem('op-edu-instructor', on ? '1' : '0');
  document.body.classList.toggle('edu-instructor', on);
  const btn = $('#btn-instructor');
  btn.classList.toggle('on', on);
  btn.textContent = on ? 'Instructor view · ON' : 'Instructor view';
  // Re-render quizzes so answer-key chips appear/disappear.
  STEPS.forEach((step, i) => {
    if (step.type !== 'quiz') return;
    const el = document.getElementById('step-' + i);
    if (el) wireQuiz(el, step);
  });
}

// ── Init ──

document.addEventListener('DOMContentLoaded', async () => {
  // Logs toggle
  const card = document.getElementById('terminal-card');
  document.getElementById('terminal-toggle').addEventListener('click', () => {
    const collapsed = card.classList.toggle('collapsed');
    document.getElementById('terminal-toggle').setAttribute('aria-expanded', String(!collapsed));
  });

  // Instructor toggle (?instructor=1 in the URL turns it on for sharing)
  const wantInstructor = new URLSearchParams(location.search).get('instructor') === '1'
    || localStorage.getItem('op-edu-instructor') === '1';
  $('#btn-instructor').addEventListener('click', () => setInstructor(!state.instructor));

  setLogCallback((msg, level) => term.push(msg, level));
  renderLesson();
  setInstructor(wantInstructor);

  try {
    await Promise.all([initWasm(), loadPhenotypesDB()]);
    if (!self.crossOriginIsolated) {
      term.push('Cross-origin isolation is not active yet — the page reloads once to enable it. If this message persists, the analysis cannot run in this browser.', 'warn');
    }
  } catch (err) {
    term.push('Init failed: ' + err.message, 'error');
  }
});
