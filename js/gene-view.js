/* gene-view.js — per-gene detail viewer built entirely from KMA's existing
 * outputs (no extra WASM): consensus sequence (.fsa), template/consensus
 * alignment (.aln), per-base depth (.mat.gz) and variant calls (.vcf.gz).
 *
 * Everything is parsed lazily and cached on the result object the first time a
 * gene is opened, so large runs don't pay the gunzip/parse cost until needed.
 */

// ── Public entry point ──

// row: a parsed .res row (Template, Template_Identity, …, _phenotype).
// result: the runAnalysis result (carries .files with the KMA outputs).
export function openGeneViewer(row, result) {
  if (!row || !row.Template) return;
  const template = row.Template;
  const viz = getViz(result);

  const consensus = viz.fsa?.get(template) || '';
  const alnBlock = viz.aln?.get(template) || '';
  const matRows = viz.mat?.get(template) || null;
  const variants = viz.vcf?.get(template) || [];

  // Depth in template coordinates: matrix rows are 1:1 with template positions
  // (insertion rows, ref="-", are dropped so POS lines up with VCF POS).
  const depth = matRows ? matRows.filter(r => r.ref !== '-').map(r => r.depth) : [];

  buildModal({ template, row, consensus, alnBlock, depth, variants });
}

// ── Lazy parse + cache ──

function getViz(result) {
  if (result._viz) return result._viz;
  const files = result.files || {};
  const viz = { fsa: null, aln: null, mat: null, vcf: null };
  try { if (files['.fsa']) viz.fsa = indexFasta(asText(files['.fsa'])); } catch (_) {}
  try { if (files['.aln']) viz.aln = indexAln(asText(files['.aln'])); } catch (_) {}
  try { if (files['.mat.gz']) viz.mat = indexMatrix(gunzipText(files['.mat.gz'])); } catch (_) {}
  try { if (files['.vcf.gz']) viz.vcf = indexVcf(gunzipText(files['.vcf.gz'])); } catch (_) {}
  result._viz = viz;
  return viz;
}

function asText(info) {
  if (typeof info.data === 'string') return info.data;
  return window.fflate.strFromU8(info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data));
}

function gunzipText(info) {
  const u8 = info.data instanceof Uint8Array ? info.data : new Uint8Array(info.data);
  return window.fflate.strFromU8(window.fflate.gunzipSync(u8));
}

function indexFasta(text) {
  const map = new Map();
  let name = null, buf = [];
  for (const line of text.split('\n')) {
    if (line[0] === '>') {
      if (name) map.set(name, buf.join(''));
      name = line.slice(1).trim();
      buf = [];
    } else if (name !== null) {
      buf.push(line.trim());
    }
  }
  if (name) map.set(name, buf.join(''));
  return map;
}

function indexAln(text) {
  const map = new Map();
  let name = null, buf = [];
  for (const line of text.split('\n')) {
    if (line[0] === '#') {
      if (name) map.set(name, buf.join('\n').replace(/\n+$/, ''));
      name = line.replace(/^#\s*/, '').trim();
      buf = [];
    } else if (name !== null) {
      buf.push(line);
    }
  }
  if (name) map.set(name, buf.join('\n').replace(/\n+$/, ''));
  return map;
}

// Matrix columns are A C G T N - ; depth = sum of all six.
function indexMatrix(text) {
  const map = new Map();
  let name = null, rows = null;
  for (const line of text.split('\n')) {
    if (line[0] === '#') {
      if (name) map.set(name, rows);
      name = line.replace(/^#\s*/, '').trim();
      rows = [];
    } else if (rows && line.trim()) {
      const p = line.split('\t');
      const c = [+p[1] || 0, +p[2] || 0, +p[3] || 0, +p[4] || 0, +p[5] || 0, +p[6] || 0];
      rows.push({ ref: p[0], c, depth: c[0] + c[1] + c[2] + c[3] + c[4] + c[5] });
    }
  }
  if (name) map.set(name, rows);
  return map;
}

// Keep only real differences (KMA writes a gVCF-style line per position).
function indexVcf(text) {
  const map = new Map();
  for (const line of text.split('\n')) {
    if (!line || line[0] === '#') continue;
    const p = line.split('\t');
    const ref = p[3] || '', alt = p[4] || '';
    if (alt === '.' || alt.toUpperCase() === ref.toUpperCase()) continue;
    const info = p[7] || '';
    const fmt = (p[8] || '').split(':');
    const smp = (p[9] || '').split(':');
    const ftIdx = fmt.indexOf('FT');
    map.has(p[0]) || map.set(p[0], []);
    map.get(p[0]).push({
      pos: +p[1] || 0,
      ref,
      alt: alt.toUpperCase(),
      dp: +((info.match(/DP=(\d+)/) || [])[1]) || 0,
      af: +((info.match(/AF=([\d.]+)/) || [])[1]) || 0,
      ft: ftIdx >= 0 ? (smp[ftIdx] || '') : '',
    });
  }
  return map;
}

function variantType(ref, alt) {
  if (alt.includes('-')) return 'del';
  if (ref.length === 1 && alt.length === 1) return 'SNP';
  return 'indel';
}

// ── Modal UI ──

let activeOverlay = null;
let lastFocus = null; // element to give focus back to when the dialog closes

function buildModal({ template, row, consensus, alnBlock, depth, variants }) {
  closeViewer();

  const haveAny = consensus || alnBlock || depth.length || variants.length;
  const meanDepth = depth.length ? depth.reduce((a, b) => a + b, 0) / depth.length : 0;
  const maxDepth = arrMax(depth);

  const overlay = document.createElement('div');
  overlay.className = 'gv-overlay';
  overlay.innerHTML = `
    <div class="gv-modal" role="dialog" aria-modal="true" aria-label="Gene detail">
      <div class="gv-head">
        <div class="gv-title">${esc(template)}</div>
        <button class="gv-close" type="button" aria-label="Close">×</button>
      </div>
      <div class="gv-stats">${renderStats(row, meanDepth)}</div>
      ${renderAnnotation(row)}
      <div class="gv-tabs">
        <button class="gv-tab active" type="button" data-tab="cov">Coverage</button>
        <button class="gv-tab" type="button" data-tab="seq">Consensus</button>
        <button class="gv-tab" type="button" data-tab="aln">Alignment</button>
        <button class="gv-tab" type="button" data-tab="var">Variants <span class="gv-count">${variants.length}</span></button>
      </div>
      <div class="gv-body">
        ${haveAny ? '' : '<p class="empty">No sequence-level data was produced for this hit.</p>'}
        <div class="gv-panel" data-panel="cov">
          ${depth.length
            ? `<div class="gv-cov-meta">Length ${depth.length} bp · mean depth ${meanDepth.toFixed(1)}× · max ${maxDepth}×</div>
               <canvas class="gv-cov"></canvas>
               <div class="gv-legend"><span class="gv-lg-area"></span> depth <span class="gv-lg-mean"></span> mean <span class="gv-lg-var"></span> variant</div>`
            : '<p class="empty">No depth matrix available.</p>'}
        </div>
        <div class="gv-panel" data-panel="seq" hidden>
          ${consensus
            ? `<div class="gv-seqtools"><span class="gv-hint"><span class="gv-lc">lowercase</span> = low-confidence base</span><button class="btn gv-copy" type="button">Copy FASTA</button></div>
               <div class="gv-seq">${renderConsensus(consensus)}</div>`
            : '<p class="empty">No consensus sequence available.</p>'}
        </div>
        <div class="gv-panel" data-panel="aln" hidden>
          ${alnBlock ? `<pre class="gv-aln">${esc(alnBlock)}</pre>` : '<p class="empty">No alignment available.</p>'}
        </div>
        <div class="gv-panel" data-panel="var" hidden>
          ${renderVariants(variants)}
        </div>
      </div>
    </div>`;

  document.body.appendChild(overlay);
  activeOverlay = overlay;
  lastFocus = document.activeElement;
  overlay.querySelector('.gv-close').focus();

  // Interactions
  overlay.addEventListener('mousedown', (e) => { if (e.target === overlay) closeViewer(); });
  overlay.querySelector('.gv-close').addEventListener('click', closeViewer);
  document.addEventListener('keydown', onKey);

  overlay.querySelectorAll('.gv-tab').forEach(tab => {
    tab.addEventListener('click', () => activateTab(overlay, tab.dataset.tab));
  });

  const copyBtn = overlay.querySelector('.gv-copy');
  if (copyBtn) copyBtn.addEventListener('click', () => {
    const fasta = `>${template}\n${consensus.replace(/(.{60})/g, '$1\n')}\n`;
    navigator.clipboard?.writeText(fasta).then(() => {
      copyBtn.textContent = 'Copied';
      setTimeout(() => { copyBtn.textContent = 'Copy FASTA'; }, 1200);
    });
  });

  // Coverage canvas is the default tab — draw once it has layout width.
  if (depth.length) requestAnimationFrame(() => drawCoverage(overlay, depth, variants, meanDepth));
}

function renderStats(row, meanDepth) {
  const items = [];
  const id = parseFloat(row.Template_Identity);
  const cov = parseFloat(row.Template_Coverage);
  if (!isNaN(id)) items.push(['Identity', id.toFixed(1) + '%']);
  if (!isNaN(cov)) items.push(['Coverage', cov.toFixed(1) + '%']);
  if (row.Template_length) items.push(['Length', row.Template_length.trim() + ' bp']);
  const d = row.Depth ? parseFloat(row.Depth) : meanDepth;
  if (!isNaN(d) && d) items.push(['Depth', d.toFixed(1) + '×']);
  if (row.p_value) items.push(['p-value', row.p_value.trim()]);
  return items.map(([k, v]) => `<span class="gv-stat"><span class="gv-k">${esc(k)}</span><span class="gv-v">${esc(v)}</span></span>`).join('');
}

// Functional annotation from the ResFinder phenotypes table (null for other DBs).
function renderAnnotation(row) {
  const p = row._phenotype;
  if (!p) return '';
  const items = [];
  if (p.class) items.push(['Drug class', esc(p.class)]);
  if (p.phenotype) items.push(['Resistance phenotype', esc(p.phenotype)]);
  if (p.mechanism) items.push(['Mechanism', esc(p.mechanism)]);
  const notes = (p.notes || '').split(';').map(s => s.trim()).filter(Boolean).join('; ');
  if (notes) items.push(['Notes', esc(notes)]);
  if (p.requiredGene) items.push(['Required gene', esc(p.requiredGene)]);
  const ref = renderPmid(p.pmid);
  if (ref) items.push(['Reference', ref]);
  if (!items.length) return '';
  return `<dl class="gv-annot">${items.map(([k, v]) => `<dt>${esc(k)}</dt><dd>${v}</dd>`).join('')}</dl>`;
}

function renderPmid(pmid) {
  const raw = (pmid || '').trim();
  if (!raw) return '';
  return raw.split(',').map(tok => {
    const t = tok.trim();
    return /^\d+$/.test(t)
      ? `<a href="https://pubmed.ncbi.nlm.nih.gov/${t}/" target="_blank" rel="noopener">PMID ${t}</a>`
      : esc(t);
  }).join(', ');
}

function renderConsensus(seq) {
  let html = '';
  for (let i = 0; i < seq.length; i += 60) {
    const chunk = esc(seq.slice(i, i + 60)).replace(/[a-z]+/g, m => `<span class="gv-lc">${m}</span>`);
    html += `<div class="gv-seqline"><span class="gv-pos">${i + 1}</span><span class="gv-bases">${chunk}</span></div>`;
  }
  return html;
}

function renderVariants(variants) {
  if (!variants.length) return '<p class="empty">No variants. The consensus matches the reference allele.</p>';
  let html = '<div class="table-wrap gv-vwrap"><table class="res-table gv-vtable"><thead><tr>'
    + '<th>Position</th><th>Change</th><th>Type</th><th>Depth</th><th>Allele freq</th><th>Filter</th>'
    + '</tr></thead><tbody>';
  for (const v of variants) {
    const ftCls = v.ft && v.ft !== 'PASS' ? ' class="gv-lc"' : '';
    html += `<tr><td>${v.pos}</td><td class="gv-change">${esc(v.ref)} → ${esc(v.alt)}</td>`
      + `<td>${variantType(v.ref, v.alt)}</td><td>${v.dp}×</td><td>${(v.af * 100).toFixed(0)}%</td>`
      + `<td${ftCls}>${esc(v.ft || '-')}</td></tr>`;
  }
  return html + '</tbody></table></div>';
}

function activateTab(overlay, name) {
  overlay.querySelectorAll('.gv-tab').forEach(t => t.classList.toggle('active', t.dataset.tab === name));
  overlay.querySelectorAll('.gv-panel').forEach(p => { p.hidden = (p.dataset.panel !== name); });
  // Canvas needs a (re)draw when it first becomes visible with real width.
  if (name === 'cov') {
    const canvas = overlay.querySelector('.gv-cov');
    if (canvas && canvas._draw) requestAnimationFrame(canvas._draw);
  }
}

// Escape closes; Tab cycles inside the dialog so keyboard focus cannot wander
// into the page behind the overlay.
function onKey(e) {
  if (e.key === 'Escape') { closeViewer(); return; }
  if (e.key !== 'Tab' || !activeOverlay) return;
  const focusable = [...activeOverlay.querySelectorAll('button,[href],input,select,textarea,[tabindex]:not([tabindex="-1"])')]
    .filter(el => !el.disabled && el.offsetParent !== null);
  if (!focusable.length) return;
  const first = focusable[0], last = focusable[focusable.length - 1];
  if (e.shiftKey && document.activeElement === first) { e.preventDefault(); last.focus(); }
  else if (!e.shiftKey && document.activeElement === last) { e.preventDefault(); first.focus(); }
}

function closeViewer() {
  if (!activeOverlay) return;
  document.removeEventListener('keydown', onKey);
  activeOverlay.remove();
  activeOverlay = null;
  if (lastFocus && typeof lastFocus.focus === 'function') lastFocus.focus();
  lastFocus = null;
}

// ── Coverage plot ──

function drawCoverage(overlay, depth, variants, mean) {
  const canvas = overlay.querySelector('.gv-cov');
  if (!canvas) return;
  const draw = () => {
    const cssW = canvas.clientWidth || 600;
    const cssH = 160;
    const dpr = window.devicePixelRatio || 1;
    canvas.width = Math.round(cssW * dpr);
    canvas.height = Math.round(cssH * dpr);
    const ctx = canvas.getContext('2d');
    ctx.setTransform(dpr, 0, 0, dpr, 0, 0);
    ctx.clearRect(0, 0, cssW, cssH);

    const padL = 38, padR = 8, padT = 14, padB = 18;
    const w = cssW - padL - padR;
    const h = cssH - padT - padB;
    const n = depth.length;
    const maxD = Math.max(arrMax(depth), 1);
    const x = i => padL + (n <= 1 ? 0 : (i / (n - 1)) * w);
    const y = d => padT + h - (d / maxD) * h;

    // baseline + frame
    ctx.strokeStyle = '#ccc';
    ctx.lineWidth = 1;
    ctx.beginPath();
    ctx.moveTo(padL, padT + h); ctx.lineTo(padL + w, padT + h);
    ctx.moveTo(padL, padT); ctx.lineTo(padL, padT + h);
    ctx.stroke();

    // depth area
    ctx.beginPath();
    ctx.moveTo(padL, padT + h);
    for (let i = 0; i < n; i++) ctx.lineTo(x(i), y(depth[i]));
    ctx.lineTo(padL + w, padT + h);
    ctx.closePath();
    ctx.fillStyle = 'rgba(0,0,0,0.12)';
    ctx.fill();
    ctx.beginPath();
    for (let i = 0; i < n; i++) (i === 0 ? ctx.moveTo(x(i), y(depth[i])) : ctx.lineTo(x(i), y(depth[i])));
    ctx.strokeStyle = '#666';
    ctx.lineWidth = 1;
    ctx.stroke();

    // mean line
    if (mean > 0) {
      ctx.save();
      ctx.setLineDash([4, 3]);
      ctx.strokeStyle = '#999';
      ctx.beginPath(); ctx.moveTo(padL, y(mean)); ctx.lineTo(padL + w, y(mean)); ctx.stroke();
      ctx.restore();
    }

    // variant ticks
    ctx.fillStyle = '#000';
    ctx.strokeStyle = '#000';
    for (const v of variants) {
      if (v.pos < 1 || v.pos > n) continue;
      const vx = x(v.pos - 1);
      ctx.beginPath(); ctx.moveTo(vx, padT); ctx.lineTo(vx, padT + h); ctx.globalAlpha = 0.35; ctx.stroke(); ctx.globalAlpha = 1;
      ctx.beginPath(); ctx.moveTo(vx - 3, padT - 1); ctx.lineTo(vx + 3, padT - 1); ctx.lineTo(vx, padT + 4); ctx.closePath(); ctx.fill();
    }

    // axis labels
    ctx.fillStyle = '#999';
    const mono = (getComputedStyle(document.documentElement).getPropertyValue('--mono') || '').trim() || 'ui-monospace,Menlo,monospace';
    ctx.font = `11px ${mono}`;
    ctx.textBaseline = 'middle';
    ctx.fillText(String(maxD), 4, padT);
    ctx.fillText('0', 4, padT + h);
    ctx.textBaseline = 'top';
    ctx.fillText('1', padL, padT + h + 4);
    ctx.textAlign = 'right';
    ctx.fillText(String(n), padL + w, padT + h + 4);
    ctx.textAlign = 'left';
  };
  canvas._draw = draw;
  draw();
}

// ── util ──

function arrMax(arr) {
  let m = 0;
  for (let i = 0; i < arr.length; i++) if (arr[i] > m) m = arr[i];
  return m;
}

function esc(s) {
  return String(s).replace(/&/g, '&amp;').replace(/</g, '&lt;').replace(/>/g, '&gt;').replace(/"/g, '&quot;');
}
