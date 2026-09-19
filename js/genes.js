/* genes.js — what a KMA hit is, in words people can read.
 *
 * KMA reports each hit by its template name, the database's FASTA header:
 *   ResFinder  mecA_1_NC_002745
 *   CARD       gb|KC243783.1|+|0-2007|ARO:3000617|mecA [Staphylococcus aureus]
 *   VFDB       VFG001276(gb|WP_024937002) (lukF-PV) Panton-Valentine leukocidin
 *              chain F precursor [PVL (VF0018) - Exotoxin (VFC0235)] [S. aureus MW2]
 *
 * This module turns those into a gene name plus what the gene does. For
 * resistance genes that is the drug classes and antibiotics from ResFinder's
 * curated phenotype table (CARD models are matched to ResFinder genes by
 * name); for virulence genes, the factor and category VFDB assigns.
 */

import { lookupGenePhenotype, resfinderNameIndex } from './resfinder-db.js';
import { antibioticClasses, geneClasses } from './antibiotic-classes.js';

// ── Template names ──

export function parseCardTemplate(template) {
  const bars = String(template || '').split('|');
  if (bars.length < 6) return { name: String(template || ''), aro: null };
  const name = bars.slice(5).join('|').split(' [')[0].trim();
  const aro = (bars[4].match(/ARO:(\d+)/) || [])[1] || null;
  return { name: name || String(template), aro };
}

// Some VFDB headers carry no "(gb|…)" accession, a few gene symbols nest
// parentheses ("msrA/B(pilB)"), and symbols ending in "*" are locus tags.
const VFDB_RE = /^(VFG\d+)(?:\([^)]*\))?\s+\(((?:[^()]|\([^()]*\))*)\)\s+(.*?)\s+\[(.+?) \((VF\d+)\) - (.+?) \((VFC\d+)\)\]\s+\[(.+)\]\s*$/;

export function parseVfdbTemplate(template) {
  const m = String(template || '').match(VFDB_RE);
  if (!m) return null;
  const gene = m[2].replace(/\*$/, '').trim();
  return {
    id: m[1],
    gene: gene && gene !== '-' ? gene : m[1],
    product: m[3],
    factor: m[4],
    category: m[6],
    organism: m[8],
  };
}

// The gene name a database uses for a template.
export function geneName(template, db) {
  if (db === 'card_homolog') return parseCardTemplate(template).name;
  if (db === 'vfdb_core') return parseVfdbTemplate(template)?.gene || String(template || '');
  return String(template || '').split('_')[0];
}

// ── Drug classes ──

// Report groups, in report order. tokens are ResFinder class names as
// geneClasses() returns them; hint names a familiar drug of the class.
export const CLASS_GROUPS = [
  { key: 'beta-lactam', label: 'Beta-lactams', tokens: ['beta-lactam'] },
  { key: 'aminoglycoside', label: 'Aminoglycosides', tokens: ['aminoglycoside'] },
  { key: 'quinolone', label: 'Quinolones', hint: 'ciprofloxacin', tokens: ['quinolone'] },
  { key: 'folate', label: 'Sulfonamides and trimethoprim', tokens: ['folate pathway antagonist'] },
  { key: 'tetracycline', label: 'Tetracyclines', tokens: ['tetracycline'] },
  { key: 'macrolide', label: 'Macrolides', hint: 'erythromycin', tokens: ['macrolide'] },
  { key: 'lincosamide', label: 'Lincosamides', hint: 'clindamycin', tokens: ['lincosamide'] },
  { key: 'streptogramin', label: 'Streptogramins', tokens: ['streptogramin a', 'streptogramin b'] },
  { key: 'amphenicol', label: 'Amphenicols', hint: 'chloramphenicol', tokens: ['amphenicol'] },
  { key: 'glycopeptide', label: 'Glycopeptides', hint: 'vancomycin', tokens: ['glycopeptide'] },
  { key: 'oxazolidinone', label: 'Oxazolidinones', hint: 'linezolid', tokens: ['oxazolidinone'] },
  { key: 'polymyxin', label: 'Polymyxins', hint: 'colistin', tokens: ['polymyxin'] },
  { key: 'fosfomycin', label: 'Fosfomycin', tokens: ['fosfomycin'] },
  { key: 'rifamycin', label: 'Rifamycins', hint: 'rifampicin', tokens: ['rifamycin'] },
  { key: 'fusidic', label: 'Fusidic acid', tokens: ['steroid antibacterial'] },
  { key: 'mupirocin', label: 'Mupirocin', tokens: ['pseudomonic acid'] },
  { key: 'nitroimidazole', label: 'Nitroimidazoles', hint: 'metronidazole', tokens: ['nitroimidazole'] },
  { key: 'pleuromutilin', label: 'Pleuromutilins', tokens: ['pleuromutilin'] },
];

const GROUP_BY_TOKEN = new Map();
for (const g of CLASS_GROUPS) for (const t of g.tokens) GROUP_BY_TOKEN.set(t, g);

// The report group of a class token. Classes outside the list (ionophores…)
// get a group of their own, which is never listed as "not found".
export function classGroup(token) {
  return GROUP_BY_TOKEN.get(token)
    || { key: token, label: token.charAt(0).toUpperCase() + token.slice(1), tokens: [token], extra: true };
}

// Antibiotics in ResFinder's phenotype table that the CABBAGE class table
// (antibiotic-classes.js) does not name. Only needed to split the antibiotics
// of a gene with several classes (erm(C): macrolide, lincosamide, streptogramin B).
const AG = 'aminoglycoside';
const EXTRA_CLASSES = {
  quinupristin: ['streptogramin b'], 'pristinamycin ia': ['streptogramin b'], 'virginiamycin s': ['streptogramin b'],
  dalfopristin: ['streptogramin a'], 'pristinamycin iia': ['streptogramin a'], 'virginiamycin m': ['streptogramin a'],
  lincomycin: ['lincosamide'], valnemulin: ['pleuromutilin'], fluoroquinolone: ['quinolone'], rifampicin: ['rifamycin'],
  spiramycin: ['macrolide'], oleandomycin: ['macrolide'], carbomycin: ['macrolide'],
  dibekacin: [AG], sisomicin: [AG], butirosin: [AG], isepamicin: [AG], lividomycin: [AG], paromomycin: [AG],
  ribostamycin: [AG], fortimicin: [AG], arbekacin: [AG], astromicin: [AG], hygromycin: [AG], kasugamycin: [AG],
  butiromycin: [AG],
  narasin: ['ionophores'], salinomycin: ['ionophores'], maduramicin: ['ionophores'],
};

const abKey = (ab) => String(ab).trim().toLowerCase().replace(/\s*\+\s*/g, '-');

function classesOfAntibiotic(ab) {
  const k = abKey(ab);
  return antibioticClasses(k) || (EXTRA_CLASSES[k] ? new Set(EXTRA_CLASSES[k]) : null);
}

// Split a gene's antibiotics over its classes. Antibiotics of a single-class
// gene, and any that cannot be placed, go to the first (main) class.
function antibioticsByClass(classes, antibiotics) {
  const out = new Map(classes.map(c => [c, []]));
  for (const ab of antibiotics) {
    let placed = false;
    if (classes.length > 1) {
      const abc = classesOfAntibiotic(ab);
      if (abc) for (const c of classes) if (abc.has(c)) { out.get(c).push(ab); placed = true; }
    }
    if (!placed && classes.length) out.get(classes[0]).push(ab);
  }
  return out;
}

// ── Matching CARD models to ResFinder genes ──
//
// "PC1_blaZ" → blaZ, "SHV-52" → blaSHV-52, "APH(3')-IIIa" → aph(3')-III,
// "ErmC" → erm(C). An allele of a beta-lactamase family ResFinder knows
// ("LEN-17") is recognised as a beta-lactam gene of unknown spectrum; vanA,
// vanB… are matched to ResFinder's van operons; other variants to their
// family's class. Apostrophes are kept, so aph(3') and aph(3'') stay
// different genes.

// A name without its allele letter: after a digit ("blaTEM-1B" →
// "blaTEM-1") or after an upper-case roman numeral ("aph(3')-IIIa" →
// "aph(3')-III"). Case matters: in "aph(3')-III" the last I is a numeral.
const dropAllele = (name) => name.replace(/(\d)[A-Za-z]$|([IVX])[a-z]$/, '$1$2');
// A name without its variant number ("fosA6" → "fosA", "AAC(6')-Ib9" → "AAC(6')-Ib").
const dropNumber = (name) => name.replace(/[-_.]?\d+$/, '');

// True when one name is the other plus a variant suffix: aph(3')-III and
// aph(3')-IIIa, aac(6')-Ib and AAC(6')-Ib9, fosA and fosA6. Different genes
// of a family (sul1 and sul2, aac(6')-Ia and aac(6')-Ib) stay apart.
export function sameVariantFamily(a, b) {
  const A = a.toLowerCase(), B = b.toLowerCase();
  const parents = (n) => [dropNumber(n), dropAllele(n)].map(x => x.toLowerCase());
  return A === B || parents(a).includes(B) || parents(b).includes(A);
}

// The family a gene's class can be read from: "aadA" and aadA1, sul4 and sul1.
const familyBase = (name) => dropAllele(dropNumber(name)).toLowerCase();

export function matchCardToResfinder(cardName) {
  const idx = resfinderNameIndex();
  if (!idx) return null;
  const n = String(cardName || '').replace(/’/g, "'").trim();
  const cands = [n, ...n.split('_').filter(t => /^bla/i.test(t))];
  const paren = n.match(/^([a-z]{3})([a-z][a-z0-9]*)$/i);
  if (paren) cands.push(`${paren[1]}(${paren[2]})`);
  for (const c of cands) {
    for (const v of [c, 'bla' + c]) {
      const hit = idx.genes.get(v.toLowerCase()) || idx.genes.get(dropAllele(v).toLowerCase());
      if (hit) return { name: hit.name, entry: hit.entry, how: 'name' };
    }
  }
  const fam = n.replace(/^bla/i, '').match(/^([a-z]+(?:-[a-z]+)?)-?\d/i);
  if (fam && idx.blaFamilies.has(fam[1].toLowerCase())) {
    return { name: 'bla' + n.replace(/^bla/i, ''), entry: null, how: 'family', classes: ['beta-lactam'] };
  }
  const van = n.match(/^van([a-z])$/i);
  if (van) {
    const L = van[1].toLowerCase();
    for (const [k, v] of idx.genes) {
      if (k === `vanh${L}x` || k === `van${L}xy` || k === `van${L}`) return { name: v.name, entry: v.entry, how: 'name' };
    }
    return { name: n, entry: null, how: 'family', classes: ['glycopeptide'] };
  }
  // Families ResFinder lists only as numbered variants ("aadA" → aadA1…,
  // "AAC(6')-Ib9" → aac(6')-Ib): the class is known, this variant's
  // antibiotics are not.
  const base = familyBase(n);
  if (base.length >= 3) {
    for (const v of idx.genes.values()) {
      if (familyBase(v.name) === base) {
        return { name: n, entry: null, how: 'family', classes: [...geneClasses(v.entry.class)] };
      }
    }
  }
  return null;
}

// ── Hit description ──

// Caveats ResFinder records in its Notes column that change what a hit means.
function caveats(notes) {
  const out = [];
  for (const m of String(notes || '').matchAll(/(Natural in [^,;]+|Not functional in [^,;]+|Chromosomal[^,;]*)/gi)) {
    const t = m[1].trim().replace(/\.$/, '');
    if (!out.includes(t)) out.push(t);
  }
  return out;
}

const MEC = /^mec[abc]$/i;

// A short mechanism label for well-known gene types, from the ResFinder
// spectrum and Ambler class. Neutral wording: it names the enzyme type.
function mechanismTag(gene, abs, notes) {
  if (MEC.test(gene)) return 'methicillin resistance';
  if (!/^bla/i.test(gene) || !abs.length) return null;
  const keys = abs.map(abKey);
  if (keys.some(k => CARBAPENEM.test(k))) return 'carbapenemase';
  if (/\bClass C\b/.test(notes || '')) return 'AmpC';
  if (keys.some(k => CEPH_GEN[k] === 3)) return 'ESBL';
  return null;
}

// Everything the report and tables show about one hit:
//   name     the gene as this database calls it
//   gene     the name the report uses (a CARD hit takes its ResFinder name)
//   classes  ResFinder class tokens; empty for virulence genes and for CARD
//            models without a ResFinder match
//   byClass  class token → antibiotics
//   spectrumKnown  false when ResFinder does not list the antibiotics
export function describeHit(template, db, row) {
  const hit = {
    db, template, name: geneName(template, db), gene: null,
    classes: [], byClass: new Map(), antibiotics: [], spectrumKnown: true,
    notes: [], tag: null, match: null, vf: null, aro: null, entry: null,
  };
  if (db === 'vfdb_core') {
    hit.vf = parseVfdbTemplate(template);
    hit.gene = hit.name;
    return hit;
  }
  let entry = null;
  if (db === 'card_homolog') {
    const card = parseCardTemplate(template);
    hit.aro = card.aro;
    const m = matchCardToResfinder(card.name);
    hit.gene = m ? m.name : card.name;
    hit.match = m ? m.how : null;
    entry = m?.entry || null;
    if (m?.classes) hit.classes = m.classes.slice();
  } else {
    entry = row?._phenotype || lookupGenePhenotype(template);
    hit.gene = hit.name;
    hit.match = entry ? 'name' : null;
  }
  hit.entry = entry;
  if (entry) {
    hit.classes = [...geneClasses(entry.class)];
    const abs = String(entry.phenotype || '').split(',').map(s => s.trim()).filter(Boolean);
    if (abs.some(a => /^unknown\b/i.test(a))) hit.spectrumKnown = false;
    hit.antibiotics = abs.filter(a => !/^unknown\b/i.test(a));
    hit.notes = caveats(entry.notes);
    hit.tag = mechanismTag(hit.gene, hit.antibiotics, entry.notes);
  } else if (MEC.test(hit.gene)) {
    hit.tag = 'methicillin resistance';
  }
  if (hit.classes.length && !hit.antibiotics.length) hit.spectrumKnown = false;
  hit.byClass = antibioticsByClass(hit.classes, hit.antibiotics);
  return hit;
}

// ── Antibiotics in words ──

const PEN = /^(amoxicillin|ampicillin|penicillin|piperacillin|ticarcillin|temocillin|mecillinam|oxacillin|cloxacillin|dicloxacillin|flucloxacillin|methicillin|carbenicillin|mezlocillin|azlocillin|nafcillin)$/;
const PEN_INH = /^(amoxicillin|ampicillin|piperacillin|ticarcillin|penicillin)-(clavulanic acid|sulbactam|tazobactam)$/;
const CEPH_GEN = {
  cephalothin: 1, cefalotin: 1, cefazolin: 1, cephalexin: 1, cefalexin: 1, cefadroxil: 1, cefatrizine: 1,
  cephradine: 1, cefapirin: 1,
  cefuroxime: 2, cefaclor: 2, cefamandole: 2, cefprozil: 2, cefonicid: 2,
  cefotaxime: 3, ceftriaxone: 3, ceftazidime: 3, cefixime: 3, cefpodoxime: 3, ceftiofur: 3, ceftizoxime: 3,
  cefoperazone: 3, ceftibuten: 3, cefdinir: 3, cefditoren: 3, cefpimizole: 3,
  cefepime: 4, cefpirome: 4,
  ceftaroline: 5, ceftobiprole: 5,
};
const CEPHAMYCIN = /^(cefoxitin|cefotetan|cefmetazole)$/;
const CEPH_INH = /^cef[a-z]+-(clavulanic acid|avibactam|tazobactam|taniborbactam|enmetazobactam)$/;
const CARBAPENEM = /^(imipenem|meropenem|ertapenem|doripenem|biapenem|faropenem)$/;
const CARB_INH = /^(imipenem|meropenem)-(relebactam|vaborbactam)$/;
const ORD = ['', '1st', '2nd', '3rd', '4th', '5th'];

function andList(items) {
  if (items.length < 2) return items.join('');
  return items.slice(0, -1).join(', ') + ' and ' + items[items.length - 1];
}

// "3rd- and 4th-generation cephalosporins"
function generationLabel(gens) {
  const parts = gens.map(g => ORD[g]);
  if (parts.length === 1) return `${parts[0]}-generation cephalosporins`;
  return andList(parts.map((p, i) => (i < parts.length - 1 ? p + '-' : p))) + '-generation cephalosporins';
}

// Beta-lactam antibiotics grouped the way they are prescribed: penicillinases,
// ESBLs, AmpC and carbapenemases differ exactly in these groups.
function betaLactamPhrase(antibiotics) {
  const keys = antibiotics.map(abKey);
  const has = (re) => keys.some(k => re.test(k));
  const parts = [];
  if (has(PEN)) parts.push('penicillins');
  if (has(PEN_INH)) parts.push('penicillin–inhibitor combinations');
  const gens = [...new Set(keys.map(k => CEPH_GEN[k]).filter(Boolean))].sort();
  if (gens.length) parts.push(generationLabel(gens));
  if (has(CEPHAMYCIN)) parts.push('cephamycins (cefoxitin)');
  if (has(CEPH_INH)) parts.push('cephalosporin–inhibitor combinations');
  if (keys.includes('cefiderocol')) parts.push('cefiderocol');
  if (keys.some(k => /^aztreonam/.test(k))) parts.push('aztreonam');
  if (has(CARBAPENEM)) parts.push('carbapenems');
  if (has(CARB_INH)) parts.push('carbapenem–inhibitor combinations');
  const known = (k) => PEN.test(k) || PEN_INH.test(k) || CEPH_GEN[k] || CEPHAMYCIN.test(k) || CEPH_INH.test(k)
    || k === 'cefiderocol' || /^aztreonam/.test(k) || CARBAPENEM.test(k) || CARB_INH.test(k);
  keys.forEach((k, i) => { if (!known(k)) parts.push(lowerFirst(antibiotics[i])); });
  return parts.join(', ');
}

// Familiar drugs first; the rest keep ResFinder's order.
const COMMON = ['gentamicin', 'tobramycin', 'amikacin', 'streptomycin', 'kanamycin', 'neomycin', 'spectinomycin',
  'netilmicin', 'ciprofloxacin', 'levofloxacin', 'nalidixic acid', 'erythromycin', 'azithromycin', 'clarithromycin',
  'clindamycin', 'lincomycin', 'tetracycline', 'doxycycline', 'minocycline', 'tigecycline', 'chloramphenicol',
  'florfenicol', 'sulfamethoxazole', 'trimethoprim', 'vancomycin', 'teicoplanin', 'linezolid', 'tedizolid',
  'colistin', 'fosfomycin', 'rifampicin', 'fusidic acid', 'mupirocin', 'metronidazole', 'tiamulin',
  'quinupristin', 'dalfopristin'];

const lowerFirst = (s) => String(s).charAt(0).toLowerCase() + String(s).slice(1);

function listPhrase(antibiotics, cap) {
  let abs = antibiotics.slice();
  const keys = abs.map(abKey);
  // "Fluoroquinolone" is a class name in the list; drop it next to real drugs.
  if (keys.includes('fluoroquinolone')) {
    abs = keys.some(k => k !== 'fluoroquinolone' && (antibioticClasses(k)?.has('quinolone')))
      ? abs.filter(a => abKey(a) !== 'fluoroquinolone')
      : abs.map(a => (abKey(a) === 'fluoroquinolone' ? 'fluoroquinolones' : a));
  }
  const rank = (a) => { const i = COMMON.indexOf(abKey(a)); return i < 0 ? COMMON.length : i; };
  abs = abs.map((a, i) => [a, i]).sort((x, y) => rank(x[0]) - rank(y[0]) || x[1] - y[1]).map(x => lowerFirst(x[0]));
  if (cap && abs.length > cap) return `${abs.slice(0, cap).join(', ')} and ${abs.length - cap} more`;
  return abs.join(', ');
}

// Antibiotics of one report group (group key from CLASS_GROUPS) in words.
export function antibioticPhrase(groupKey, antibiotics, { cap = 6 } = {}) {
  const uniq = [...new Set(antibiotics)];
  if (!uniq.length) return '';
  return groupKey === 'beta-lactam' ? betaLactamPhrase(uniq) : listPhrase(uniq, cap);
}

// One-line "what it does" for the gene tables.
export function describeFunction(hit) {
  if (hit.vf) return [hit.vf.product, `${hit.vf.factor} (${hit.vf.category.toLowerCase()})`].filter(Boolean).join(' · ');
  if (!hit.classes.length) return '';
  const groups = new Map();
  for (const [token, abs] of hit.byClass) {
    const g = classGroup(token);
    if (!groups.has(g.key)) groups.set(g.key, { g, abs: [] });
    groups.get(g.key).abs.push(...abs);
  }
  const parts = [...groups.values()].map(({ g, abs }) => {
    const words = antibioticPhrase(g.key, abs, { cap: 4 });
    return words ? `${g.label}: ${words}` : g.label;
  });
  let text = parts.join('; ');
  if (!hit.spectrumKnown) text += ' (antibiotics not listed for this allele)';
  return text;
}

export function cardAroUrl(aro) {
  return aro ? `https://card.mcmaster.ca/aro/${aro}` : null;
}
