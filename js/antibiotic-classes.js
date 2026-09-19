/* antibiotic-classes.js — drug class for every antibiotic name in the CABBAGE
 * snapshot, expressed in ResFinder's class vocabulary (the "Class" column of
 * resfinder_db/phenotypes.txt). The report uses it to split the antibiotics
 * of a gene with several classes (erm(C): macrolides, lincosamides,
 * streptogramins) and to place each CABBAGE rate under the right class, and
 * to tell a portal link that names a drug (mecA → methicillin) from one that
 * names a class (blaZ → "beta-lactam antibiotic").
 *
 * Tokens are lower-case ResFinder class names. Combination products list the
 * class of every active component (an inhibitor combination is still a
 * beta-lactam). Drugs with no acquired-gene class in ResFinder (first-line
 * tuberculosis drugs, antifungals, lipopeptides…) carry their own class name,
 * which never matches a ResFinder gene.
 */

const BL = 'beta-lactam';
const AG = 'aminoglycoside';
const MAC = 'macrolide';
const LIN = 'lincosamide';
const SGA = 'streptogramin a';
const SGB = 'streptogramin b';
const TET = 'tetracycline';
const QUI = 'quinolone';
const AMP = 'amphenicol';
const FOL = 'folate pathway antagonist';
const GLY = 'glycopeptide';
const OXA = 'oxazolidinone';
const POL = 'polymyxin';
const FOS = 'fosfomycin';
const RIF = 'rifamycin';
const FUS = 'steroid antibacterial';
const MUP = 'pseudomonic acid';
const NIM = 'nitroimidazole';
const PLE = 'pleuromutilin';

const TABLE = {
  // beta-lactams: penicillins, cephalosporins, carbapenems, monobactams, inhibitor combinations
  'amoxicillin': [BL], 'amoxicillin-clavulanic acid': [BL], 'ampicillin': [BL], 'ampicillin-sulbactam': [BL],
  'aztreonam': [BL], 'carbenicillin': [BL], 'cefaclor': [BL], 'cefamandole': [BL], 'cefatrizine': [BL],
  'cefazolin': [BL], 'cefepime': [BL], 'cefepime-taniborbactam': [BL], 'cefiderocol': [BL], 'cefixime': [BL],
  'cefmetazole': [BL], 'cefoperazone': [BL], 'cefotaxime': [BL], 'cefotaxime-clavulanic acid': [BL],
  'cefotetan': [BL], 'cefoxitin': [BL], 'cefpimizole': [BL], 'cefpirome': [BL], 'cefpodoxime': [BL],
  'cefpodoxime-clavulanic acid': [BL], 'ceftaroline': [BL], 'ceftazidime': [BL], 'ceftazidime-avibactam': [BL],
  'ceftazidime-clavulanic acid': [BL], 'ceftibuten': [BL], 'ceftiofur': [BL], 'ceftizoxime': [BL],
  'ceftolozane-tazobactam': [BL], 'ceftriaxone': [BL], 'ceftriaxone-cefpodoxime': [BL], 'cefuroxime': [BL],
  'cephalexin': [BL], 'cephalothin': [BL], 'dicloxacillin': [BL], 'doripenem': [BL], 'ertapenem': [BL],
  'imipenem': [BL], 'imipenem-relebactam': [BL], 'mecillinam': [BL], 'meropenem': [BL],
  'meropenem-vaborbactam': [BL], 'methicillin': [BL], 'mezlocillin': [BL], 'oxacillin': [BL], 'penicillin': [BL],
  'penicillin-clavulanic acid': [BL], 'piperacillin': [BL], 'piperacillin-tazobactam': [BL], 'sulbactam': [BL],
  'temocillin': [BL], 'ticarcillin': [BL], 'ticarcillin-clavulanic acid': [BL],
  // aminoglycosides (spectinomycin: aadA-type genes are classed as aminoglycoside in ResFinder)
  'amikacin': [AG], 'apramycin': [AG], 'gentamicin': [AG], 'kanamycin': [AG], 'neomycin': [AG],
  'netilmicin': [AG], 'plazomicin': [AG], 'spectinomycin': [AG], 'streptomycin': [AG], 'tobramycin': [AG],
  // macrolides, lincosamides, streptogramins
  'azithromycin': [MAC], 'clarithromycin': [MAC], 'erythromycin': [MAC], 'telithromycin': [MAC],
  'tilmicosin': [MAC], 'tylosin': [MAC], 'clindamycin': [LIN], 'quinupristin-dalfopristin': [SGA, SGB],
  // tetracyclines and glycylcyclines
  'tetracycline': [TET], 'doxycycline': [TET], 'minocycline': [TET], 'chlortetracycline': [TET],
  'oxytetracycline': [TET], 'tigecycline': [TET], 'eravacycline': [TET], 'omadacycline': [TET],
  // quinolones
  'ciprofloxacin': [QUI], 'levofloxacin': [QUI], 'moxifloxacin': [QUI], 'norfloxacin': [QUI], 'ofloxacin': [QUI],
  'nalidixic acid': [QUI], 'enrofloxacin': [QUI], 'gatifloxacin': [QUI], 'pefloxacin': [QUI],
  'sparfloxacin': [QUI], 'trovafloxacin': [QUI], 'delafloxacin': [QUI],
  // amphenicols
  'chloramphenicol': [AMP], 'florfenicol': [AMP],
  // folate pathway antagonists
  'sulfamethoxazole': [FOL], 'sulfisoxazole': [FOL], 'sulfafurazole': [FOL], 'sulfadimethoxine': [FOL],
  'sulphathiazole': [FOL], 'trimethoprim': [FOL], 'trimethoprim-sulfamethoxazole': [FOL],
  'trimethoprim-sulfobactam': [FOL],
  // glycopeptides, oxazolidinones, polymyxins, others with ResFinder classes
  'vancomycin': [GLY], 'teicoplanin': [GLY], 'dalbavancin': [GLY], 'telavancin': [GLY],
  'linezolid': [OXA], 'tedizolid': [OXA],
  'colistin': [POL], 'polymyxin b': [POL],
  'fosfomycin': [FOS], 'rifampin': [RIF], 'rifabutin': [RIF], 'fusidic acid': [FUS], 'mupirocin': [MUP],
  'metronidazole': [NIM], 'tiamulin': [PLE],
  // no acquired-gene class in ResFinder
  'daptomycin': ['lipopeptide'], 'bacitracin zinc': ['bacitracin'], 'nitrofurantoin': ['nitrofuran'],
  'furazolidone': ['nitrofuran'], 'fidaxomicin': ['macrocyclic'], 'avilamycin': ['orthosomycin'],
  'zoliflodacin': ['spiropyrimidinetrione'], 'azidothymidine': ['nucleoside analogue'],
  'isoniazid': ['antituberculosis'], 'ethionamide': ['antituberculosis'], 'prothionamide': ['antituberculosis'],
  'ethambutol': ['antituberculosis'], 'pyrazinamide': ['antituberculosis'], 'nicotinamide': ['antituberculosis'],
  'pyrazinamide-nicotinamide': ['antituberculosis'], 'bedaquiline': ['antituberculosis'],
  'delamanid': ['antituberculosis'], 'clofazimine': ['antituberculosis'], 'capreomycin': ['antituberculosis'],
  'cycloserine': ['antituberculosis'], 'para-aminosalicylic acid': ['antituberculosis'],
  'pentizidone': ['antituberculosis'],
  'fluconazole': ['azole antifungal'], 'voriconazole': ['azole antifungal'], 'clotrimazole': ['azole antifungal'],
};

// Set of class tokens for a CABBAGE antibiotic name, or null when unknown.
export function antibioticClasses(name) {
  const entry = TABLE[String(name || '').trim().toLowerCase()];
  return entry ? new Set(entry) : null;
}

// Set of class tokens from a ResFinder "Class" string ("Macrolide, Lincosamide,
// Streptogramin B" → three tokens). Empty when no annotation is known.
export function geneClasses(classStr) {
  const out = new Set();
  for (let tok of String(classStr || '').split(',')) {
    tok = tok.trim().toLowerCase();
    if (!tok) continue;
    if (tok === 'tetracyclines') tok = TET;
    if (tok === 'phenicol') tok = AMP;
    out.add(tok);
  }
  return out;
}

// True when a gene's classes and an antibiotic's classes overlap.
export function sameClass(geneCls, abCls) {
  if (!geneCls || !abCls) return false;
  for (const c of geneCls) if (abCls.has(c)) return true;
  return false;
}

// ── Portal gene→antibiotic links ──
//
// CABBAGE attaches to every AMRFinderPlus call the antibiotic names the
// gene confers resistance to (genotype-antibiotic_name). Most are specific
// drugs; some are class names. This expands one link name to the phenotype
// antibiotic names it covers: a drug covers itself and combination products
// that contain it; a class name covers the drugs of that class.
const LINK_CLASS = {
  // mecA/mecC carry the AMRFinderPlus subclass METHICILLIN. Oxacillin and
  // cefoxitin are the phenotypic tests for mecA-mediated resistance (CLSI,
  // EUCAST), so the methicillin link covers them and dicloxacillin.
  'methicillin': (n) => /^(methicillin|oxacillin|dicloxacillin|cefoxitin)$/.test(n),
  'beta-lactam antibiotic': (n, c) => c.has(BL),
  'cephalosporin': (n, c) => c.has(BL) && /^(cef|ceph)/.test(n),
  'carbapenem': (n, c) => c.has(BL) && /penem/.test(n),
  'fluoroquinolone antibiotic': (n, c) => c.has(QUI),
  'sulfonamide antibiotic': (n, c) => c.has(FOL) && /sul[fp]/.test(n),
  'aminoglycoside antibiotic': (n, c) => c.has(AG),
  'streptogramin b antibiotic': (n, c) => c.has(SGB),
  'streptogramin antibiotic': (n, c) => c.has(SGA) || c.has(SGB),
  'rifamycin antibiotic': (n, c) => c.has(RIF),
  'lincosamide antibiotic': (n, c) => c.has(LIN),
  'oxazolidinone antibiotic': (n, c) => c.has(OXA),
  'tetracycline': (n) => /^(tetracycline|doxycycline|minocycline|chlortetracycline|oxytetracycline)$/.test(n),
  'pleuromutilin': (n, c) => c.has(PLE),
  'kanamycin a': (n) => n === 'kanamycin',
};
const linkCache = new Map();

export function linkedAntibiotics(linkName) {
  const key = String(linkName || '').trim().toLowerCase();
  let set = linkCache.get(key);
  if (set) return set;
  set = new Set();
  const rule = LINK_CLASS[key];
  if (rule) {
    for (const [name, classes] of Object.entries(TABLE)) if (rule(name, new Set(classes))) set.add(name);
  } else if (key) {
    for (const name of Object.keys(TABLE)) {
      if (name === key || name.split('-').includes(key)) set.add(name);
    }
  }
  linkCache.set(key, set);
  return set;
}
