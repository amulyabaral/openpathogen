let phenotypesMap = null;
let phenotypesPromise = null;

export async function loadPhenotypesDB() {
  if (phenotypesMap) return phenotypesMap;
  if (phenotypesPromise) return phenotypesPromise;

  phenotypesPromise = fetch('resfinder_db/phenotypes.txt')
    .then(r => r.text())
    .then(text => {
      phenotypesMap = parsePhenotypes(text);
      console.log(`[ResFinder DB] Loaded ${Object.keys(phenotypesMap).length} gene entries`);
      return phenotypesMap;
    })
    .catch(err => {
      console.error('[ResFinder DB] Failed to load phenotypes:', err);
      phenotypesMap = {};
      return phenotypesMap;
    });

  return phenotypesPromise;
}

function parsePhenotypes(text) {
  const map = {};
  const lines = text.trim().split('\n');
  if (lines.length < 2) return map;

  const headers = lines[0].split('\t');

  for (let i = 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line) continue;

    const cols = line.split('\t');
    const geneAccession = cols[0];
    if (!geneAccession) continue;

    // The gene name is the first part of the accession (e.g., "aac(2')-Ic_1_U72714" -> "aac(2')-Ic")
    const geneName = geneAccession.split('_')[0];

    const entry = {
      geneAccession,
      class: cols[1] || '',
      phenotype: cols[2] || '',
      pmid: cols[3] || '',
      mechanism: cols[4] || '',
      notes: cols[5] || '',
      requiredGene: cols[6] || '',
    };

    // Store by full accession and by gene name
    map[geneAccession] = entry;
    if (geneName && geneName !== geneAccession) {
      if (!map[geneName]) map[geneName] = [];
      if (Array.isArray(map[geneName])) {
        map[geneName].push(entry);
      } else {
        map[geneName] = [map[geneName], entry];
      }
    }
  }

  return map;
}

// Case-insensitive gene-name index, built on first use: "blashv-52" →
// { name: 'blaSHV-52', entry }. blaFamilies maps a beta-lactamase family
// ("shv", "ctx-m") to its ResFinder prefix, so an allele ResFinder does not
// list can still be recognised as a beta-lactamase.
let nameIndex = null;

export function resfinderNameIndex() {
  if (nameIndex || !phenotypesMap) return nameIndex;
  const genes = new Map();
  const blaFamilies = new Map();
  for (const [key, val] of Object.entries(phenotypesMap)) {
    if (!Array.isArray(val)) continue; // accession keys hold single entries
    const k = key.toLowerCase();
    if (!genes.has(k)) genes.set(k, { name: key, entry: val[0] });
    const fam = key.match(/^bla([A-Za-z]+(?:-[A-Za-z]+)?)-?\d/);
    if (fam && !blaFamilies.has(fam[1].toLowerCase())) blaFamilies.set(fam[1].toLowerCase(), 'bla' + fam[1]);
  }
  nameIndex = { genes, blaFamilies };
  return nameIndex;
}

export function lookupGenePhenotype(geneNameOrAccession) {
  if (!phenotypesMap) return null;

  // Try exact match first
  if (phenotypesMap[geneNameOrAccession]) {
    const entry = phenotypesMap[geneNameOrAccession];
    return Array.isArray(entry) ? entry[0] : entry;
  }

  // Try matching by gene name (first part before underscore)
  const geneName = geneNameOrAccession.split('_')[0];
  if (geneName && phenotypesMap[geneName]) {
    const entry = phenotypesMap[geneName];
    return Array.isArray(entry) ? entry[0] : entry;
  }

  return null;
}

export function parseResFile(tsvContent) {
  const lines = tsvContent.trim().split('\n');
  if (lines.length < 1) return { headers: [], rows: [] };

  const headers = lines[0].replace(/^#/, '').split('\t').map(h => h.trim());
  const rows = [];

  for (let i = 1; i < lines.length; i++) {
    const line = lines[i].trim();
    if (!line) continue;

    const cols = line.split('\t');
    const row = {};
    headers.forEach((h, idx) => {
      row[h] = cols[idx] || '';
    });

    // Get the template/gene name
    const template = row['Template'] || '';
    if (template) {
      row._phenotype = lookupGenePhenotype(template);
    }

    rows.push(row);
  }

  return { headers, rows };
}
