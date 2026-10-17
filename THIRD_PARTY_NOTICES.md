# Third-party notices

The openpathogen application code is released under the MIT license (see
LICENSE). The components below are distributed with or used by the app and
remain under their own licenses and terms — the MIT license does not apply
to them. Please cite the publications below if you use openpathogen results
in published work.

## Tools (WebAssembly ports)

### fastp — MIT license
- Source: https://github.com/OpenGene/fastp
- Quality filtering and adapter trimming.
- Cite: Chen S., Zhou Y., Chen Y., Gu J. (2018). fastp: an ultra-fast
  all-in-one FASTQ preprocessor. *Bioinformatics* 34(17): i884–i890.

### KMA 1.5.1 — Apache License 2.0
- Source: https://bitbucket.org/genomicepidemiology/kma
- K-mer-based read mapping for gene detection.
- Cite: Clausen P.T.L.C., Aarestrup F.M., Lund O. (2018). Rapid and precise
  alignment of raw reads against redundant databases with KMA. *BMC
  Bioinformatics* 19: 307.
- Full Apache 2.0 license text: https://www.apache.org/licenses/LICENSE-2.0

### Deacon 0.18.0 (deacon-wasm) — MIT license
- Source: https://github.com/bede/deacon (unmodified; built by
  deacon/build.sh, see deacon/WASM_BUILD.md)
- Human read removal by minimizer comparison with a host index. Its
  minimizers come from SimdMinimizers (packed-seq, simd-minimizers; MIT).
- Cite: Constantinides B., Lees J., Crook D.W. (2025). Deacon: fast sequence
  filtering and contaminant depletion. *bioRxiv* 2025.06.09.658732.
  Also: Groot Koerkamp R., Martayan I. (2025). SimdMinimizers: computing
  random minimizers, fast. *bioRxiv* 2025.01.27.634998.

## Databases

### Deacon panhuman-1 human index (k=31, w=61) — CC0 1.0
- Source: https://github.com/bede/deacon-indexes (panhuman-1 workflow;
  the w=15 build is Zenodo record 17288185). Built from HPRC year 1
  assemblies, CHM13v2.0 and GRCh38.p14, minus minimizers shared with
  FDA-ARGOS bacteria and RefSeq viruses. openpathogen ships it frozen into
  a 16-bit binary fuse filter (deacon index freeze -b 16).

### ResFinder 2.6.0 — Apache License 2.0
- Source: https://bitbucket.org/genomicepidemiology/resfinder_db
- Acquired antimicrobial resistance genes.
- Cite: Zankari E., et al. (2012). Identification of acquired antimicrobial
  resistance genes. *Journal of Antimicrobial Chemotherapy* 67(11): 2640–2644.

### CARD 4.0.1 (Protein Homolog model) — CARD Terms of Use
- Source: https://card.mcmaster.ca
- Antimicrobial resistance determinants. CARD data are freely available for
  personal and public **non-commercial use only**, under the terms of use at
  https://card.mcmaster.ca/about; commercial use requires a written license
  from McMaster University.
- Cite: Alcock B.P., et al. (2023). CARD 2023: expanded curation, support
  for machine learning, and resistome prediction at the Comprehensive
  Antibiotic Resistance Database. *Nucleic Acids Research* 51(D1): D690–D699.

### VFDB core set A — CC BY-NC 4.0
- Source: https://www.mgc.ac.cn/VFs/
- Virulence factors. Licensed under the Creative Commons
  Attribution-NonCommercial 4.0 International license; commercial use
  requires a separate license from the VFDB team.
- Cite: Chen L., Yang J., Yu J., Yao Z., Sun L., Shen Y., Jin Q. (2005).
  VFDB: a reference database for bacterial virulence factors. *Nucleic Acids
  Research* 33(S1): D325–D328.

The KMA indexes built from these databases are derived works and carry the
same terms, including the copies served from the Zenodo index record
(https://zenodo.org/records/22102687). The record's own licence field does
not extend CARD's or VFDB's terms; for commercial use, obtain a licence from
McMaster University and the VFDB team.

## Formats

### hAMRonization — GNU LGPL v3
- Source: https://github.com/pha4ge/hAMRonization
- openpathogen writes its gene hits in the hAMRonization 1.3.1 format (the
  PHA4GE AMR gene and variant specification). No hAMRonization code is
  included; js/hamronize.js is an independent implementation of the format.

## Vendored JavaScript

### fflate — MIT license
- Source: https://github.com/101arrowz/fflate
- Decompression of downloaded database archives.

### plotly.js 1.2.0 — MIT license
- Source: https://github.com/plotly/plotly.js (the minified build that
  fastp's HTML report loads from opengene.org, served here from js/vendor/
  instead so the report needs nothing beyond this site)
- Copyright 2012-2015 Plotly, Inc. Draws the plots in the fastp quality report.

### CABBAGE (AMR genotype–phenotype database) — see data-source terms
- Source: https://www.ebi.ac.uk/amr (EMBL-EBI Antimicrobial Resistance
  portal); processing scripts https://github.com/Leonardini/CABBAGE
- Experimental antimicrobial susceptibility phenotypes linked to sequenced
  isolates, plus NCBI AMRFinderPlus in silico genotype predictions. The app
  ships a compact snapshot pinned to a portal release (release id in
  databases/cabbage/manifest.json). Data originates from BV-BRC/PATRIC,
  NCBI Pathogen Detection, PubMLST, EUCAST and literature curation; please
  respect the terms of the portal and the original data submitters, and
  use the live portal for the current release.
- Cite: Dickens E., et al. (2026). A comprehensive AMR genotype–phenotype
  database (CABBAGE). *Nucleic Acids Research* 54(16): gkag780.
