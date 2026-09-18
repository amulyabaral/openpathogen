# openpathogen

AMR and virulence gene detection in the browser. Reads never leave your machine.

Open [openpathogen.org](https://openpathogen.org), load FASTQ files (Illumina or Nanopore) or fetch a public run by accession, and run. No install, no account, no upload. Background, citations and usage notes are on the [About page](https://openpathogen.org/about.html).

**One click runs the pipeline:** optional quality control with [fastp](https://github.com/OpenGene/fastp) (compiled to wasm64 WebAssembly from source; see fastp/WASM_BUILD.md), then read mapping with [KMA](https://bitbucket.org/genomicepidemiology/kma) against any combination of ResFinder, CARD and VFDB, with adjustable identity and coverage thresholds. Results are tables of detected genes with a viewer per gene (coverage, consensus, alignment, variants), the fastp reports, a ZIP with the full KMA output, and a CABBAGE phenotype card.

All tools run locally as WebAssembly. Databases download on first use, are checked against pinned SHA-256 hashes (js/integrity.js) and are cached in the browser.

**Example isolate:** *S. aureus* USA300_TCH1516 ([SRR10341524](https://www.ebi.ac.uk/ena/browser/view/SRR10341524)), a 48 MB community-associated MRSA run downloaded from ENA and cached in your browser. Any other public run can be fetched by accession (SRR, ERR, DRR).

**[CABBAGE](https://doi.org/10.1093/nar/gkag780) phenotype associations.** After each run, a *Phenotype associations (CABBAGE)* card looks the detected resistance genes up in CABBAGE ([Dickens et al., NAR 2026](https://doi.org/10.1093/nar/gkag780)), which links 170,750 sequenced isolates to 1.7 million susceptibility test results. Genotypes are AMRFinderPlus calls; phenotypes use the 2025 CLSI and EUCAST breakpoints. For the sample's species, each gene is matched to its CABBAGE symbol, so *mecA* reads "97.7% of 1,791 *S. aureus* carriers were methicillin resistant", next to the rate in all sequenced *S. aureus*. Antibiotics that AMRFinderPlus links to the gene are listed first; other antibiotics are marked as co-occurrence. A species antibiogram (resistance rates from all CABBAGE records, filterable by country, year and source) gives context. The table is exported as `cabbage_predictions.csv` in the results ZIP. The snapshot downloads once (under 1 MB) and is cached. To browse the full database use the [EMBL-EBI AMR portal](https://www.ebi.ac.uk/amr). The snapshot is rebuilt with `scripts/fetch-cabbage.py` and `scripts/build-cabbage.py`; `scripts/pin-hashes.sh` updates the hash pins after any database change.

**Run locally:** any static file server at the repo root works, for example `python3 -m http.server 8000`. The bundled service worker adds the cross-origin isolation headers that WebAssembly threads need and reloads the page once on first visit.

Cite [KMA](https://doi.org/10.1186/s12859-018-2336-6), [fastp](https://doi.org/10.1093/bioinformatics/bty560) and [CABBAGE](https://doi.org/10.1093/nar/gkag780).
