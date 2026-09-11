# openpathogen

AMR and virulence gene detection in the browser. Reads never leave your machine.

Open [openpathogen.org](https://openpathogen.org), load FASTQ (Illumina or Nanopore), and run. No install, no account, no upload. Background, citations and usage notes are on the [About page](https://openpathogen.org/about.html).

**Simple mode (default)** runs the whole pipeline in one click: quality control with [fastp](https://github.com/OpenGene/fastp) (compiled to wasm64 WebAssembly from source; see fastp/WASM_BUILD.md) → resistance genes (ResFinder + CARD) → virulence factors (VFDB) → a plain-language summary of what was found and what it means, with explicit "genotype ≠ phenotype" caveats. **Advanced mode** exposes each database individually with full KMA tables and tunable thresholds.

All tools run locally as WebAssembly; databases download on first use and cache in the browser.

**Example isolate** is *S. aureus* USA300_TCH1516 ([SRR10341524](https://www.ebi.ac.uk/ena/browser/view/SRR10341524)): a 48 MB community-associated MRSA run downloaded from ENA and cached in your browser. Any other public run can be fetched by accession (SRR…, ERR…, DRR…).

**[CABBAGE](https://doi.org/10.1093/nar/gkag780) phenotype predictions.** After every run, a *Phenotype prediction (CABBAGE)* card interprets the detected genes against the comprehensive AMR genotype–phenotype database ([Dickens et al., NAR 2026](https://doi.org/10.1093/nar/gkag780)): for the sample's species, each gene is looked up among ~165,000 isolates that have both a genome and an antibiogram, so *mecA* becomes "97.7% of 1,791 *S. aureus* isolates carrying mecA were methicillin-resistant" — with the species' background rate alongside for comparison. A species antibiogram (per-antibiotic resistance rates from all 1.7M CABBAGE AST records, refinable by country/year/isolation source) provides context. Predictions are exported as `cabbage_predictions.csv` in the results ZIP. The association table (~0.3 MB) downloads once and caches; to browse the full database use the [EMBL-EBI AMR portal](https://www.ebi.ac.uk/amr). Snapshots rebuild with `scripts/fetch-cabbage.py` + `scripts/build-cabbage.py`.

Cite [KMA](https://doi.org/10.1186/s12859-018-2336-6), [fastp](https://doi.org/10.1093/bioinformatics/bty560) and [CABBAGE](https://doi.org/10.1093/nar/gkag780).
