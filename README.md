# openpathogen

AMR and virulence gene detection in the browser. Reads never leave your machine.

Open [openpathogen.org](https://openpathogen.org), load FASTQ (Illumina or Nanopore), and run. No install, no account, no upload.

**Simple mode (default)** runs the whole pipeline in one click: quality control with [fastp](https://github.com/OpenGene/fastp) (via [biowasm](https://biowasm.com)) → resistance genes (ResFinder + CARD) → virulence factors (VFDB) → a plain-language summary of what was found and what it means, with explicit "genotype ≠ phenotype" caveats. **Advanced mode** exposes each database individually with full KMA tables and tunable thresholds.

All tools run locally as WebAssembly; databases download on first use and cache in the browser.

**Example isolate** is *S. aureus* JKD6159 ([SRR21386014](https://www.ebi.ac.uk/ena/browser/view/SRR21386014)): a ~50 MB subsample ships with the app; the full run can be fetched from ENA.

Cite [KMA](https://doi.org/10.1186/s12859-018-2336-6) and [fastp](https://doi.org/10.1093/bioinformatics/bty560).
