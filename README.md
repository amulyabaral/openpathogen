# openpathogen

Find antimicrobial resistance genes, resistance mutations and virulence genes in your sequencing reads, in your browser. The reads never leave your computer.

Open [openpathogen.org](https://openpathogen.org), drop in your FASTQ files (Illumina or Nanopore) or fetch a public run by its accession, and press Run. That's it. There is nothing to install and no account to make. The [About page](https://openpathogen.org/about.html) has the background, the usage notes and the citations.

## How does openpathogen work?

openpathogen uses [fastp](https://github.com/OpenGene/fastp), [KMA](https://bitbucket.org/genomicepidemiology/kma) and [Deacon](https://github.com/bede/deacon) compiled from source to WebAssembly, and provides indexed database files. Because of this, your sequencing data is processed entirely client-side (in your own computer - which means there's no big, bulky servers anywhere!). No need to worry about accounts, queueing, or data compliance. Once you load the website and fetch the databases, you can even plug it off the internet and it'll still work, even on a smartphone.

A run is fastp for quality control if you want it, Deacon to take out human reads if you want that too, and then KMA mapping your reads against whichever of ResFinder, CARD and VFDB are ticked, with adjustable identity and coverage thresholds, plus PointFinder's genes for your species when it has them. Chrome, Edge and Firefox run on a 64-bit build of the tools. Safari, and every browser on iPhone and iPad (they all run on Safari's engine underneath), run on a 32-bit build of the same sources. The gene calls are the same either way, the 32-bit build just has at most 4 GB of memory to work with, and a good deal less on a phone. The build notes are in fastp/WASM_BUILD.md, kma/WASM_BUILD.md and deacon/WASM_BUILD.md.

### Databases

ResFinder and CARD for acquired resistance genes, VFDB for virulence factors, and [PointFinder](https://bitbucket.org/genomicepidemiology/pointfinder_db) for the known resistance mutations of ten species and genera (E. coli and Shigella, the *K. pneumoniae* complex, Salmonella, *S. aureus*, Campylobacter, both enterococci, *H. pylori*, *M. tuberculosis* and *N. gonorrhoeae*). Each index is downloaded once, from this site or from [Zenodo](https://zenodo.org/records/22102687), checked against a SHA-256 hash pinned in js/integrity.js and cached in the browser. A corrupted or altered download is simply rejected.

### Results by drug class

The results open with a report. Resistance genes are grouped by the drug classes and antibiotics they affect, from ResFinder's curated phenotype table, followed by the classes with no acquired gene and the virulence factors by VFDB category. You can print it or save it as a PDF. Below the report you get a gene table per database, a viewer per gene (coverage plot, consensus sequence, alignment and variants), the fastp and Deacon reports, and a ZIP with the full KMA output.

When the species is known (it is filled in for public runs, or you type it in), each gene also shows how often isolates of that species carrying it tested resistant in [CABBAGE](https://doi.org/10.1093/nar/gkag780), a database that links 170,750 sequenced isolates to 1.7 million susceptibility test results ([Dickens et al., NAR 2026](https://doi.org/10.1093/nar/gkag780)). So *mecA* in an *S. aureus* sample reads "In CABBAGE, *S. aureus* with mecA tested resistant to methicillin in 97.7% of 1,792 isolates". Rates are only shown for antibiotics the gene is known to affect, and only when at least 20 isolates were tested. They are in the results ZIP as `cabbage_predictions.csv`, and the whole snapshot is under 1 MB. If you want to dig through the full database, the [EMBL-EBI AMR portal](https://www.ebi.ac.uk/amr) is the place.

### Gene calls

KMA reports database templates, and a template is not always a gene in your sample. So before anything goes into the report, openpathogen looks at the reads behind each hit (KMA's own per-base output, nothing extra is run). A stretch of a gene with no reads, at a depth where that can't be chance, means that part of the gene is missing, and the gene is listed as partial with the missing part named. Reads that KMA divided between two sister alleles (which happens when your allele isn't in the database) are put back together into one gene, and the closest allele is named with the number of differences. A stop codon or frameshift that the reads carry takes a gene out of the counts. The defaults, 90% identity and 90% coverage, call the genes ResFinder calls, and you can move both sliders. The details are in the notes on the [About page](https://openpathogen.org/about.html).

How well does that work? 34 sequenced isolates of 32 species have genome assemblies in NCBI, most of them complete, so their genes are known. Run through openpathogen at 40× of their reads, 119 of the 120 intact genes in those genomes were called, the allele was named right for 114 of 116, and duplicate calls went from 33 (KMA's numbers alone) to 1. Most of the calls with no gene in the genome turned out to be plasmid genes at 60 to 860× depth that the assemblies themselves left out, and the two that were really contamination, at under 2×, are marked in the report. PointFinder's mutations agreed with NCBI's AMRFinderPlus on the same genomes wherever PointFinder's table holds the mutation.

### Point mutations

Give the species and, if PointFinder covers it, the reads are also mapped to PointFinder's wild-type genes for that species (gyrA, parC, the 23S rRNA, promoters and so on), with the same KMA command. Every position in PointFinder's mutation table is then read straight from KMA's per-base counts, so the report can say *gyrA p.S83L* and *parC p.S80I* next to the genes, with PointFinder's antibiotics and notes and the CABBAGE rates of isolates carrying the same mutation. A few things are done differently from ResFinder on purpose: a position with no reads is reported as not read instead of as a deletion, the copies of an rRNA gene are read together (with the share of reads when only some carry the change), an insertion is found wherever the reads place it in a repeat, and a mutation that PointFinder only counts alongside another one is listed apart when the other one isn't there. Every position as read, with the reason when it couldn't be, is in the ZIP. scripts/build-pointfinder.py builds the indexes from pointfinder_db and checks every row against its gene first. Doing that against real genomes turned up a surprise: PointFinder's Klebsiella *gyrA* and *parC* are *E. coli* genes, so *K. pneumoniae* reads never reach the 80% identity a gene needs and its quinolone mutations go unread (in ResFinder too). openpathogen reads them on the *K. pneumoniae* type strain's genes instead, and since PointFinder marks its whole Klebsiella table as under development, a Klebsiella row only counts when NCBI's AMRFinderPlus curates the same mutation.

### Human reads

Swabs, tissue and other clinical samples carry human DNA. Tick "Remove human reads" and Deacon 0.18.0 (its deacon-wasm build, unmodified) takes out the reads that match the human pangenome before the databases run, on your device like everything else. It uses Deacon's panhuman-1 index (k=31, w=61), frozen into a 239 MB binary fuse filter that downloads once and is hash-pinned like the gene indexes. The report tells you how many reads were removed, and the remaining reads can be downloaded as gzipped FASTQ. In testing, the browser build gave byte-identical output to native Deacon, removing 98.5% of a human WGS run. The details are in deacon/WASM_BUILD.md. Deacon catches most human reads, but not every single one, so check the rules that apply to your samples before you share the leftover reads.

### hAMRonization

Every gene hit (and every PointFinder mutation, as a variant row) can also be downloaded in the [PHA4GE hAMRonization](https://github.com/pha4ge/hAMRonization) format (TSV from the results page, TSV and JSON in the ZIP), so you can put openpathogen results together with the output of other AMR tools, for example with `hamronize summarize`. The files match what the hAMRonization 1.3.1 package itself writes for the same records. The field mapping is documented in js/hamronize.js.

## Frequently asked questions

### Is my sequencing data uploaded anywhere?

No. FASTQ files are read into browser memory, processed by WebAssembly and discarded when you close the tab. The only downloads are the app, the databases and, if you ask for one, a public run from ENA. The [privacy section](https://openpathogen.org/about.html#privacy) of the About page lists every request the site makes, so you can check it against your browser's network inspector.

### Which inputs are supported?

Illumina paired-end, single-end and Oxford Nanopore reads in FASTQ, plain or gzipped. You can also fetch any public run by its SRA or ENA accession (SRR, ERR or DRR). Files up to 1 GB each are handled in memory.

### Is there an example I can try?

Yes. The example is *S. aureus* USA300_TCH1516 ([SRR10341524](https://www.ebi.ac.uk/ena/browser/view/SRR10341524)), a 48 MB MiSeq run of the community-associated MRSA reference strain, downloaded from ENA and cached in your browser.

### Can I run it on my own machine?

Yes. Clone the repository and serve it with any static file server, for example:

```
python3 -m http.server 8000
```

The bundled service worker adds the cross-origin isolation headers that WebAssembly threads need, and reloads the page once on your first visit. One thing is not in the repository: the human index for Deacon is over GitHub's file size limit, so `deacon/build.sh --index` downloads Deacon's index and freezes it into `databases/deacon/` for you (the result matches the pinned hash). If you change a database, `scripts/pin-hashes.sh` updates the hash pins, the CABBAGE snapshot is rebuilt with `scripts/fetch-cabbage.py` and `scripts/build-cabbage.py`, and the PointFinder files with `scripts/build-pointfinder.py <pointfinder_db checkout>`.

### How do I cite openpathogen?

A short communication is in preparation. Until then, cite [KMA](https://doi.org/10.1186/s12859-018-2336-6), [fastp](https://doi.org/10.1093/bioinformatics/bty560) and [CABBAGE](https://doi.org/10.1093/nar/gkag780), [Deacon](https://doi.org/10.1101/2025.06.09.658732) if you removed human reads, [PointFinder](https://doi.org/10.1093/jac/dkx217) if it ran, and the databases you used. Database citations and component licences are in [THIRD_PARTY_NOTICES.md](THIRD_PARTY_NOTICES.md). openpathogen itself is MIT licensed.
