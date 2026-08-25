# hv — hypervector species profiler

Metagenome → **one 32 kB signed vector** → deconvolution against an int8
reference matrix. Companion engine to the app's KMA gene detection: KMA
answers *which genes*, hv answers *which organisms, at what abundance*.

- **Design of record**: [`PLAN.md`](./PLAN.md) — theory (with the
  spectroscopy / radio-astronomy / CDMA parallels), formats, algorithms,
  milestones, validation protocol, honest limitations.
- **Core**: Rust, zero dependencies (`crates/hvcore`), compiled natively for
  validation and to `wasm32` for the browser — same code produces both sets
  of numbers.
- **CLI**: `hypv` (`crates/hvcli`).

## Quick start

```sh
cd hv
cargo build --release

# synthetic community with ground truth (6 species + a 3%-divergent novel organism)
./target/release/hypv simulate --out /tmp/hv/demo --n 6 --glen 300000 \
    --cov "8,5,4,3,2,0.4" --novel-rate 0.03 --novel-cov 1.5
./target/release/hypv build-db /tmp/hv/demo/manifest.tsv --out /tmp/hv/demo/db --spike
./target/release/hypv encode /tmp/hv/demo/reads_1.fq.gz --out /tmp/hv/demo/sample.hvf
./target/release/hypv profile /tmp/hv/demo/sample.hvf /tmp/hv/demo/db

# LOD sweep (detection floor vs coverage)
./target/release/hypv sweep /tmp/hv/sweep --dim 13   # D=8192  → floor ≈ 0.5×
./target/release/hypv sweep /tmp/hv/sweep --dim 15   # D=32768 → floor ≈ 0.1×
```

Real reference DB (WHO priority pathogens + commensals, ~25 genomes):

```sh
python3 scripts/fetch_genomes.py          # → genomes_ref/, manifest_ref.tsv
./target/release/hypv build-db manifest_ref.tsv --out ../hvprof/hv_patho_v1
./target/release/hypv encode ../SRR21386014_sub_1.fastq.gz ../SRR21386014_sub_2.fastq.gz --out /tmp/hv/saureus.hvf
./target/release/hypv profile /tmp/hv/saureus.hvf ../hvprof/hv_patho_v1
```

## Cross-validation vs sylph (real data)

Same 25-genome DB, same isolate, two independent methods:

| | HV profiler | sylph 0.9.0 |
|---|---|---|
| S. aureus JKD6159 | 100% of assignable mass, eff. cov 6.02, unexplained 47.9% | 100% abundance, ANI 98.53%, containment 0.633 |
| mystery sample ERR15682208 | 0 taxa, 100% uncharacterised | 0 taxa |

Consistency check: sylph's ANI 98.53% implies shared k-31 fraction 0.9853³¹ ≈
0.633 — exactly sylph's measured containment, and our fitted effective
coverage (base-coverage × window survival × identity × Q-survival) lands at
the observed 6.0×. The two observables HV reports (coverage + unexplained
mass) are linear combinations of the two sylph reports (coverage + ANI); v0
cannot separate them, M1's dual accumulator can.

Reference points from the published sylph (Nat Biotechnol 43:1348, 2025):
its real-data detection crossover vs MetaPhlAn4 is ~0.03× effective
coverage — the honest benchmark for our v0 floor of ~0.2× (D=32768,
σ₀-calibrated). The gap is structural (sylph keeps per-k-mer sketch data; we
compress into D buckets) and is the M1/M4 target. See PLAN.md §9.1 for the
full list of lessons imported from the paper, including PCR-duplicate
handling (their Illumina accuracy hinges on it — same physics as our §1.6
multiplicity term) and reference k-mer spacing rules.

## Measured v0 numbers (synthetic ground truth, σ₀-calibrated statistics)

| metric | value |
|---|---|
| detection floor, D=8192 (32 kB fp) | ~0.5× target k-mer coverage vs 10× background |
| detection floor, D=32768 (128 kB fp) | ~0.2× target k-mer coverage |
| floor vs DB size | flat from N=26 to N=2000 (EV penalty ≈ 1σ) |
| false discovery rate | 0 FP in ~26 000 null tests at N=2000 (q<0.05) |
| abundance accuracy (cov ≥ 0.5) | within ~10% (bounded by Q-masking survival ≈ 0.93) |
| false positives on novel organism | 0 (mass reported as "unexplained" instead) |
| DP release | ε=3 near-lossless; ε=1 ≈ halves z; ε≤0.3 kills weak signals |
| fingerprint aggregation | sum of 2 site fingerprints → correct pooled community |
| encode throughput (native) | ~90–170 MB/s |
| profile query | 0.03 s at N=2000 (sub-second at any planned DB size) |
| reference matrix | N × D int8 ≈ 0.3 MB for 26 pathogens |

Note: an earlier (pre-calibration) measurement quoted 0.1× at D=32768; that
statistic divided raw dots by the *depleted* noise scale, inflating z for
absent rows ~2–3×. The corrected floor is 0.2×. See PLAN.md §8/T6 and
`hypv bench fdr`.

## Web integration

`hvprof/hvprof.wasm` + `hvprof/hv_patho_v1.{hvd,tsv}` are served as static
files; `js/hv-worker.js` streams the user's FASTQ (never uploaded) through
the wasm encoder and posts the profile + fingerprint. Select
**HV pathogen profiler (β)** in the database dropdown on the main page.

## Tests

```sh
cargo test -p hvcore        # 17 unit + 7 scientific integration tests
```

Integration tests are the calibration claims of PLAN.md §8: exact recovery
(T4), vector linearity (T5), null calibration on novel background (T6),
LOD (T7), ANI strain/relative behaviour (T8), spike-in null.
