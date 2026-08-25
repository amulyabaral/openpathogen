# HyperVector profiler (HV) — full implementation plan

*A species-profiling engine for OpenPathogen: metagenome → one fixed-length
hypervector → deconvolution against a reference matrix. Companion to the
existing KMA gene-detection engine: KMA answers "which genes", HV answers
"which organisms, at what abundance".*

**Status of this document**: design-of-record. Sections 1–4 are settled
theory and architecture; §5–§7 are the build plan; §8 is the validation
protocol that gates each milestone; §9 records risks and honest limitations.

---

## 0. One-paragraph summary

Stream a FASTQ once; for each canonical 31-mer that passes an i.i.d. density
filter, add ±1 into one of D=8192 buckets of an Int32 accumulator chosen by
hash bits. The result — a 32 KB vector — is, in expectation, a **random
signed linear projection of the sample's k-mer count spectrum**, so it obeys

```
v(sample) ≈ Σᵢ μᵢ · v(genomeᵢ) + noise        (μᵢ = coverage of genome i)
```

Profiling = sparse non-negative recovery of μ against a matrix of reference
projections, followed by per-taxon significance tests in projection space,
iterative depletion (subtract what is explained, re-test the residual), and a
bounded-size fingerprint the user can keep, share, or aggregate. The browser
is the natural client because the encoder is sequential, branchless,
O(1)-memory streaming — the exact shape WASM is best at — and the reference
matrix is int8-quantized, appendable without rebuilds.

Positioning vs prior art (must be cited and benchmarked, not skipped):
Quikr/WGSQuikr (2013/14) invented compressed-sensing abundance estimation
over short-k-mer spectra; MetaPalette (2018) did NNLS over k-mer palettes;
sourmash gather (2024) does greedy subtractive sketch decomposition; sylph
(Nat Biotech 2024) is the accuracy/speed SOTA with zero-inflated-Poisson
sketch statistics; HyperGen (Bioinformatics 2024) established HDC genome
sketches for ANI. **Our claim is not "first deconvolution"** — it is: *the
first profiler whose entire state is a bounded, fixed-length, strand-symmetric
signed projection that supports (a) joint NNLS deconvolution, (b) calibrated
detection statistics, (c) exact aggregation across samples, (d) sub-second
queries against an append-only int8 matrix, in O(1) memory — deployable as a
static web page.* Novelty is in the combination and the access model, and the
paper must be written that way.

---

## 1. Theory — and what each borrowed idea buys us

HV sits at the junction of five fields. Each parallel below has a **concrete
design consequence** (marked ⇒), which is why it is listed.

### 1.1 Compressed sensing / sparse recovery (mathematics)

The encoder is a random linear operator S (one ±1 entry per selected k-mer,
i.i.d. Bernoulli(ρ) row selection) applied to the k-mer count vector c:
`v = S·c`. Since `c = Σᵢ μᵢ·1[Gᵢ]` in expectation, `v = Σᵢ μᵢ·S·1[Gᵢ]`.
Communities are sparse (10²–10³ taxa vs N=10³–10⁴ references), so recovery of
μ is a standard non-negative sparse least-squares problem. RIP-style heuristics
for Bernoulli ±1 projections require roughly `D ≳ C·s·log(N/s)`; with s≈300,
N≈4000, C≈3–5 → D≈3–8k. **D=8192 default, 32768 "high-sensitivity" mode.**

Mutual incoherence of the reference matrix is bounded below by collision
noise: unrelated genomes have empirical column correlation ≈ D^(-1/2)
(≈0.011 at D=8192). Within-genus columns correlate by shared k-mer fraction
(ANI ~90% → k=31 sharing of a few %–tens of %) — the dictionary is a
**union of near-orthogonal blocks with coherent sub-blocks** (species
complexes). ⇒ hierarchical reporting: solve jointly, but collapse
identifiability-limited species complexes to genus when their post-fit
coefficients are unstable (§4.6). Donoho–Tanner-style phase-transition plots
(recoverable sparsity vs D) are a planned validation figure.

### 1.2 Spectroscopy & chemometrics (chemistry)

A mixture's absorbance spectrum is the concentration-weighted sum of pure
spectra (Beer–Lambert law); ours is the coverage-weighted sum of genome
projections. A century of mixture spectroscopy already solved our problems:

- **Classical least squares with non-negativity and closure** ⇒ NNLS, and
  (v1.1) a closure constraint Σxᵢ ≤ total mass.
- **MCR-ALS (multivariate curve resolution by alternating least squares)**
  ⇒ our depletion loop: fit → subtract explained part → re-test residual →
  admit new candidates → refit (≤3 rounds).
- **Hadamard spectrometry / Fellgett's multiplex advantage**: measuring with
  a ±1 encoded mask instead of a 0/1 mask converts systematic overlap *bias*
  into zero-mean *variance*. ⇒ **signs are mandatory** — an unsigned
  CountMin-style sketch has a positive collision bias that scales with
  background mass and would make abundance estimates systematically
  over-positive; signed accumulation cancels it exactly in expectation.
- **Isotope dilution / standard addition** ⇒ a synthetic spike-in genome row
  (random 100 kb sequence, ships in the DB) at known copy number verifies
  encoder health and calibrates coverage per run — the ERCC-spike trick.

### 1.3 Radio astronomy (physics)

Aperture synthesis images a sparse sky from correlation "visibilities"
(baselines); Högbom CLEAN is literally orthogonal matching pursuit with
successive subtraction, invented in 1974 for radio images. Our buckets are
baselines, taxa are point sources, OMP+depletion is CLEAN. Consequences:
- **CLEAN's lesson: source confusion** ⇒ when two sources are within a
  synthesised beam (two strains within a projection's coherence), report the
  *blend* at the right total flux, not a spurious split — same as our
  genus-collapse rule.
- **Radiometer equation / matched filtering** ⇒ detection statistic
  zᵢ = ⟨v, aᵢ⟩ / (σ̂·‖aᵢ‖) with σ̂ a robust per-bucket noise scale
  (MAD-based); SNR grows as √(coherent events), i.e. with √coverage.

### 1.4 CDMA multiuser detection (engineering)

Deconvolving superposed users with known spreading codes = our problem with
genomes as users and projections as codes. OMP is successive interference
cancellation; the **near–far problem** (a strong user masking weak ones) is
exactly "dominant commensal masks low-abundance pathogen". ⇒ the depletion
loop is not a demo gimmick, it is the standard near–far remedy and must be
on-by-default.

### 1.5 Extreme value theory / statistical mechanics (mathematics)

The screening stage takes the max z over N tests. Max of N ≈ 4·10³ Gaussian
nulls concentrates near √(2 ln N) ≈ 4.3σ with Gumbel fluctuations — random
energy model mathematics. ⇒ screening thresholds and FDR control must be
calibrated to the extreme-value distribution, not a single-test Gaussian
tail; we implement Benjamini–Hochberg on permutation/exchangeable-null
p-values plus a Gumbel sanity check (§4.5, §8 test T6).

### 1.6 The multiplicity problem (the subtle one — read twice)

References are built from **distinct** k-mers (one copy per genome); reads
arrive as **occurrences**. If the encoder naively accumulates occurrences,
background variance carries an E[count²] = μ+μ² term: at 30× coverage the
noise floor is ~√31 ≈ 5.6× higher than at 1×, and *deeper sequencing makes
detection worse*. Two remedies, both implemented:

- **Distinct-mode**: a ~2^27–2^29-bit Bloom "seen" filter (17–68 MB) makes
  the encoder count each distinct selected k-mer once, matching the reference
  construction; the floor then *improves* with depth (Poisson thinning of
  distinct k-mers) until saturation. Trade-off: distinct-count signal
  saturates at high coverage (1−e^(−μ)) — fine for detection, wrong for
  abundance.
- **Occurrence-mode with count-aware variance**: keep occurrences (linear in
  μ — right for abundance), and model Var = Σ(μ+μ²) in the z-test
  (sylph's zero-inflated-Poisson insight, ported to projection space).

⇒ production default: **dual accumulator** — occurrence vector (NNLS /
abundance) + distinct vector (detection z-tests), one Bloom filter shared.
v0 (this repo, first milestone) ships occurrence-mode with count-aware
variance to keep memory at ~50 MB; distinct-mode lands in v0.2 and is
expected to buy ~3–5× in LOD at high coverage (to be measured, §8).

**Measured 2026-08 (prototype, HashSet in lieu of Bloom — `hypv bench
distinct`)**: at 30× per-genome background, occurrence-mode detects NOTHING
even at 2× target coverage (the μ² tax is fatal), while distinct-mode at
D=8192 detects to 0.5× (z=11.5) and at D=32768 to **0.25×** (z=6.95), with
μ̂ = −ln(1−ĉ) recovering true coverage (c=0.184 → 0.20 vs true 0.25).
The two factors multiply as predicted (√μ_bg from dedup × √2 from D).
M1 architecture is validated; the Bloom filter is now an engineering task,
not a research risk.

### 1.6.2 The rich-background wall and the hybrid answer (validated 2026-08)

Acid test: 200-genome community × 2× coverage (D_bg ≈ 3×10⁶ distinct
selected k-mers — gut-metagenome-like; ≈ 91 background k-mers/bucket at
D=32768). Floors for a target genome, σ₀-calibrated, kmax > community size:

| D | occurrence floor | distinct floor |
|---|---|---|
| 32 768 | >2× | 1× |
| 131 072 | 1× | **0.25×** (z=5.75) |
| 262 144 | 0.5× | 0.25× (z=9.0; 0.1× still below floor) |

Empirical law: z ∝ μ·√D in the collision regime; **selection density ρ
cancels** (row mass and background mass both scale with ρ). Reaching sylph's
practical 0.03× floor by D alone needs ~10⁷ buckets (≈ 2–14 MB fingerprint):
the projection-only path stops making sense there. This is the
information-theoretic wall: a bounded destructive projection pays collision
noise √(D_bg/D) that identity-preserving sketches do not; the ~100× state
ratio (tens of kB vs sylph's 1–2 MB sketch) buys a ~10× floor ratio, in
both directions.

**Found the hard way**: with 200 true community members, the top-K admission
(kmax=96) silently evicted true positives — the support filled with the
luckiest background rows and the "floor" was an artifact. kmax must scale
with community richness (or admission must be pure-threshold). Fixed in the
bench (`kmax = n + 8`); production default raised in M1.

**The hybrid architecture that IS good enough** (M1/M2 design of record):

1. **Distinct-mode dual accumulator** (Bloom seen-set) — validated above.
2. **Nested dual-width projection**: one accumulation pass, two bucket
   widths (screen D_s=16384 = top bits, confirm D_c=131072). Screen matrix
   genus-collapsed (≈300 rows × D_s — coherent sums detect any member;
   analytic √(L_genus/L_species) ≈ 2–4× discovery sensitivity), confirm
   rows for candidate genera streamed from OPFS on demand.
3. **Exact confirm stage** (M2): for candidates + a pre-registered
   watchlist (WHO BPPL), keep per-genome fine-density (ρ=0.005, nested
   subset of our selection ✓) k-mer sets in memory (64 genomes × 25k × 8 B
   ≈ 13 MB) and count exact containment on a second streaming pass —
   sylph's own statistic (ZIP, ratio-of-multiplicity λ̂, ANI). Watchlist
   detection then runs at full sylph parity (~0.03×) with the projection
   providing discovery ≥0.25×, fingerprint/aggregation/DP always.
4. **Skellam/GLS weights** on the confirm-stage z (1.1–1.3×).
5. **Pooling**: federation z ∝ √M across sites.

Expected system floors (gut metagenome): watchlist 0.03× (sylph parity),
de-novo discovery 0.25× at 256 kB fingerprint (D_c=131072 int16-ish),
16 kB federation fingerprint for sharing at 0.5–1×. "Good enough" is
redefined honestly: parity where the clinical question lives (known
pathogens), projection-grade discovery, and bounded shareable state —
rather than pretending a 32 kB sketch beats a 2 MB one at its own game.

### 1.6.1 Further optimizations (brainstormed + triaged 2026-08)

- **Skellam/GLS variance weighting** — each bucket is a difference of two
  Poissons with bucket rate λ_b ≈ Σⱼ x̂ⱼ|aⱼ[b]|; NNLS assumes homoscedastic
  Gaussian. IRLS with w_b = 1/λ̂_b is the true MLE. Expected 1.1–1.3× in z;
  cheap once the fit exists. M1.
- **Early-exit encoding (Wald SPRT)** — the accumulator is an online
  statistic; screen periodically and stop when composition is determined.
  Measured on the real S. aureus isolate: **2% of the file already gives a
  confident call** (z=16.5, correct 47.3% unexplained); 10% gives z=27.8.
  Browser UX win: provisional results in <1 s, refine on demand.
- **Two-resolution screen-then-confirm** — encode dual accumulators
  (D=8192 + D=65536; 288 KB total) at encode time; screen permissively
  (z>1.5) against the full small matrix, then load only the ~100 candidate
  rows of the big matrix on demand from OPFS for confirmation. The floor
  becomes the confirmation floor (~2× better than D=32768) while resident
  memory stays ~32 MB — the answer to the WASM memory ceiling.
- **Pooling gain (federation)** — sum of M fingerprints: signal M·μ, noise
  √M ⇒ z ∝ √M. 100 participating sites = 10× floor improvement for shared
  targets; no per-sample optimization competes with this.
- Rejected with reasons: antithetic variates (<√2 for 2× events);
  multi-sketch averaging (identical to one larger D — only median-robustness
  vs duplicate outliers, which dedup fixes better); deterministic/equiangular
  measurement matrices (matrix is induced by hashing genomes, not chosen);
  lower k (8% more events, large specificity cost); quality-weighted
  accumulation (masking already captures most signal).

### 1.7 Sequencing-error chemistry

Per-base error ε creates ~k·ε novel k-mers per error (each poisons a 31-wide
window). At ε=10⁻³, a 3 Gbp run carries ~10⁸ error k-mer occurrences —
comparable to the real distinct signal. They are unique, so they land as
zero-mean random noise (no coherent bias toward any reference) but inflate
σ̂. ⇒ quality masking (Q<20 → N, skip affected k-mers) is **on by default**
and removes most of this at negligible cost. ε is also estimated per-run
from the "k-mers seen exactly once" tail where a Bloom+counter is available
(v0.2), feeding the variance model.

### 1.8 Differential privacy (only promise what holds)

The fingerprint is information-destroying (10⁷–10⁸ k-mers → 8192 numbers,
non-invertible in any practical sense) but *composition-revealing by design*.
Correct claim: **human-DNA-private, not composition-private** (a human row
can even be depleted). DP release is an optional mode: Gaussian noise on the
int32 vector calibrated to per-read sensitivity (one read contributes ≤ ⌈ρ·(r−k+1)⌉
signed events ⇒ L2 sensitivity ≈ √(ρ·r)); the ε-vs-LOD trade-off is measured
empirically (§8 F4) and reported honestly. Linearity gives one free gift:
**a server can sum member fingerprints and deconvolve the aggregate** —
population prevalence surveillance from additive 32 KB blobs.

---

## 2. Notation, parameters, defaults

| symbol | meaning | default | range |
|---|---|---|---|
| k | k-mer length | 31 | fixed 31 v1 |
| D | projection dimension | 8192 | 4096/8192/16384/32768 |
| ρ | selection density (i.i.d. mod-hash) | 0.05 | 0.02–0.25 |
| seed | hash seed (bakes DB & fingerprints) | fixed constant v1 | any u64 |
| Qmin | Phred mask threshold | 20 | 0–30 |
| N | reference rows | ~40 (demo DB) → 4k+ | — |
| K | screen→NNLS candidate cap | 96 | 32–256 |
| R | depletion rounds | 3 | 1–5 |
| z_min | per-taxon report threshold | 5.0 | post-BH q<0.05 governs |
| row store | int8 + per-row f32 scale | — | — |

Mod-hash (i.i.d. Bernoulli) selection — **not** minimizers: minimizer density
is GC-skewed and window-correlated, which violates the i.i.d. sampling the
§1.1 theory wants; mod-hash is uniform, cheaper (no window state), and
error-robust (a base error cannot re-select a neighboring window). Cost:
finalize a hash for every k-mer (~1 splitmix64 per k-mer, ~20–40 s/Gbp
native, ~2–3× in WASM — encode stays I/O-bound).

---

## 3. Architecture

One Rust core, three faces. **The numbers in the paper and the numbers in
the browser come from the same code** (the repo's existing native↔WASM
concordance culture, extended to HV).

```
hv/
  Cargo.toml                 # workspace
  crates/hvcore/             # pure algorithms, zero deps, no I/O
  crates/hvcli/              # bin `hypv`: build-db, encode, profile,
                             #   simulate, sweep, bench  (flate2 only)
  crates/hvwasm/             # cdylib → hvprof/hvprof.wasm (no-bindgen,
                             #   hand-rolled exports + tiny bump allocator)
  PLAN.md                    # this file
  python/                    # (later) figures, CAMI glue, baselines
hvprof/                      # browser-deployable artifacts (like kma/)
  hvprof.wasm                # built wasm
  hv_patho_v1.hvd            # int8 matrix DB + sidecar
  hv_patho_v1.tsv            # row metadata: id, name, genus, len, L_sel
js/hv-runtime.js             # wasm loader + profiling API (main thread)
js/hv-worker.js              # streaming encoder worker (module worker)
```

Core module map (`hvcore/src/`):

| module | contents |
|---|---|
| `hash.rs` | canonical rolling hash: 2-bit Rabin–Karp forward `h_f` and reverse `h_r` (odd multipliers), strand-symmetric `u = splitmix64(h_f ^ h_r)`; `derive(u) → (selected, bucket, sign)` via a second mix on disjoint bits |
| `encode.rs` | FASTQ/FASTA byte-stream state machine, Q-mask, occurrence & distinct accumulators, (v0.2) Bloom seen-filter |
| `db.rs` | DB build from FASTA (per-genome distinct selected set → f64 row → int8 + scale), `.hvd` read/write, append row |
| `stats.rs` | robust σ̂ (MAD), z-tests, count-aware variance, BH-FDR, Gumbel screen bound |
| `solve.rs` | screening pass, OMP-lite candidate expansion, Lawson–Hanson NNLS on Gram, depletion loop, genus collapse, spike-in QC |
| `fp.rs` | fingerprint (.hvf) v1: header + i16-quantized vector + counters; optional Gaussian noise (DP) |
| `sim.rs` | seeded synthetic genomes, ANI-controlled mutants, Illumina-style read simulator (Poisson coverage, per-base errors, Q strings) |

**Format `.hvd` (v1)** — little-endian: magic `HVP1`, u32 version, u32 k,
u32 D, f64 ρ, u64 seed, u32 N, N×D i8 rows, N×f32 scales, N×u32 L_sel,
N×u32 genome_len. Sidecar TSV: `row_idx, id, name, genus, species, flags`.
**Format `.hvf` (v1)**: magic `HVF1`, same param header, u64 n_bases,
u64 n_events, f32 scale, D×i16. ~16.5 KB at D=8192. Parameters live in the
file; mismatched DB↔fingerprint params are a hard error (loud, not silent).

**Browser pipeline** (`js/hv-worker.js`): `File.stream()` →
`DecompressionStream('gzip')` when magic is 1f 8b → 4 MiB chunks →
`hv_encoder_feed(ptr, len)` → progress posts (bytes, events, MB/s) → on EOF
copy accumulator out → fetch/cached `.hvd` (IndexedDB via existing `db.js`)
→ `hv_profile()` → results array + `.hvf` fingerprint → main thread renders.
No SharedArrayBuffer, no COOP/COEP requirement, single-threaded worker;
WebGPU is explicitly **not** used (query is sub-second; the bottleneck is
decode+encode, where GPU helps nothing).

---

## 4. Algorithms (pseudocode of record)

### 4.1 Encoder (per k-mer, occurrence mode)

```
roll hf, hr over 2-bit codes; if any base non-ACGT or Q<Qmin: window invalid for k-1 bases
u = splitmix64(hf ^ hr)                 # strand-symmetric, avalanched
if u < T (= ρ·2^64):                    # i.i.d. selection
    w = splitmix64(u ^ KAPPA)           # independent bits
    acc[(w >> 51) & (D-1)] += (w >> 50) & 1 ? 1 : -1
    n_events += 1
```

Invariant tests: hash(seq) == hash(revcomp(seq)) (T1); selection/bucket/sign
uniformity (χ², T2); int32 saturation impossible at practical depths
(expected |acc[b]| ≤ ~10⁴; guard with saturating add anyway).

### 4.2 DB row build (per reference genome)

same encoder over FASTA, but each **distinct** selected hash contributes once
(per-genome HashSet<u64>); row normalised to f64, then `scale=127/max|v|`,
int8 store; record L_sel and genome length. Human/commensal/viral rows are
ordinary rows (flags differ) — depletion needs no special casing.

### 4.3 Screening

`dot_i = ⟨v_f32, row_i⟩·scale_i` (batched, int8×f32); `z_i = dot_i/(σ̂‖a_i‖₂)`
with σ̂ from MAD of v buckets. Keep top-K by z **plus** spike-in row always.
Gumbel bound: warn if max null z expected under H0 (√(2 ln N_eff)) is within
1.5 of observed top z (N_eff discounts row correlations).

### 4.4 Joint solve

Candidates C (|C| ≤ K): Lawson–Hanson NNLS on the Gram system
(Aᶜᵀ Aᶜ) x = Aᶜᵀ v, columns pre-scaled to unit norm. Output x̂ ≥ 0.
Unit tests pin exact recovery on synthetic incoherent mixtures (T4).

### 4.5 Depletion loop (rounds r = 1..R)

```
resid = v − Σ x̂ᵢ aᵢ
σ̂_r  = MAD(resid)                       # noise floor after explaining
z_i^r = ⟨resid, a_i⟩/(σ̂_r ‖a_i‖) for all i ∉ support
admit { i : z_i^r > z_min } up to K; refit NNLS
stop when no admissions or max rounds
```

Per-taxon final stats: coverage μ̂ᵢ = x̂ᵢ (occurrence mode ⇒ linear in
coverage), DNA fraction fᵢ = μ̂ᵢ gᵢ / Σⱼ μ̂ⱼ gⱼ, unexplained mass
1 − Σⱼ μ̂ⱼ Lⱼ / n_events (the "novel content" headline number), BH q-values
over tested taxa, flags: `spike_qc`, `genus_complex`, `low_identity_support`
(support present but per-bucket fit poor ⇒ likely a close relative of the
row, not the row itself — v1 heuristic, documented as such).

### 4.6 Genus collapse

If ≥2 selected rows share a genus and pairwise row cosine > 0.5 (computed
at runtime over the tiny candidate sub-matrix), report the genus as one
aggregate row ("S. aureus complex") with summed mass, flag `genus_complex`.
Strain splitting inside a complex is **out of scope v1** and the UI says so.

### 4.7 Fingerprint + DP

Serialize accumulator as i16×scale (quantization RMSE ≤ 0.4% of σ̂ — tested
T3). Optional `--dp-eps`, `--dp-delta`: add N(0, s²·2 ln(1.25/δ)/ε²) with
s = √(ρ·readlen) before quantization; record noise params in header.

---

## 5. Milestones

| M | scope | exit criterion |
|---|---|---|
| **M0** (this session) | hvcore+cli+tests, occurrence mode, sim sweeps, real-isolate smoke test, wasm+worker+UI panel | all T-tests green; LOD curve produced; S. aureus isolate profiles as S. aureus |
| **M1** | distinct-mode dual accumulator + Bloom; count-aware variance; error-rate estimator | measured LOD improvement on sweep; calibration T6 |
| **M2** | real DB buildout: WHO priority + GTDB reps (~3–5k rows); genus metadata; download/OPFS pipeline ≤ 40 MB | 3 Gbp simulated gut sample < 5 min encode in-browser |
| **M3** | AMR second accumulator (full density, D₂=2048) over CARD/VFDB gene rows — reuse of §4.3–4.5 machinery | gene-family detection on the app's synthetic ResFinder benchmark |
| **M4** | baselines + CAMI II (marine, strain-madness) vs Kraken2+Bracken, MetaPhlAn4, sylph, sourmash gather; phase-transition & EV-calibration figures | species F1 within striking distance of sylph at equal DB; honest table |
| **M5** | federation demo (sum-of-fingerprints prevalence), DP ε-LOD curve, paper figures | — |

Timeline estimate: M0 done now; M1 ≈ days; M2 ≈ a week (mostly DB plumbing);
M3–M4 ≈ 2–4 weeks (simulation + writing dominate); M5 follows the paper.

---

## 6. What gets built in M0, file by file

1. `crates/hvcore/src/{hash,encode,db,stats,solve,fp,sim}.rs` + unit tests.
2. `crates/hvcli/src/main.rs` — subcommands:
   - `hypv build-db --fasta-dir DIR --out DB --meta META.tsv`
   - `hypv encode READS.fq[.gz]... --out X.hvf [--dp-eps ε]`
   - `hypv profile X.hvf --db DB --out report.tsv [--json]`
   - `hypv simulate --out DIR --config TOML` (community, abundances,
     coverage, errors, ANI structure)
   - `hypv sweep --db DB --out DIR` (LOD × coverage × ANI grid)
3. `crates/hvwasm/src/lib.rs` — exports `hv_alloc/hv_free/hv_encoder_reset/
   hv_encoder_feed/hv_acc_ptr/hv_acc_f32/hv_profile/hv_fp_bytes`, results as
   packed f64 array; build.rs-free, `--target wasm32-unknown-unknown
   -O3 --export-dynamic`… (optimization flags in §7).
4. `js/hv-worker.js`, `js/hv-runtime.js`, index.html panel, `hvprof/`
   artifacts, README section.

## 7. Build & perf budget

- native: `cargo build --release` (aarch64) — encode ≥ 60 MB/s/min expected;
  profile < 50 ms for N=40, < 1 s for N=4k.
- wasm: `cargo build --release --target wasm32-unknown-unknown` (opt-level=3,
  lto, `mut-globals` default off) → expect ≥ 20–40 MB/s encode in-worker
  ⇒ 51 MB S. aureus subsample ≈ 1–3 s; 3 Gbp ≈ 2–4 min.
- memory: wasm heap = acc (32 KB) + DB matrix (N·D bytes) + workspace ≈
  50 MB at N=4k; chunked feeding means the FASTQ never resides in memory.

## 8. Validation protocol (tests gate milestones)

| id | test | passes when |
|---|---|---|
| T1 | canonical hash symmetry + avalanche (NIST-ish bucket uniformity) | hash(s)==hash(rc(s)) ∀random s; χ² p>1e-3 on 10⁷ draws |
| T2 | selection/bucket/sign independence | mutual χ² independence p>1e-3; density within ±2% of ρ |
| T3 | quantization (i16 fp, int8 DB) vs float64 | profile outputs differ < 1% on fixed mixtures |
| T4 | exact recovery, incoherent case | NNLS recovers mixture, abundances within 25% rel. error (floor-limited members excluded) |
| T5 | linearity | ‖v(A⊕B) − v(A) − v(B)‖ / ‖v(A⊕B)‖ small and shrinking with coverage |
| T6 | null calibration | pure-random sample (or leave-one-out): z nulls ≈ N(0,1) (KS p>1e-2), BH q≥0.05 for all absent taxa (FDR ≤ 5% over 100 trials) |
| T7 | LOD sweep | detection vs abundance × coverage grid → LOD curve; measured v0: ≈0.5× at D=8192, ≈0.2× at D=32768 (10× background) |
| T8 | ANI behaviour | strain at ANI 99.9% → genus_complex flag; relative at 95% → low_identity_support, no species call |
| T9 | real isolates | SRR21386014 → S. aureus dominant (>80% mass); Klebsiella run → Klebsiella; ERR15682208 (mystery) → sensible composition, large unexplained mass reported |
| T10 | native↔wasm concordance | identical fingerprints (bit-exact) for the same FASTQ through CLI and browser worker (repo's existing concordance-harness pattern) |
| F1–F4 | figures (M1+) | phase transition; EV calibration; depletion waterfall; DP ε-vs-LOD |

**Calibration post-mortem (found in review, 2026-08)**: the first
implementation computed final detection z as (raw dot)/(σ̂ *depleted*). Under
the row-absent null a raw dot's variance scales with the **pre-depletion**
σ̂₀ — dividing by the depleted σ̂ inflates every absent row's z by
σ̂₀/σ̂_resid ≈ 2–3× on isolate-like samples, precisely where FDR matters.
The T6 global-null test could not catch this (nothing fits ⇒ σ̂ never
shrinks). Correction: all BH p-values and detection z use σ̂₀; the depleted
σ̂ remains a QC number. Measured floor moved 0.1×→0.2× at D=32768 — the
honest number. **Lesson: null-calibration tests must include the
mixed presence/absence case, not only the global null** (now `hypv bench
fdr`: 0 FP / ~26 000 nulls at N=2000, flat floor from N=26→2000).

## 9. Risks & honest limitations (keep this list in the paper)

1. **Strain resolution**: near-identical rows are collinear by construction;
   v1 reports genus-level complexes with `genus_complex`. (Azeotrope
   metaphor from §1.2: some mixtures are inseparable by this "distillation";
   a second pass at higher D over the complex only is the v2 cure.)
2. **Novel relatives**: a 95%-ANI relative deposits ~ρ·a^k ≈ 20% of its
   k-31 mass into the nearest row. Detection stays calibrated (T6) but the
   *identity* claim weakens — hence `low_identity_support` and the paper's
   framing: HV is a **detection + mass** instrument, not a strain typer.
3. **Occurrence-mode noise at high coverage** (§1.6) — v0 limitation,
   measured and then fixed by dual accumulators in M1.
4. **PCR duplicates and overlapping mate pairs** double-count k-mer
   occurrences, inflating fitted coverage and multiplicity noise — exactly
   the Σcount² term of §1.6. The published sylph (Nat Biotechnol 2025)
   credits its Illumina accuracy to dedicated duplicate handling: an LSH
   tuple scheme (FracMinHash k-mer + alternating-mask 16-mers at read
   starts) in a scalable cuckoo filter, plus counting a k-mer seen in both
   mates of an overlapping pair only once. Port for M1 (read-level dedup
   before accumulation; cheap in the encoder).
5. **Reference-row k-mer correlation**: adjacent selected k-mers share
   reads, violating the independence behind our variance model. sylph drops
   reference k-mers <30 bp apart; at ρ=0.05 our selections average ~20 bp
   apart. M1 DB builder should enforce a minimum spacing.
6. **DB breadth vs coherence**: adding thousands of rows raises N_eff and
   the EV screen bound; coherence within genera grows. Mitigation: genus
   collapse + measured FDR at DB size (T6 run at N=4k in M4). Measured so
   far: floor flat N=26→2000, 0 FP/~26k nulls.
7. **AMR layer density budget** (M3): gene rows at ρ=0.05 carry only ~10²
   events/gene — below the floor for low-abundance organisms; hence the
   separate full-density accumulator. We will not claim single-sketch AMR.
8. **Federated DP**: composition privacy is not provided (by design);
   membership-inference against the fingerprint is possible and is literally
   the product's purpose — the paper must say this in one sentence.
9. **Browser variance**: iOS Safari memory ceilings; mitigate by chunk size
   and keeping matrix in wasm heap only.

## 9.1 Lessons imported from the published sylph (Nat Biotechnol 2025)

- Terminology: our fitted "coverage" x̂ is exactly their **effective
  coverage λ** (δ·E[(1−ε)^k]·(L−k+1)/L) — sequence depth after error and
  read-length attenuation. Align output naming and docs to "effective
  coverage" for comparability.
- Their λ estimator is a **ratio-of-multiplicity histogram**
  (λ̂ = (N_{a+1}/N_a)·(a+1)), which is identity-independent — the trick that
  makes zero-inflation separable from coverage. This is the statistic our
  occurrence-mode vector cannot form (we discard multiplicity structure in
  the projection); the M1 dual accumulator's distinct/occurrence ratio is
  our route to the same separation.
- Their **E-estimator** (E = 1 − n₁/Σ_{a>1} a·n_a, error attenuation from
  the singleton fraction of the sample sketch) ports directly to our encoder
  stats and feeds the count-aware variance model (M1).
- Their real-data detection crossover vs MetaPhlAn4 sits at **~0.03×
  effective coverage** (Fig. 4d) — the honest reference point for our floor
  claims (v0 calibrated: ~0.2× at D=32768; sylph retains per-k-mer
  information we compress away).
- Their **ANI-as-covariate MWAS** validates our analogous claim that a
  non-compositional statistic (our z) can support association testing on
  aggregated fingerprints — the federation pitch has a published precedent
  in spirit.
- Their profiling step is **winner-take-all k-mer reassignment**, not joint
  deconvolution — the NNLS joint solve remains a genuine HV
  differentiator vs the published sylph.
