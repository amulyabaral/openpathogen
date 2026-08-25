//! hypv — hypervector profiler CLI.
//!
//! Subcommands:
//!   build-db   FASTA dir + manifest → .hvd reference matrix + sidecar TSV
//!   encode     FASTQ(.gz)… → .hvf fingerprint (± DP noise)
//!   profile    .hvf + .hvd → TSV report on stdout
//!   simulate   synthetic community + reads (ground truth by construction)
//!   sweep      LOD grid over coverage/abundance/ANI
//!   agg        sum fingerprints (federation primitive)
//!   spike      emit the synthetic spike-in genome as FASTA

use flate2::read::MultiGzDecoder;
use hvcore::db::{HvDb, RowMeta, FLAG_SPIKE};
use hvcore::encode::{Accumulator, EncodeParams, FastqScanner};
use hvcore::fp::Fingerprint;
use hvcore::sim::*;
use hvcore::solve::{profile, ProfileParams};
use std::fs;
use std::io::{Read, Write};
use std::path::{Path, PathBuf};
use std::time::Instant;

fn main() {
    let args: Vec<String> = std::env::args().collect();
    if args.len() < 2 {
        usage();
    }
    let cmd = args[1].as_str();
    let rest = &args[2..];
    let t0 = Instant::now();
    let result = match cmd {
        "build-db" => cmd_build_db(rest),
        "encode" => cmd_encode(rest),
        "profile" => cmd_profile(rest),
        "simulate" => cmd_simulate(rest),
        "sweep" => cmd_sweep(rest),
        "bench" => cmd_bench(rest),
        "agg" => cmd_agg(rest),
        "spike" => cmd_spike(rest),
        _ => usage(),
    };
    if let Err(e) = result {
        eprintln!("error: {e}");
        std::process::exit(1);
    }
    eprintln!("[hypv] done in {:.1}s", t0.elapsed().as_secs_f64());
}

fn usage() -> ! {
    eprintln!("usage: hypv <build-db|encode|profile|simulate|sweep|bench|agg|spike> …");
    std::process::exit(2);
}

// ── arg helpers ──────────────────────────────────────────────────────────

fn arg_val(args: &[String], name: &str) -> Option<String> {
    args.iter()
        .position(|a| a == name)
        .and_then(|i| args.get(i + 1).cloned())
}

fn arg_f64(args: &[String], name: &str, default: f64) -> f64 {
    arg_val(args, name).and_then(|v| v.parse().ok()).unwrap_or(default)
}

fn arg_u64(args: &[String], name: &str, default: u64) -> u64 {
    arg_val(args, name).and_then(|v| v.parse().ok()).unwrap_or(default)
}

/// Positional args: skips `--flag` tokens AND their following value.
fn positional(args: &[String]) -> Vec<String> {
    let mut out = Vec::new();
    let mut i = 0;
    while i < args.len() {
        if args[i].starts_with("--") {
            i += 2; // flag + value
        } else {
            out.push(args[i].clone());
            i += 1;
        }
    }
    out
}

fn encode_params(args: &[String]) -> EncodeParams {
    let mut p = EncodeParams::default_params();
    if let Some(v) = arg_val(args, "--dim") {
        p.dim_bits = v.parse().unwrap();
    }
    if let Some(v) = arg_val(args, "--density") {
        let d: f64 = v.parse().unwrap();
        p.dens_num = (d * 1000.0).round() as u32;
        p.dens_den = 1000;
    }
    if let Some(v) = arg_val(args, "--qmin") {
        p.qmin = v.parse().unwrap();
    }
    if let Some(v) = arg_val(args, "--seed") {
        p.seed = v.parse().unwrap();
    }
    p
}

fn read_maybe_gz(path: &Path) -> std::io::Result<Vec<u8>> {
    let raw = fs::read(path)?;
    if raw.len() >= 2 && raw[0] == 0x1f && raw[1] == 0x8b {
        let mut out = Vec::new();
        MultiGzDecoder::new(&raw[..]).read_to_end(&mut out)?;
        Ok(out)
    } else {
        Ok(raw)
    }
}

// ── build-db ─────────────────────────────────────────────────────────────
//
// manifest TSV (one row per genome, in row order):
//   fasta_path<TAB>id<TAB>name<TAB>genus[<TAB>flags]

fn cmd_build_db(args: &[String]) -> Result<(), String> {
    let pos = positional(args);
    let manifest = pos
        .first()
        .ok_or("usage: hypv build-db MANIFEST.tsv --out DBPREFIX [--spike]")?;
    let out = arg_val(args, "--out").ok_or("missing --out")?;
    let params = encode_params(args);

    let text = fs::read_to_string(manifest).map_err(|e| e.to_string())?;
    let mut genomes = Vec::new();
    for (ln, line) in text.lines().enumerate() {
        if line.trim().is_empty() || line.starts_with('#') {
            continue;
        }
        let f: Vec<&str> = line.split('\t').collect();
        if f.len() < 4 {
            return Err(format!("manifest line {}: need ≥4 fields", ln + 1));
        }
        let bytes = read_maybe_gz(Path::new(f[0].trim())).map_err(|e| format!("{}: {e}", f[0]))?;
        let meta = RowMeta {
            id: f[1].trim().to_string(),
            name: f[2].trim().to_string(),
            genus: f[3].trim().to_string(),
            len: 0,
            l_sel: 0,
            flags: f.get(4).and_then(|s| s.trim().parse().ok()).unwrap_or(0),
        };
        eprintln!("[build-db] {} ({:.1} kb)", meta.id, bytes.len() as f64 / 1024.0);
        genomes.push((bytes, meta));
    }
    if positional(args).iter().any(|_| arg_val(args, "--spike").is_some())
        || args.iter().any(|a| a == "--spike")
    {
        let meta = RowMeta {
            id: "hypv_spike".into(),
            name: "synthetic spike-in (QC)".into(),
            genus: "HypvSpike".into(),
            len: 0,
            l_sel: 0,
            flags: FLAG_SPIKE,
        };
        genomes.push((hvcore::db::spike_genome(1, 100_000), meta));
    }
    let t = Instant::now();
    let db = HvDb::build(&params, genomes);
    fs::write(format!("{out}.hvd"), db.to_bytes()).map_err(|e| e.to_string())?;
    fs::write(format!("{out}.tsv"), db.meta_tsv()).map_err(|e| e.to_string())?;
    eprintln!(
        "[build-db] {} rows, D={}, int8 {:.2} MB, build {:.1}s",
        db.n(),
        db.dim(),
        db.rows.len() as f64 / 1e6,
        t.elapsed().as_secs_f64()
    );
    Ok(())
}

// ── encode ───────────────────────────────────────────────────────────────

fn cmd_encode(args: &[String]) -> Result<(), String> {
    let pos = positional(args);
    if pos.is_empty() {
        return Err("usage: hypv encode READS.fq.gz… --out X.hvf [--dp-eps e]".into());
    }
    let out = arg_val(args, "--out").ok_or("missing --out")?;
    let params = encode_params(args);
    let mut acc = Accumulator::new(params);
    let t = Instant::now();
    for path in &pos {
        let bytes = read_maybe_gz(Path::new(path)).map_err(|e| e.to_string())?;
        let mut sc = FastqScanner::new();
        sc.feed(&bytes, &mut acc);
        sc.finish(&mut acc);
        eprintln!(
            "[encode] {} → {} reads, {} events ({:.1} MB/s)",
            path,
            acc.stats.n_reads,
            acc.stats.n_events,
            bytes.len() as f64 / 1e6 / t.elapsed().as_secs_f64().max(0.001)
        );
    }
    let dp = arg_val(args, "--dp-eps").map(|eps| {
        let eps: f64 = eps.parse().unwrap();
        (eps, 1e-6, params.seed ^ 0xD00D)
    });
    let fp = Fingerprint::from_accumulator(&acc, dp);
    fs::write(&out, fp.to_bytes()).map_err(|e| e.to_string())?;
    eprintln!(
        "[encode] fingerprint {} ({} bases, {} events, {:.1} kB)",
        out,
        fp.n_bases,
        fp.n_events,
        fp.to_bytes().len() as f64 / 1024.0
    );
    Ok(())
}

// ── profile ──────────────────────────────────────────────────────────────

fn cmd_profile(args: &[String]) -> Result<(), String> {
    let pos = positional(args);
    if pos.len() < 2 {
        return Err("usage: hypv profile X.hvf DBPREFIX [--zmin 5] [--json]".into());
    }
    let fp_bytes = fs::read(&pos[0]).map_err(|e| e.to_string())?;
    let fp = Fingerprint::from_bytes(&fp_bytes)?;
    let db_bytes = fs::read(format!("{}.hvd", pos[1])).map_err(|e| e.to_string())?;
    let mut db = HvDb::from_bytes(&db_bytes)?;
    let meta = fs::read(format!("{}.tsv", pos[1])).unwrap_or_default();
    db.attach_meta_tsv(&meta)?;

    let v = fp.vector();
    let pp = ProfileParams {
        zmin: arg_f64(args, "--zmin", 5.0),
        ..Default::default()
    };
    let report = profile(&v, fp.n_events, &db, pp);

    let json = args.iter().any(|a| a == "--json");
    if json {
        println!("{{");
        println!("  \"n_events\": {},", report.n_events);
        println!("  \"sigma\": {:.4},", report.sigma);
        println!("  \"unexplained\": {:.4},", report.unexplained);
        println!("  \"spike_z\": {:.3},", report.spike_z);
        println!("  \"taxa\": [");
        for (i, t) in report.taxa.iter().enumerate() {
            let m = &db.metas[t.idx];
            println!(
                "    {{\"id\": \"{}\", \"name\": \"{}\", \"coverage\": {:.4}, \"z\": {:.2}, \"q\": {:.2e}, \"dna_frac\": {:.5}, \"flags\": {}}}{}",
                m.id, m.name, t.x, t.z, t.q, t.dna_frac, t.flags,
                if i + 1 < report.taxa.len() { "," } else { "" }
            );
        }
        println!("  ]");
        println!("}}");
    } else {
        println!(
            "# taxa={} events={} sigma={:.3} unexplained={:.4} spike_z={:.2}",
            report.taxa.len(),
            report.n_events,
            report.sigma,
            report.unexplained,
            report.spike_z
        );
        println!("id\tname\tcoverage\tz\tq\tdna_frac\tflags");
        for t in &report.taxa {
            let m = &db.metas[t.idx];
            println!(
                "{}\t{}\t{:.4}\t{:.2}\t{:.3e}\t{:.5}\t{}",
                m.id, m.name, t.x, t.z, t.q, t.dna_frac, t.flags
            );
        }
    }
    Ok(())
}

// ── simulate ─────────────────────────────────────────────────────────────
//
// hypv simulate --out DIR [--n 8] [--glen 300000] [--cov "5,5,…|5"]
//               [--err 0.001] [--seed 1] [--novel-rate 0.0] [--read-len 150]
// Writes genomes/, manifest.tsv (for build-db), truth.tsv, reads_1/2.fq.gz.

fn cmd_simulate(args: &[String]) -> Result<(), String> {
    let out_dir = PathBuf::from(arg_val(args, "--out").ok_or("missing --out")?);
    let n = arg_u64(args, "--n", 8) as usize;
    let glen = arg_u64(args, "--glen", 300_000) as usize;
    let err = arg_f64(args, "--err", 0.001);
    let seed = arg_u64(args, "--seed", 1);
    let cov_spec = arg_val(args, "--cov").unwrap_or_else(|| "5".to_string());
    let novel_rate = arg_f64(args, "--novel-rate", 0.0);
    let novel_cov = arg_f64(args, "--novel-cov", 1.0);
    let read_len = arg_u64(args, "--read-len", 150) as usize;

    fs::create_dir_all(out_dir.join("genomes")).map_err(|e| e.to_string())?;
    let mut rng = Rng::new(seed);

    let covs: Vec<f64> = if cov_spec.contains(',') {
        cov_spec.split(',').map(|s| s.parse().unwrap()).collect()
    } else {
        vec![cov_spec.parse().unwrap(); n]
    };
    if covs.len() != n {
        return Err(format!("--cov has {} entries, need {n}", covs.len()));
    }

    let mut manifest = String::from("#fasta\tid\tname\tgenus\tflags\n");
    let mut truth = String::from("id\tcoverage\n");
    let mut genome_seqs = Vec::new();
    let mut all_covs: Vec<f64> = Vec::new();
    for i in 0..n {
        let g = synth_genome(&mut rng, glen, 0.5);
        let id = format!("sp{i}");
        let fasta = genome_to_fasta(&id, &g);
        let path = out_dir.join("genomes").join(format!("{id}.fasta"));
        fs::write(&path, fasta).map_err(|e| e.to_string())?;
        manifest.push_str(&format!(
            "{}\t{}\tGenus{i} species{i}\tGenus{i}\t0\n",
            path.display(),
            id
        ));
        truth.push_str(&format!("{id}\t{}\n", covs[i]));
        genome_seqs.push(g);
        all_covs.push(covs[i]);
    }
    if novel_rate > 0.0 {
        let novel = mutate(&mut rng, &genome_seqs[0], novel_rate);
        genome_seqs.push(novel);
        truth.push_str(&format!("NOVEL(rate={novel_rate})\t{novel_cov}\n"));
        all_covs.push(novel_cov);
    }
    fs::write(out_dir.join("manifest.tsv"), manifest).map_err(|e| e.to_string())?;
    fs::write(out_dir.join("truth.tsv"), truth).map_err(|e| e.to_string())?;

    let rs = ReadSimParams {
        read_len,
        err_rate: err,
        ..Default::default()
    };
    use flate2::write::GzEncoder;
    use flate2::Compression;
    let write_gz = |data: &[u8], path: &Path| -> Result<(), String> {
        let mut e = GzEncoder::new(Vec::new(), Compression::fast());
        e.write_all(data).map_err(|e| e.to_string())?;
        fs::write(path, e.finish().map_err(|e| e.to_string())?).map_err(|e| e.to_string())
    };
    let mut r1 = Vec::new();
    let mut r2 = Vec::new();
    for (gi, g) in genome_seqs.iter().enumerate() {
        let cov = all_covs[gi];
        if cov <= 0.0 {
            continue;
        }
        let reads = sim_reads(&mut rng, g, cov, &rs, gi);
        for r in &reads {
            let tgt = if r.name.ends_with('1') { &mut r1 } else { &mut r2 };
            tgt.extend_from_slice(&reads_to_fastq(&[r.clone()]));
        }
    }
    write_gz(&r1, &out_dir.join("reads_1.fq.gz"))?;
    write_gz(&r2, &out_dir.join("reads_2.fq.gz"))?;
    eprintln!(
        "[simulate] {} genomes × {} bp → {} (reads_1/2.fq.gz + manifest + truth)",
        genome_seqs.len(),
        glen,
        out_dir.display()
    );
    Ok(())
}

// ── sweep ────────────────────────────────────────────────────────────────
//
// LOD sweep: a target genome at descending coverage over a fixed background
// community, all rows in-DB. Measures the detection floor (T7 en masse).

fn cmd_sweep(args: &[String]) -> Result<(), String> {
    let pos = positional(args);
    if pos.is_empty() {
        return Err("usage: hypv sweep OUT_DIR [--n 6] [--glen 300000]".into());
    }
    let out_dir = PathBuf::from(&pos[0]);
    let n = arg_u64(args, "--n", 6) as usize;
    let glen = arg_u64(args, "--glen", 300_000) as usize;
    let bg = arg_f64(args, "--bg", 2.0);
    fs::create_dir_all(&out_dir).map_err(|e| e.to_string())?;

    let params = encode_params(args);
    let pp = ProfileParams::default_params();
    let mut rng = Rng::new(arg_u64(args, "--seed", 7));

    // genomes first, then the DB from those exact genomes
    let genome_seqs: Vec<Vec<u8>> = (0..n)
        .map(|_| synth_genome(&mut rng, glen, 0.5))
        .collect();
    let genomes: Vec<(Vec<u8>, RowMeta)> = genome_seqs
        .iter()
        .enumerate()
        .map(|(i, g)| {
            (
                genome_to_fasta(&format!("sp{i}"), g),
                RowMeta {
                    id: format!("sp{i}"),
                    name: format!("Genus{i} species{i}"),
                    genus: format!("Genus{i}"),
                    len: g.len() as u32,
                    l_sel: 0,
                    flags: 0,
                },
            )
        })
        .collect();
    let db = HvDb::build(&params, genomes);

    let rs = ReadSimParams::default();
    let mut tsv = String::from("target_cov\tbg_cov\tdetected\tz\tq\tcoverage_est\tn_taxa\n");
    for &target_cov in &[3.0f64, 1.0, 0.5, 0.2, 0.1, 0.05, 0.02, 0.01, 0.005] {
        let mut acc = Accumulator::new(params);
        for (i, g) in genome_seqs.iter().enumerate().skip(1) {
            let reads = sim_reads(&mut rng, g, bg, &rs, i);
            let fq = reads_to_fastq(&reads);
            let mut sc = FastqScanner::new();
            sc.feed(&fq, &mut acc);
        }
        let reads = sim_reads(&mut rng, &genome_seqs[0], target_cov, &rs, 0);
        let fq = reads_to_fastq(&reads);
        let mut sc = FastqScanner::new();
        sc.feed(&fq, &mut acc);
        let report = profile(&acc.acc, acc.stats.n_events, &db, pp);
        let hit = report.taxa.iter().find(|t| t.idx == 0);
        tsv.push_str(&format!(
            "{}\t{bg}\t{}\t{:.2}\t{:.3e}\t{:.4}\t{}\n",
            target_cov,
            hit.is_some(),
            hit.map(|h| h.z).unwrap_or(0.0),
            hit.map(|h| h.q).unwrap_or(1.0),
            hit.map(|h| h.x).unwrap_or(0.0),
            report.taxa.len()
        ));
        println!("{target_cov}\t{}", hit.is_some());
    }
    let out = out_dir.join("sweep.tsv");
    fs::write(&out, tsv).map_err(|e| e.to_string())?;
    eprintln!("[sweep] wrote {}", out.display());
    Ok(())
}

// ── bench ────────────────────────────────────────────────────────────────
//
// Benchmark battery behind PLAN.md §8. All modes build synthetic ground
// truth in-memory (no FASTQ files) at the requested DB size.
//
//   hypv bench fdr  --n 500 --present 8 --trials 15   → false-discovery rate
//   hypv bench lodn --ns 26,500,2000                  → detection floor vs N
//   hypv bench dp   --cov 5 --eps 0.1,0.3,1,3,10      → DP utility curve
//   hypv bench agg                                      → fingerprint aggregation
//   hypv bench runtime --n 2000                        → throughput scaling

fn cmd_bench(args: &[String]) -> Result<(), String> {
    let pos = positional(args);
    let mode = pos.first().cloned().unwrap_or_default();
    match mode.as_str() {
        "fdr" => bench_fdr(args),
        "lodn" => bench_lodn(args),
        "dp" => bench_dp(args),
        "agg" => bench_agg(args),
        "runtime" => bench_runtime(args),
        "distinct" => bench_distinct(args),
        "earlyexit" => bench_earlyexit(args),
        _ => Err("bench modes: fdr | lodn | dp | agg | runtime | distinct | earlyexit".into()),
    }
}

struct BenchWorld {
    db: HvDb,
    genomes: Vec<Vec<u8>>,
    novel: Vec<u8>, // background organism absent from the DB
    params: EncodeParams,
}

fn bench_world(n: usize, glen: usize, seed: u64) -> BenchWorld {
    let mut params = EncodeParams::default_params();
    if let Ok(v) = std::env::var("HV_DIM_BITS") {
        params.dim_bits = v.parse().unwrap_or(13);
    }
    if let Ok(v) = std::env::var("HV_DENS_NUM") {
        params.dens_num = v.parse().unwrap_or(5);
    }
    let mut rng = Rng::new(seed);
    let mut genomes = Vec::with_capacity(n);
    let mut fastas = Vec::with_capacity(n);
    for i in 0..n {
        let g = synth_genome(&mut rng, glen, 0.5);
        genomes.push(g.clone());
        fastas.push((
            genome_to_fasta(&format!("sp{i}"), &g),
            RowMeta {
                id: format!("sp{i}"),
                name: format!("Genus{i} species{i}"),
                genus: format!("Genus{i}"),
                len: glen as u32,
                l_sel: 0,
                flags: 0,
            },
        ));
    }
    let db = HvDb::build(&params, fastas);
    let novel = synth_genome(&mut rng, glen, 0.5);
    BenchWorld { db, genomes, novel, params }
}

fn encode_comm(
    w: &BenchWorld,
    members: &[(usize, f64)],   // (row idx, coverage)
    novel_cov: f64,
    err: f64,
    rng: &mut Rng,
) -> Accumulator {
    let mut acc = Accumulator::new(w.params);
    let rs = ReadSimParams { err_rate: err, ..Default::default() };
    let mut sc = FastqScanner::new();
    for &(idx, cov) in members {
        let reads = sim_reads(rng, &w.genomes[idx], cov, &rs, idx);
        let fq = reads_to_fastq(&reads);
        sc.feed(&fq, &mut acc);
    }
    if novel_cov > 0.0 {
        let reads = sim_reads(rng, &w.novel, novel_cov, &rs, 999);
        let fq = reads_to_fastq(&reads);
        sc.feed(&fq, &mut acc);
    }
    sc.finish(&mut acc);
    acc
}

/// FDR: present members + absent rows + novel background, many trials.
/// The metric that would expose any anti-conservative z inflation.
fn bench_fdr(args: &[String]) -> Result<(), String> {
    let n = arg_u64(args, "--n", 500) as usize;
    let present = arg_u64(args, "--present", 8) as usize;
    let trials = arg_u64(args, "--trials", 15) as usize;
    let cov = arg_f64(args, "--cov", 4.0);
    let novel_cov = arg_f64(args, "--novel-cov", 2.0);
    let glen = arg_u64(args, "--glen", 200_000) as usize;
    let seed = arg_u64(args, "--seed", 1);
    eprintln!("[bench fdr] building N={n} DB…");
    let w = bench_world(n, glen, seed);
    let pp = ProfileParams::default_params();

    let mut rng = Rng::new(seed ^ 0xF00D);
    let mut tp = 0u32;
    let mut fp = 0u32;
    let mut fnn = 0u32;
    println!("trial\tpresent_reported\tabsent_reported\tmissed\tscreen_max_z");
    for t in 0..trials {
        // random distinct present rows
        let mut members: Vec<(usize, f64)> = Vec::new();
        let mut seen = std::collections::HashSet::new();
        while members.len() < present {
            let i = rng.below(n as u64) as usize;
            if seen.insert(i) {
                members.push((i, cov * (0.5 + rng.unit() * 1.5)));
            }
        }
        let acc = encode_comm(&w, &members, novel_cov, 0.001, &mut rng);
        let rep = profile(&acc.acc, acc.stats.n_events, &w.db, pp);
        let present_set: std::collections::HashSet<usize> =
            members.iter().map(|m| m.0).collect();
        let mut t_tpr = 0u32;
        let mut t_fp = 0u32;
        for tax in &rep.taxa {
            if present_set.contains(&tax.idx) {
                t_tpr += 1;
            } else {
                t_fp += 1;
                eprintln!(
                    "  trial {t}: FALSE POSITIVE row {} z={:.2} q={:.2e}",
                    tax.idx, tax.z, tax.q
                );
            }
        }
        tp += t_tpr;
        fp += t_fp;
        fnn += (present - t_tpr as usize) as u32;
        // screen max (includes present rows — a scale diagnostic, not a null)
        let max_null = rep.screen_max;
        println!("{t}\t{t_tpr}\t{t_fp}\t{}\t{max_null:.2}", present - t_tpr as usize);
    }
    let fdr = fp as f64 / (fp + tp).max(1) as f64;
    println!("# N={n} present={present} trials={trials} cov≈{cov} novel={novel_cov}");
    println!(
        "# TP={tp} FP={fp} FN={fnn} FDR={fdr:.4} (target <= 0.05) sensitivity={:.3}",
        tp as f64 / (tp + fnn) as f64
    );
    if fdr > 0.05 {
        eprintln!("[bench fdr] FDR ABOVE TARGET — calibration problem");
    }
    Ok(())
}

/// Detection floor as a function of DB size N (EV bound grows with N).
fn bench_lodn(args: &[String]) -> Result<(), String> {
    let ns_spec = arg_val(args, "--ns").unwrap_or_else(|| "26,500,2000".into());
    let glen = arg_u64(args, "--glen", 200_000) as usize;
    let bg = arg_f64(args, "--bg", 2.0);
    let seed = arg_u64(args, "--seed", 5);
    let pp = ProfileParams::default_params();
    println!("n\tfloor_cov\tbg_total");
    for n_spec in ns_spec.split(',') {
        let n: usize = n_spec.trim().parse().map_err(|_| "bad --ns")?;
        eprintln!("[bench lodn] N={n}…");
        let w = bench_world(n, glen, seed + n as u64);
        let mut rng = Rng::new(seed ^ n as u64);
        // background: 5 fixed DB genomes + target ladder on row 0
        let bg_rows: Vec<usize> = (1..=5.min(n - 1)).collect();
        let mut floor_cov = None;
        for &target in &[2.0f64, 1.0, 0.5, 0.25, 0.1, 0.05] {
            let mut members: Vec<(usize, f64)> =
                bg_rows.iter().map(|&i| (i, bg)).collect();
            members.push((0, target));
            let acc = encode_comm(&w, &members, 0.0, 0.001, &mut rng);
            let rep = profile(&acc.acc, acc.stats.n_events, &w.db, pp);
            if rep.taxa.iter().any(|t| t.idx == 0 && t.q < 0.05) {
                floor_cov = Some(target);
                // continue down to find the true floor
            }
        }
        // floor = smallest target still detected; re-scan ascending
        let mut found = None;
        for &target in &[0.05f64, 0.1, 0.25, 0.5, 1.0, 2.0] {
            let mut members: Vec<(usize, f64)> = bg_rows.iter().map(|&i| (i, bg)).collect();
            members.push((0, target));
            let acc = encode_comm(&w, &members, 0.0, 0.001, &mut rng);
            let rep = profile(&acc.acc, acc.stats.n_events, &w.db, pp);
            if rep.taxa.iter().any(|t| t.idx == 0 && t.q < 0.05) {
                found = Some(target);
                break;
            }
        }
        let _ = floor_cov;
        println!("{n}\t{}\t{}", found.unwrap_or(f64::NAN), bg * bg_rows.len() as f64);
    }
    Ok(())
}

/// DP utility: detection z of a true member vs ε in the released fingerprint.
fn bench_dp(args: &[String]) -> Result<(), String> {
    let cov = arg_f64(args, "--cov", 5.0);
    let eps_spec = arg_val(args, "--eps").unwrap_or_else(|| "0.1,0.3,1,3,10".into());
    let n = arg_u64(args, "--n", 26) as usize;
    let w = bench_world(n, 200_000, 77);
    let pp = ProfileParams::default_params();
    let mut rng = Rng::new(3);
    let members: Vec<(usize, f64)> = vec![(0, cov), (3, cov / 2.0), (7, cov / 5.0)];
    let acc = encode_comm(&w, &members, 1.0, 0.001, &mut rng);

    let z_of_row = |fp: &Fingerprint, row: usize| -> f64 {
        let v = fp.vector();
        let rep = profile(&v, fp.n_events, &w.db, pp);
        rep.taxa.iter().find(|t| t.idx == row).map(|t| t.z).unwrap_or(0.0)
    };
    let clean = Fingerprint::from_accumulator(&acc, None);
    println!("eps\tz_row0\tz_row3\tz_row7");
    println!("none\t{:.1}\t{:.1}\t{:.1}", z_of_row(&clean, 0), z_of_row(&clean, 3), z_of_row(&clean, 7));
    for eps_s in eps_spec.split(',') {
        let eps: f64 = eps_s.trim().parse().map_err(|_| "bad --eps")?;
        let fp = Fingerprint::from_accumulator(&acc, Some((eps, 1e-6, 42)));
        println!(
            "{eps_s}\t{:.1}\t{:.1}\t{:.1}",
            z_of_row(&fp, 0),
            z_of_row(&fp, 3),
            z_of_row(&fp, 7)
        );
    }
    Ok(())
}

/// Aggregation: sum of member fingerprints deconvolves to the pooled
/// community — the federation primitive (PLAN §1.8).
fn bench_agg(args: &[String]) -> Result<(), String> {
    let n = arg_u64(args, "--n", 12) as usize;
    let w = bench_world(n, 200_000, 99);
    let pp = ProfileParams::default_params();
    let mut rng = Rng::new(8);
    // site A: rows 0,1,2 heavy; site B: rows 3,4 heavy + shared row 0 light
    let site_a = vec![(0usize, 6.0f64), (1, 4.0), (2, 2.0)];
    let site_b = vec![(3usize, 5.0f64), (4, 3.0), (0, 1.0)];
    let mut fps = Vec::new();
    for members in [site_a, site_b] {
        let acc = encode_comm(&w, &members, 0.5, 0.001, &mut rng);
        fps.push(Fingerprint::from_accumulator(&acc, None));
    }
    let refs: Vec<&Fingerprint> = fps.iter().collect();
    let agg = Fingerprint::aggregate(&refs)?;
    let v = agg.vector();
    let rep = profile(&v, agg.n_events, &w.db, pp);
    println!("pooled community (sum of 2 fingerprints):");
    for t in &rep.taxa {
        println!(
            "  {} {}\tcoverage {:.2}\tdna_frac {:.3}\tz {:.1}",
            w.db.metas[t.idx].id,
            if t.flags & 256 != 0 { "(complex)" } else { "" },
            t.x,
            t.dna_frac,
            t.z
        );
    }
    // truth: pooled cov = sum of per-site coverages
    println!("truth: sp0=7.0 sp1=4.0 sp2=2.0 sp3=5.0 sp4=3.0");
    Ok(())
}

/// Distinct-mode prototype (M1 validation): encode the sample counting each
/// DISTINCT selected k-mer once (HashSet — the production plan replaces this
/// with a Bloom filter). Background noise then scales with Σcount instead of
/// Σcount², so at high background coverage the detection floor should drop
/// by ≈ √μ_bg. Reference rows are already distinct-based, so the same NNLS
/// machinery applies; coefficients read as cⱼ = 1−e^(−μⱼ).
fn bench_distinct(args: &[String]) -> Result<(), String> {
    use std::collections::HashSet;
    let n = arg_u64(args, "--n", 6) as usize;
    let glen = arg_u64(args, "--glen", 300_000) as usize;
    let bg = arg_f64(args, "--bg", 30.0);
    let seed = arg_u64(args, "--seed", 13);
    let w = bench_world(n, glen, seed);
    // kmax must exceed the true community size, else the top-K admission
    // silently evicts true members (found the hard way in rich backgrounds)
    let pp = ProfileParams {
        kmax: n + 8,
        ..Default::default()
    };
    let mut rng = Rng::new(seed ^ 0x0D1571C7);
    let rs = ReadSimParams { err_rate: 0.001, ..Default::default() };

    let mut run = |target: f64, distinct: bool| -> (bool, f64, f64) {
        // build the read set once per (target, mode)
        let mut acc = Accumulator::new(w.params);
        let mut seen: HashSet<u64> = HashSet::new();
        let mut scanner = FastqScanner::new();
        for (i, g) in w.genomes.iter().enumerate() {
            let cov = if i == 0 { target } else { bg };
            let reads = sim_reads(&mut rng, g, cov, &rs, i);
            let fq = reads_to_fastq(&reads);
            if distinct {
                // distinct-mode: re-parse with dedup accumulation
                distinct_encode(&fq, &w.params, &mut acc, &mut seen);
            } else {
                scanner.feed(&fq, &mut acc);
            }
        }
        let rep = profile(&acc.acc, acc.stats.n_events, &w.db, pp);
        if std::env::var("HV_DEBUG").is_ok() {
            eprintln!(
                "  [dbg] target={target} distinct={distinct} sigma0={:.2} screen_max={:.2} ev={:.2} n_taxa={} n_events={}",
                rep.sigma,
                rep.screen_max,
                rep.ev_bound,
                rep.taxa.len(),
                rep.n_events
            );
        }
        let hit = rep.taxa.iter().find(|t| t.idx == 0);
        (
            hit.map(|h| h.q < 0.05).unwrap_or(false),
            hit.map(|h| h.z).unwrap_or(0.0),
            hit.map(|h| h.x).unwrap_or(0.0),
        )
    };

    println!("target\tmode\tdetected\tz\tcoef");
    for &target in &[2.0f64, 1.0, 0.5, 0.25, 0.1, 0.05, 0.02, 0.01] {
        let (d_o, z_o, x_o) = run(target, false);
        let (d_d, z_d, x_d) = run(target, true);
        println!("{target}\toccurrence\t{d_o}\t{z_o:.2}\t{x_o:.4}");
        println!("{target}\tdistinct\t{d_d}\t{z_d:.2}\t{x_d:.4}");
    }
    Ok(())
}

fn distinct_encode(fq: &[u8], params: &EncodeParams, acc: &mut Accumulator, seen: &mut std::collections::HashSet<u64>) {
    // minimal FASTQ walk: process lines 2 and 4 of each record
    let mut lines: Vec<&[u8]> = fq.split(|&c| c == b'\n').collect();
    if lines.last() == Some(&&[][..]) {
        lines.pop();
    }
    let mut rh = acc.contig_hasher();
    let proj = hvcore::hash::Projector::new(params.dim_bits, params.dens_num, params.dens_den, params.seed);
    let mut idx = 0usize;
    while idx + 3 < lines.len() + 1 {
        // records of 4: header, seq, +, qual
        if idx + 1 < lines.len() {
            let seq = lines[idx + 1];
            let qual = if idx + 3 < lines.len() { lines[idx + 3] } else { &[] };
            rh.reset();
            for (i, &b) in seq.iter().enumerate() {
                let code = if i < qual.len() && qual[i].saturating_sub(33) < params.qmin {
                    4
                } else {
                    hvcore::hash::code_of(b)
                };
                if let Some(u) = rh.push(code) {
                    let (sel, bucket, sign) = proj.derive(u);
                    if sel && seen.insert(u) {
                        let v = acc.acc[bucket].saturating_add(sign);
                        acc.acc[bucket] = v;
                        acc.stats.n_events += 1;
                    }
                }
            }
        }
        idx += 4;
    }
}

/// Early-exit encoding (sequential analysis): how much of the file must be
/// read before the composition call is stable? Demonstrates the Wald-SPRT
/// style optimization on the real S. aureus isolate.
fn bench_earlyexit(args: &[String]) -> Result<(), String> {
    let fq_path = positional(args).get(1).cloned();
    let fq_path = match fq_path {
        Some(p) => p,
        None => return Err("usage: hypv bench earlyexit READS.fq.gz".into()),
    };
    let raw = read_maybe_gz(Path::new(&fq_path)).map_err(|e| e.to_string())?;
    let db_bytes = fs::read("../hvprof/hv_patho_v1.hvd").map_err(|e| e.to_string())?;
    let mut db = HvDb::from_bytes(&db_bytes)?;
    let meta = fs::read("../hvprof/hv_patho_v1.tsv").unwrap_or_default();
    db.attach_meta_tsv(&meta)?;

    println!("fraction\tz_saureus\tdetected\tunexplained");
    for &frac in &[0.02f64, 0.05, 0.10, 0.25, 0.50, 1.0] {
        let cut = ((raw.len() as f64 * frac) as usize / 4) * 4; // record-aligned
        let mut acc = Accumulator::new(db.params);
        let mut sc = FastqScanner::new();
        sc.feed(&raw[..cut.min(raw.len())], &mut acc);
        sc.finish(&mut acc);
        let rep = profile(&acc.acc, acc.stats.n_events, &db, ProfileParams::default_params());
        let hit = rep.taxa.iter().find(|t| db.metas[t.idx].id == "staaur");
        println!(
            "{frac:.2}\t{:.1}\t{}\t{:.3}",
            hit.map(|h| h.z).unwrap_or(0.0),
            hit.map(|h| h.q < 0.05).unwrap_or(false),
            rep.unexplained
        );
    }
    Ok(())
}

/// Runtime scaling: DB build, encode throughput, profile time vs N.
fn bench_runtime(args: &[String]) -> Result<(), String> {
    let n = arg_u64(args, "--n", 2000) as usize;
    let glen = arg_u64(args, "--glen", 200_000) as usize;
    let seed = arg_u64(args, "--seed", 1);
    let t0 = Instant::now();
    let w = bench_world(n, glen, seed);
    eprintln!("[bench runtime] world built ({} genomes × {} bp) in {:.1}s", n, glen, t0.elapsed().as_secs_f64());

    let mut rng = Rng::new(seed);
    let members: Vec<(usize, f64)> = (0..10).map(|i| (i, 3.0)).collect();
    let t1 = Instant::now();
    let acc = encode_comm(&w, &members, 1.0, 0.001, &mut rng);
    let enc_s = t1.elapsed().as_secs_f64();
    let bases = acc.stats.n_bases;
    println!("encode\t{:.2} MB/s\t({bases} bases)", bases as f64 / 1e6 / enc_s);

    let pp = ProfileParams::default_params();
    let t2 = Instant::now();
    let rep = profile(&acc.acc, acc.stats.n_events, &w.db, pp);
    println!("profile_N{}\t{:.2} s\t({} taxa)", n, t2.elapsed().as_secs_f64(), rep.taxa.len());
    Ok(())
}



fn cmd_agg(args: &[String]) -> Result<(), String> {
    let pos = positional(args);
    if pos.len() < 2 {
        return Err("usage: hypv agg OUT.hvf IN1.hvf IN2.hvf…".into());
    }
    let mut fps = Vec::new();
    for p in &pos[1..] {
        let b = fs::read(p).map_err(|e| e.to_string())?;
        fps.push(Fingerprint::from_bytes(&b)?);
    }
    let refs: Vec<&Fingerprint> = fps.iter().collect();
    let agg = Fingerprint::aggregate(&refs)?;
    fs::write(&pos[0], agg.to_bytes()).map_err(|e| e.to_string())?;
    eprintln!("[agg] {} → {}", pos.len() - 1, pos[0]);
    Ok(())
}

// ── spike ────────────────────────────────────────────────────────────────

fn cmd_spike(args: &[String]) -> Result<(), String> {
    let out = positional(args).first().cloned().unwrap_or_else(|| "spike.fasta".into());
    fs::write(&out, hvcore::db::spike_genome(1, 100_000)).map_err(|e| e.to_string())?;
    eprintln!("[spike] wrote {out}");
    Ok(())
}
