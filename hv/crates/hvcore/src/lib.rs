//! hvcore — hypervector metagenome profiling core.
//!
//! Pure algorithms over byte slices and typed arrays; no filesystem, no
//! clock, no dependencies. Compiled natively for validation (CLI) and to
//! wasm32 for the browser — the same code produces both sets of numbers.

pub mod db;
pub mod encode;
pub mod fp;
pub mod hash;
pub mod sim;
pub mod solve;
pub mod stats;

pub use db::HvDb;
pub use encode::{Accumulator, EncodeParams, EncodeStats, FastqScanner, FastaScanner};
pub use fp::Fingerprint;
pub use solve::{ProfileParams, ProfileReport, TaxonResult};
