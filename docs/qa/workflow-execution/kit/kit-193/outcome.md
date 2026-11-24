# KIT-193 outcome (in progress: statements gate passed with corrections S01-S07 applied; waiting for its closure check)

Ticket `tickets/kit/KIT-193.md` r1.2. Branch `kit-193-lean-proofs` from `v0.9-staging` at 51959495. Lean runs on the
laptop by the owner's exception (VAL-95, VAL-103); no R was run anywhere for this stage.

## S0: toolchain and Mathlib (acceptance criterion A1)

Observed on 2026-11-24, laptop, every Lean/lake call pinned to four cores (`taskset -c 0-3`, `LEAN_NUM_THREADS=4`;
Lake 5.0.0 has no `-j`/`--jobs` option, read from `lake --help`; `lean --help` lists `-j, --threads`), one call at a
time. CPU time per call: `cpu-time.tsv` (bash's `time`, user and system time over the call and every child it
waited for; GNU `time` is not installed on the laptop).

- Toolchain `leanprover/lean4:v4.36.0` (`proofs/lean/lean-toolchain`), fetched by elan 4.4.0 on the first `lake`
  call: archive 607,649,028 bytes (`lean-4.36.0-linux.tar.zst`, Content-Length), installed 3.1 GB; 60.1 s.
  `lake --version`: `Lake version 5.0.0-src+611c117 (Lean version 4.36.0)`.
- Mathlib `rev = "v4.36.0"` (`proofs/lean/lakefile.toml`) resolved by `lake update` to
  **fc4f59e577feb1a4956671bc04ac18aaaf8cd842**; `git ls-remote https://github.com/leanprover-community/mathlib4
  refs/tags/v4.36.0` returns the same commit (a lightweight tag). Transitive dependencies are pinned by commit in
  `proofs/lean/lake-manifest.json` (plausible, LeanSearchClient, importGraph, proofwidgets, aesop, Qq, batteries,
  Cli).
- Mathlib's cache: `lake update` ran Mathlib's update hook, which built the `cache` executable (61 small files of
  Mathlib's `Cache/`, not Mathlib itself) and downloaded 9,117 files, 460,999,207 bytes (`~/.cache/mathlib`),
  unpacked to 6,817,999,615 bytes (`.lake/packages/mathlib/.lake/build`). Git checkouts of the nine packages:
  172.6 MB. The whole `lake update` took 98.6 s elapsed. `lake exe cache get` afterwards: "No files to download",
  4.9 s.
- `lake build --no-build Mathlib`: "All targets up-to-date (18249 jobs)": every Mathlib module is served from the
  cache; nothing of Mathlib was compiled on the laptop.
- A trivial root module builds with `lake build --wfail` (third attempt; see below).

Two facts of Mathlib v4.36.0 found at S0 (not in the ticket):
1. `Mathlib.Logic.Basic` no longer exists (first trivial build: "no such file").
2. Mathlib uses Lean's module system. A file that is not a `module` and imports Mathlib gets a Lake warning
   ("designed for use with the module system; consider adding `module`"), which `--wfail` turns into a failure.
   Every project file is therefore a `module` with `public import` and `@[expose] public section` (so the Proof
   files can unfold the Spec definitions). The CI record (4.4) inherits this.

## Statements (section 5.3 item 1)

Spec modules, each built alone with `lake build --wfail MultiScholaRProofs.Spec.<name>` and together through the
root (`lake build --wfail`, 6,136 jobs, all from the cache except the project's 9 modules); no Proof file exists.
`grep -rnE 'sorry|admit|axiom |native_decide|implemented_by|extern'` over `MultiScholaRProofs*`: none.

| File | Statements |
|---|---|
| `Spec/FilterCount.lean` | 8 statements, 6 witnesses |
| `Spec/ReplayRule.lean` | 8 statements, 4 witnesses |
| `Spec/GoClosure.lean` | 12 statements, 4 witnesses |
| `Spec/QvalueBh.lean` | 16 statements, 5 witnesses |
| `Spec/Permutation.lean` | 9 statements, 4 witnesses |
| `Spec/RuvAlgebra.lean` | 12 statements (R-A1 to R-E2, R-K1), 5 witnesses |
| `Spec/RuvDownstream.lean` | 11 statements (R-F1 to R-H1, R-M1), 6 witnesses |
| `Spec/RuvIIIC.lean` | 6 statements (R-J1, R-J2, R-L1, R-L2), 4 witnesses |

Statements gate: `tickets/kit/audits/AUDIT-KIT-193-statements-20261124.md` (Fable). Its corrections S01 to S07
are applied as written (82 statements, 38 witnesses after them). One is a result: the ticket's 4d-ii formula is false
at a = 1 (`kept_family_unit_a_counterexample_stmt`; the formula is now stated for 0 < a < 1).

**Statement not yet frozen:** R-I1, R-J3, R-N1. Each gets its own Spec file and statements review before its proof.

`models.md` holds the model-to-code table and each pair's kind and planted defect. `tests/testdata/kit193/bindings.tsv` lists the bound objects; every digest is `pending` until
its pair first runs on merged code (8.2). `SPEC.sha256` is not written: it follows the gate.

## Freeze

`proofs/lean/SPEC.sha256` written from the Spec files at 4183eb10, equal byte for byte to the values in the
statements register's closure check; committed alone (c783e1ae). From here a Spec file changes only through 5.3 item 4.

## Target 3 (GO closure): proved

Observed on the laptop, 2026-11-24 (four cores, one Lean process at a time). Proof file
`MultiScholaRProofs/Proof/GoClosure.lean`; every statement of `Spec/GoClosure.lean` (12) and every witness (4) has a
theorem whose type is the frozen statement (`Audit.lean`'s `example` lines compile), and each depends on
`propext`, `Classical.choice` and `Quot.sound` only (`#print axioms`; `axioms.tsv`).

- **Laptop CPU time for target 3: 39.8 CPU-seconds** (30.5 s user + 9.3 s system) over 11 Lean calls (the oracle
  module's build and run, six compiles of the proof file, the module build, the root build, the audit), **32.5 s of
  Lean wall time**; labels `t3-*` in `cpu-time.tsv`.
- Agent time: 0.12 h for the proofs (11:41 to 11:48), against 14 expected and 24 ceiling; the oracle table
  (`go-closure.tsv`) 0.07 h before it (`effort.tsv`).
- How the search is proved: a loop invariant (every seen pair is a (term, ancestor) pair; every start pair is seen;
  every seen pair off the frontier has all its successors seen; every seen pair is in terms × nodes) and a potential
  (pairs not yet seen, plus one while the frontier is not empty) that falls by at least one on every pass that does not
  stop; so the pass budget never binds and the end state is exactly the closure.

## Not checked

- Targets other than 3 are not proved yet; nothing was run in R.
- The gate's own reads (RUVIIIC at c534151f, `ruv` 0.9.7.1 as locked, RUVSeq) were not done.
- The CI workflow (4.4) and the lifecycle L0 are not written.
