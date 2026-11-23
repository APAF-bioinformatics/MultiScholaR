# KIT-159 G1D evidence (ticket section 2.1, with the KIT-159-C13 items 1 and 2 applied to the script)

- Run: 45121359461 (https://github.com/APAF-bioinformatics/MultiScholaR/actions/runs/45121359461), attempt 1, green; script exit 0 after 403 s.
- Commit G: c9c4a412d239738e412dcce2ea7d035a217f1b95 (branch `run-k159-g1d`, parent 12d68ca7; adds only `docs/qa/workflow-execution/kit/kit-159/g1d/g1d.R`).
- Request file: `.github/run-tests/requests/k159-g1d-c9c4a412.json` on `run-tests` (request commit e796de66); tools_ref 76b57c0a, reference BLAS, ubuntu-22.04, AMD EPYC 9V74, 4 cores.
- Sides: base 28d5a4be (oracle's tree), head 12d68ca7 (run 44953187101's base), each in its own worktree and child process under `pkgload::load_all`.
- Files: `status.tsv` (mechanism checks), `g1d-compare.tsv` (DA tables compared, aggregates only), `g1d-captures.tsv` (reproduction checks, per-step capture key counts).

Observed: both sides wrote the same DA table (masked digest 8a4c2dd5... on both; 393 rows, 37 columns, headers, key order, masked lines, every value identical; 19 calls on each side, none changed). Reproduction: head holds (8a4c2dd5 = run 44953187101); base fails (8a4c2dd5, not the oracle's 80cfedaa).
