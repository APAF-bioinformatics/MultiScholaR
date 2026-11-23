# KIT-159 gates (section 2)

Each gate is answered in the lifecycle before anything is marked done. Entries are added as gates
are answered; a gate without an entry is not yet answered.

## G1(r): the reference table

Not yet built (it waits for the release commit's work; section 11 order). One input is recorded:

### G1D: the `prot_dia_limpa` diagnosis (section 2.1)

- Run 45121359461 (request `k159-g1d-c9c4a412`, request commit e796de66 on `run-tests`), script
  commit c9c4a412d239738e412dcce2ea7d035a217f1b95 on `run-k159-g1d` (parent 12d68ca7; adds only
  `g1d/g1d.R`), tools_ref 76b57c0a, reference BLAS, ubuntu-22.04. Run by another agent; copied here
  with its returned tables: `g1d/status.tsv`, `g1d/g1d-compare.tsv`, `g1d/g1d-captures.tsv`,
  `g1d/README.md`, and the script `g1d/g1d.R` (byte-identical to c9c4a412's). The tables hold
  counts, column names, largest differences and digests only; no feature row.
- Mechanism: valid (every `status.tsv` row `ok`; script exit 0 after 403 s; no `incompatible`).
- Observed: 28d5a4be (oracle's tree) and 12d68ca7 (run 44953187101's base) wrote the same DA table
  on the runner: masked digest 8a4c2dd5… on both, 393 rows, 37 columns, headers equal, key order
  equal, 0 non-numeric differences, 0 differing masked lines, every numeric difference 0, 0 keys on
  one side only, 19 calls on each side and 0 changed.
- Reproduction checks: head holds (12d68ca7 gives 8a4c2dd5…, run 44953187101's digest); **base
  fails** (28d5a4be gives 8a4c2dd5… on the runner, the stored oracle is 80cfedaa…).
- Decision under 2.1: a reproduction check fails, so the class `limpa_da_last_digits` is **not
  admitted on G1D alone**; decision D19 is the owner's (put to them; no answer recorded here).
- What it shows, within 2.1's limits: on this runner and library, the tree change from 28d5a4be to
  12d68ca7 has no effect on this lane's DA table. It does not attribute the oracle's difference to
  anything (not to decision M29, the engine pin, the CPU, BLAS, the recording machine or other
  inputs), and says nothing about the failing M29 test (`test-dia-limpa-ruv-rules.R`).
- Status of `prot_dia_limpa` for this ticket's work: **awaiting D19**. G1D is not re-sent; nothing
  is done on D19 either way.
