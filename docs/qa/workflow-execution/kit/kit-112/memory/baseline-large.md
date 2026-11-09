# KIT-112 W0: the start, the large baseline and the reachability check (2026-11-10)

**Start condition met.**
- KIT-111 was closed at `59e4a08c`.
- At that commit, `R/`, `tests/`, `DESCRIPTION`, `NAMESPACE`, `inst/` and `src/` are identical to `kit_111`'s source `c568dfa5`. `git diff c568dfa5 59e4a08c` on those paths is empty.

**The large baseline is `111q`, so no `112p` is run.**
- KIT-111 recorded it valid (largest interval 0.21 s, swap 0, slice stopped) and complete.
- Its store is identical to `111p`'s: 0 of 98 comparisons differ (`kit-111/memory/compare-111p-111q.tsv`).
- Its verifier mismatch set is KIT-065's 18.
- The build is `kit_111` (`MlQroA`, `c568dfa5`).

**The screen baselines** are `111l`, `111m` and `111t` ([baseline-screen.tsv](baseline-screen.tsv)). `111m` is the measured member; `111n` (the rerun) is within 0.05 GiB of it on every step except the floor episode that KIT-111 recorded.

## 111q per step

| Step | Controller peak (GiB) | Slice peak (GiB) | Seconds |
|---|---|---|---|
| import | 2.435 | 2.779 | 41.3 |
| design | 2.907 | 3.253 | 107.1 |
| qvalue | 2.913 | 5.078 | 178.2 |
| protein_peptide | 2.172 | 4.212 | 120.9 |
| rollup | 2.173 | 5.813 | 43.7 |
| intensity | 2.173 | 4.101 | 132.0 |
| sample | 2.173 | 3.214 | 35.0 |
| replicate | 2.173 | 4.087 | 116.6 |
| imputation | 2.173 | 4.538 | 468.0 |
| iq | 2.177 | 4.504 | 160.8 |
| accession | 2.182 | 3.039 | 15.2 |
| protein_intensity | 2.179 | 2.989 | 16.0 |
| duplicates | 2.184 | 3.193 | 21.0 |
| protein_replicate | 2.189 | 3.116 | 12.4 |
| normalisation_ruv | 2.626 | 2.987 | 390.4 |
| correlation | 1.806 | 2.449 | 19.0 |
| session_export | 1.708 | 2.070 | 29.9 |
| da | 2.208 | 3.461 | 100.9 |
| report | 2.319 | 2.680 | 27.1 |

## The reachability check (D8)

- **The ticket's rule:** `da` peak 2.208 GiB against `session_export` peak 1.708 GiB. That is +0.500, at least 0.25, so **"better" is not ruled out.**
- **The trace** (`kit-journey-111q/candidate-1/memory-trace.tsv`, 50 ms samples), because KIT-111 showed a step's peak can sit outside the moved work:
  - the controller is flat at 1.70 to 1.71 GiB from `session_export` to the DA click;
  - within about 2 s of the click (the `da` step's start; the driver's `run_da_analysis` action spans the whole 99.6 s step), it rises to 2.02 to 2.03 GiB;
  - about 28 s in, it rises to 2.19 to 2.21, and stays there to the end;
  - from 84 to 96 s, the worker's `da_publication` raises only the slice (to 3.46 GiB), not the controller.
- **Read:** both rises (+0.31 at the click, +0.17 at about 28 s; +0.50 in all) come while the run handler works in the controller. So the "better" target (`da` peak ≤ 1.958 GiB) is reachable only if the controller's request half (D3: the gate, the validation, the commit of the DA settings) does not hydrate the protein matrix or the design. If it does, the +0.31 at the click stays, and the most the move could save is about 0.19 GiB, below the 0.25 threshold.
- **Handed to W2** as a design constraint: the request half reads references and parameters only. Whatever it must hydrate is recorded with its cost.
