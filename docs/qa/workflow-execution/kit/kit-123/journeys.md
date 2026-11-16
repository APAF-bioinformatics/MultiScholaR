# KIT-123 journeys on kit_940 (`Xm75ve`, `652eda20`), 2026-11-16

The KIT-124 simulations were paused in place (SIGSTOP) for the whole run; the profile was changed from balanced to power saver by the user during 940m (heat), and stayed there.

| Member | Route | Judged against | Content | Measurement |
|---|---|---|---|---|
| 940m | DIA_LIMPA, covariates | the kit_940 record (`kit-123/m`) | exact | valid |
| 940r | DIA_LIMPA, correction | the kit_940 record (`kit-123/mc`) | exact | invalid (largest sample gap 1.06 s) |
| 940l | LFQ screen | its oracles; store against store with 939l | exact; 57 of 57 fields identical | invalid (1.27 s) |
| 940t | TMT screen | its oracles; store against store with 939t | exact; 57 of 57 fields identical | invalid (2.86 s) |
| 940q | DIA full | store against store with 939q | 13 of 105 fields differ, all from `ruv_corrected` on | invalid (6.68 s) |

- **Measurement:** the four invalid members fail only `measurement_valid` (sampling gaps over 0.5 s, all after the switch to power saver; 940m, under balanced, has 0.21 s). Every content check passed. KIT-123 claims correctness, not memory, so the content results stand; no memory figure is taken from these members, and no rerun was made (power saver would give the same gaps).
- **LFQ and TMT** are identical to kit_939: the screen journeys' stores are unchanged by R1 and R2 (KIT-123's revision 1.3 note expected their corrected values to move; the store comparison shows they do not on these corpora).
- **DIA full (940q against 939q):** every state before `ruv_corrected` is identical. The sweep's k is 4 (was 45): R1, with `max_acceptable_k` = 5. The table after the correlation filter keeps 2,500 proteins (was 2,476): with the levels restored (R2) the low-intensity filter after the correction no longer removes the 24 proteins the per-protein offsets had pushed under its cutoff. The DA follows from both. Tables: `compare-939q-940q.tsv`, `compare-939l-940l.tsv`, `compare-939t-940t.tsv`.
- **DIA_LIMPA records:** against KIT-118's records every difference starts at `ruv_covariates` (`compare-m.tsv`: the design's W columns; `compare-mc.tsv`: the adjusted table): R2.
