# KIT-123: proteomics RUV, k bounded and levels kept (2026-11-16)

Scope: two bugs of the proteomics RUV code (parked row 14), in both proteomics RUV-III computations. R3 (how the DA should account for RUV) moved to KIT-124 (revision 1.1).

**Result: both fixes landed, every moved value is attributed to one of them, and the tree has no regression.** Correctness only: KIT-123 makes no memory claim.

| Item | Change | Commit |
|---|---|---|
| R1 | `findBestNegCtrlPercentage` chooses k among the values up to `max_acceptable_k`; before, the limit applied only to the penalty and a larger k could be returned. | `9e8db5a7` |
| R2 | `ruvIII_C_Varying` for protein data keeps each protein's mean (`restoreRuvFeatureLevels`); `computeRuvCovariates` computes W from controls centred on their means. | `bee0e705` |
| W1 | Reference tests against the `ruv` package (correction equal up to a constant per protein). | `0b6fcb6d` |
| Re-freeze | Two norm test files (`9e9f4f9c`); the DIA_LIMPA entry oracle `prot_dia_limpa.json` (this close-out). | |
| NEWS | Both under "Scientific changes". | `b468a885` |

## Attribution of every moved value

- **Unit fixtures:** every moved frozen value is R2's; R1 moves none. Pre-RUV states, the sweep's k, percentage, controls and scores, and the lipid and metab digests are identical. R2 can change which protein the low-intensity filter after the classic correction removes (fixture: `PG_000019` before, `PG_000053` now).
- **DIA_LIMPA records** (`m/`, `mc/`, `compare-m.tsv`, `compare-mc.tsv`): every difference from `ruv_covariates` on (W columns; the corrected table in correction mode).
- **Journeys on kit_940** (`journeys.md`): 940m and 940r exact against the new records; 940l and 940t identical to 939 (screen journeys skip RUV); 940q against 939q identical up to `ruv_corrected`, then k 45 to 4 (R1) and 2,476 to 2,500 proteins after the correlation filter (R2: the old offsets pushed 24 under the cutoff), DA following.
- **GUI entry lanes:** DIA, TMT and both LFQ lanes identical at every step; DIA_LIMPA identical before `ruv_covariates` and different from it on (`final/rerun/oracle-diff.txt`).

## Suite ([final/suite.md](final/suite.md))

Full suite both modes, `R CMD check` (same six findings), seven e2e lanes: no regression. The one failing file (`test-gui-entry-proteomics.R`, DIA_LIMPA lane) failed only at R2's step; its oracle is re-frozen and the file passes, 193 of 193.

## Measurement validity

940r, 940l, 940t and 940q were measured after the switch to the power-saver profile; their sampling gaps (1.06 to 6.68 s) exceed the 0.5 s limit, so no memory figure is taken from them. Their content checks hold. 940m (balanced profile) is valid. Not rerun: the profile would repeat the gaps, and the ticket claims no memory change.

## What stays open

- Lipid and metab keep the level shift (out of scope as written); it goes to KIT-124's decision table.
- Whether the DA should test corrected values, and with which controls and k: KIT-124, up to the user's gate. Publication numbers that use RUV wait for it.
