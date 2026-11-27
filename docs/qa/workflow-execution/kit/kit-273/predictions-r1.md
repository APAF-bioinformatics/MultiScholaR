# KIT-273 R1 predictions (tests-only commit; written before the run is sent)

Purpose: the new file `tests/testthat/test-kit273-limpa-ruv-delivery.R` on the base tree (staging 7b4ae7dad + the file; nothing
under R/), one job, ceiling 45 job-minutes. Decides: A1 (failing rows fail at the stated value; guards pass) and stop rule 1.

Eight test_that blocks; expectations per block (ids are in the messages).

| Block | Rows | Predicted on base |
|---|---|---|
| T1, T2b, T3, T7, T8, T9, T13 covariates, in process and worker | setup/DA/status expectations; T13 no file before the DA (pass); T1 file exists (FAIL: absent); T3 no matrix (pass); T2b normalised_results.tsv equals E (pass); T9 (pass); T7 file exists after QC files removed (FAIL: absent) | fails at T1 and T7 existence per path; the value rows inside `if (file.exists(...))` do not run on base |
| T2, T3, T7, T8, T13 correction, in process and worker | T2 file exists (FAIL: absent); T7 (FAIL: absent); T13, T3, T9 pass | fails at T2 and T7 per path |
| T4 skip | all pass | passes |
| T5 stale | all pass | passes |
| T6 mode switch | design delivered (FAIL: absent); matrix not (pass) | one failure |
| T10 classic, LFQ, TMT | all pass | passes |
| T11 writer fails | `ruv_delivery` entry (FAIL: absent); design absent, other deliveries (pass) | one failure |
| T12 dropped sample | matrix exists (FAIL: absent) | one failure |

Failing rows by value, never by error (setup runs the unchanged code). Stop rule 1: T3, T4, T5, T9, T10 or T13 failing on base.
If a failure shows an error class instead of a value, it is read from the log first (instrument defect rule, VAL-54).
