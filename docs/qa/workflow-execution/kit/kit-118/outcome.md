# KIT-118: the combined re-record (2026-11-15)

Scope: the executive calls of 2026-11-15 ([../parked-20261113.md](../parked-20261113.md)), four mechanisms landed one at a time. Decisions left for the user are rows 10 to 15 of that file's last section.

| Mechanism | Result | Commits |
|---|---|---|
| M1: DIA-NN reports an exact 0 for "not detected"; it is read as missing | Landed. Results as declared. **The import memory gate was not met** (row 15). | `91e81f91` (kit_935) |
| M2: RUV negative controls chosen among complete proteins | **Not landed.** The rule written before the run failed (rows 10, 11). | measurement `ab28eea2` |
| M3: lipid and metab zeros read as missing, with a warning and a notice | Landed. | `7cfdfe7e` |
| M4: P-09, P-14, P-18, P-10, the DIA-NN `Precursor.Normalised` checkbox removed | Landed. P-27 and KIT-115 item 10 not done (rows 12, 13). | `f1d39e8b`, `5b6cc3b9`, `7a804d01`, `bdbc0286`, `7c6f2f9c` |

## M1

- **DIA_LIMPA:** exact against records made on kit_935 (935m covariates, 935r correction).
- **Other routes:** exact against today's oracles (935l, t, a, b, c).
- **DIA (IQ) full journey, before and after, store against store** (`m1-dia/`): 107,704 precursor cells became missing. No value and no row changed in the peptide tables. Protein tables before RUV change by a median 0.004 log2 (largest 1.78); 279 cells became missing and 15 became observed. No protein is significant in either build.
- **Not aligned:** the `imputed`, `ruv_corrected` and `correlation_filtered` tables could not be compared cell by cell. After RUV, 19 proteins differ, and RUV-corrected values carry arbitrary per-protein offsets on this corpus (k = 44). That is today's behaviour, not M1's (row 14).
- **Memory:** the application's import peak is 2.42 GiB against 2.30 GiB, twice (+0.12; the gate allowed +0.10). After design it is no higher. The declared consequence of a miss was to move the rule to the store write or to park and revert. Moving it is not contained, and the code was not reverted: that departure is the user's to confirm (row 15).

## M2

`kit-114/sim/d_ruv_controls.tsv`: at n = 12 with 50% missing, the achieved false discovery rate is 0.066 with today's selection and 0.120 with complete-case controls, against a limit of 0.088. The switch is kept as a patch (`kit-114/sim/d_ruv_controls-switch.patch`) and is not in the package. The same simulator shows RUV-III correction mode at 0.22 to 0.48 against 0.06 to 0.09 without RUV (row 11).

## M3 and M4

- No corpus holds a zero in a sample column, so no journey result moved. Six frozen digests in `test-workflow-omics-import.R` were re-frozen, each with an attribution comment.
- P-18 changes the metab normalisation export to the design's samples (one capture re-frozen). The checkbox removal changes one UI hash in two entry oracles.
- kit_936 (`da3320b6`): 936m, l, t, a, b, c, s, v all accepted (`exit/journeys/member-results.txt`).

## Exit

Full suite in both modes, `R CMD check` and seven e2e lanes at `e94b92b1`: no regression ([exit/suite.md](exit/suite.md)). One new test errored in check mode and was corrected (`96129e0c`).

## Acceptance

Met for M3 and M4. M1 is met on results; its import memory gate was missed by 0.02 GiB and the user accepted that on 2026-11-16 (parked row 15). M2 did not land: its question is taken up by KIT-123 (the control rule and the DA model of the classic RUV path).
