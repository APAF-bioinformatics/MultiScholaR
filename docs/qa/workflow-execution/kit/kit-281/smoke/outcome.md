# S09 send, run 45759150472 (package commit a352632a): result against smoke/declaration.md (recorded, nothing repaired, no second send)

| Declared | Observed |
|---|---|
| three parses | met (sim, score, l1 all parse) |
| `SELFTEST: 2007 checks, 0 differ` | NOT met: 2007 checks, 4 differ |
| probe stopped as instrument | met (Stop, instrument) |
| cal105 Calibrated | met |
| seed 9301 task exit 0; E2, E3, E7, E10 3 of 3; marker assertion passed; provenance smoke true | met (exit 0; 3 of 3 each; provenance smoke true, complete) |
| scorer refuses the directory ("smoke flag") | met (exit 1; the driver's check passed) |

The script exited 1 because of the self-test line only (102 s; run 14.1 min created to updated = 0.94 core-hours, ceiling 2.00).

The four differences, all in case `fields`, column `upper_jack` of decision.tsv (rho_g, G1a, G1b, G2): got 0.97004, 0.06692, 0.00528, -0.02900; predicted 0.95127, 0.06667, 0.00333, -0.03333.
Reading (not yet audited): the prediction in `selftest/make-selftest.py` (`predict_decision`) sets the jackknife upper limit equal to the estimate ("degenerate"), which holds only when deleting a replicate leaves every statistic unchanged. In `fields` the cells differ in their rates, so deleting one replicate of one cell changes the pooled statistic and the jackknife SE is not 0. The scorer's values are consistent with a non-zero SE; the prediction (written in Python before the scorer ran) is the defect, an instrument defect of the self-test, not a scorer result. Every bootstrap-limit, state and boundary prediction (including all q_* cases) matched. A repair needs the coordinator's ruling.

## Instrument defect 1 (self-test item), repair after the audit's Re-check 1 (RC1-1 to RC1-4) and the coordinator's ruling
- Defect: `predict_decision` in `selftest/make-selftest.py` set the jackknife upper limit of the decision rows equal to the estimate (SE = 0 assumed). True only when no (cell, replicate) deletion changes the pooled statistic; in `fields` the cells differ. The scorer was right (the auditor reproduced its four values to 17 digits from the case inputs by a separate script; `upper_jack` is reported only and decides nothing).
- Fix: `predict_decision` computes the jackknife by the delete-one formula over the pool's (cell, replicate) units for EVERY case (`pool_jack_stats`, `jack_se_of`, z = inverse normal of 0.975 from the generator's own code); no number taken from the run or the audit. No R script changed.
- Diff of `expected.tsv` (old against new, `python3 -I`): 1,937 lines both; exactly four values change, all `fields` / decision.tsv / `upper_jack`: rho_g 0.9512715977480101 -> 0.9700430222170371; G1a 0.06666666666666667 -> 0.06692495429001503; G1b 0.0033269353806781865 -> 0.005281320797762732; G2 -0.033333333333333326 -> -0.02899713199741668. They equal the auditor's RC1-1 values (0.97004302221703709, 0.066924954290015035, 0.0052813207977627317, -0.028997131997416679) to 16 digits. No other value changed.
- No second parse-and-self-test send is used; L1's own E4 carries the self-test (`SELFTEST: 2007 checks, 0 differ`).
