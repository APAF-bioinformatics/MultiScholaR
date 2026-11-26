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
