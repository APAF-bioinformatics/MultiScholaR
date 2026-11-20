# KIT-151 c, task 1 (SCI-06): training seeds

Jobs: 1; CPU: AMD EPYC 9V74 80-Core Processor; cells: 80, failed: 0.

| data set | arm | mode | tier | identical pairs | min pair-min abs r (Z) | min (W) | max abs dZ | iterations | ELBO range (rel) | start: max abs d | PCA solver |
|---|---|---|---|---|---|---|---|---|---|---|---|
| n20weak@7 | scale_views | fast | D different solutions | 0/45 | 0.258758 | 0.258325 | 1.08 | 21,26,31 | 21.2 (2.73e-05) |    2 | randomized (20, 22300) |
| n20weak@7 | scale_views | slow | D different solutions | 0/45 | 0.0777853 | 0.0746876 | 1.53 | 516,526,536,541,561,596,641,646,746,756 | 4.67 (6.03e-06) |    2 | randomized (20, 22300) |
| n20weak@7 | cap1000 | fast | D different solutions | 0/45 | 0.547353 | 0.51565 | 1.12 | 76,81 | 5.42 (6.83e-05) |    2 | randomized (20, 2300) |
| n20weak@7 | cap1000 | slow | C same solution | 0/45 | 0.999008 | 0.998721 | 0.137 | 546,551,556,561,571,576 | 0.0161 (2.03e-07) |    2 | randomized (20, 2300) |
| n60strong@1 | scale_views | fast | D different solutions | 0/45 | 0.00201829 | 0.00226912 | 1.33 | 26,31,41,46,51 |  331 (0.000163) |    2 | randomized (60, 22300) |
| n60strong@1 | scale_views | slow | D different solutions | 0/45 | 0.0181378 | 0.00208048 | 2.41 | 801,836,951,956,961,966,971 |  259 (0.000127) |    2 | randomized (60, 22300) |
| n60strong@1 | cap1000 | fast | D different solutions | 0/45 | 0.111776 | 0.103792 | 3.57 | 76,81,86,91 | 15.9 (7.9e-05) |    2 | randomized (60, 2300) |
| n60strong@1 | cap1000 | slow | D different solutions | 0/45 | 0.0774322 | 0.00347531 |  3.3 | 796,806,811,831,1001 | 26.5 (0.000132) |    2 | randomized (60, 2300) |

Tiers (declared): A every pair bit-identical; B max abs factor difference after matching <= 1e-8; C every matched factor |r| >= 0.99; D otherwise.

Reading: seeds differ in 8 of 8 cells (tiers C, D); magnitudes in the table and seed-factors.tsv.
