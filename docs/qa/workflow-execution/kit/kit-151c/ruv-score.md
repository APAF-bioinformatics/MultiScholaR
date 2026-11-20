# KIT-151 c, task 3 (SCI-04, INT-004): group-corrected inputs

Jobs: 8 (CPU: INTEL(R) XEON(R) PLATINUM 8573C, AMD EPYC 9V74 80-Core Processor, AMD EPYC 7763 64-Core Processor); cells: 200, failed: 0. RUVIII_C_Varying against ruv::RUVIII, largest abs difference: 9.01e-10.

| N | arm | S recovered (abs r >= 0.9) | median abs r S | G recovered | group factor: median VE % (mean of views) | median R2 vs group | median -log10 p | S kept in data | G kept | U kept |
|---|---|---|---|---|---|---|---|---|---|---|
| 20 | raw | 11/20 | 0.933 | 13/20 | 14.2 | 0.882 | 9.08 |    1 |    1 |    1 |
| 20 | oracle | 14/20 | 0.956 | 17/20 | 19.3 | 0.945 | 12.3 | 0.954 |    1 | 4.46e-16 |
| 20 | ruv_k1 | 12/20 | 0.932 | 18/20 | 19.3 | 0.992 | 19.4 | 0.849 | 0.999 | 0.0684 |
| 20 | ruv_k2 | 0/20 | 0.0172 | 20/20 | 23.6 | 0.995 | 21.6 | 0.0162 | 0.997 | 0.0102 |
| 20 | ruv_k3 | 0/20 | 0.0176 | 20/20 | 25.2 | 0.996 | 22.3 | 0.0154 | 0.996 | 0.00973 |
| 60 | raw | 12/20 | 0.928 | 12/20 | 12.1 | 0.855 | 25.3 |    1 |    1 |    1 |
| 60 | oracle | 11/20 | 0.911 | 11/20 | 15.9 | 0.838 | 23.9 | 0.981 |    1 | 3.13e-16 |
| 60 | ruv_k1 | 14/20 | 0.948 | 14/20 |   16 | 0.96 | 42.7 | 0.919 | 0.999 | 0.0472 |
| 60 | ruv_k2 | 0/20 | 0.00929 | 20/20 | 18.7 | 0.996 | 71.2 | 0.00922 | 0.998 | 0.00611 |
| 60 | ruv_k3 | 0/20 | 0.00961 | 20/20 | 19.2 | 0.996 | 71.8 | 0.00894 | 0.997 | 0.00605 |

Readings (declared rules):
- N = 20: the oracle arm recovers S in 14 of 20 seeds (< 15): the design cannot show removal at this N.
- N = 60: the oracle arm recovers S in 11 of 20 seeds (< 15): the design cannot show removal at this N.
