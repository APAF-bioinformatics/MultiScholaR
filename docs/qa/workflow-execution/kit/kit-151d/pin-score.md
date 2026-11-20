# KIT-151 d, Q3: the pin on an Intel processor

Jobs: 20 (Intel 8, AMD 12); ran fully: 12. Models: AMD EPYC 7763 64-Core Processor (1); AMD EPYC 9V45 96-Core Processor (5); AMD EPYC 9V74 80-Core Processor (6); Intel(R) Xeon(R) 6973P-C (1); INTEL(R) XEON(R) PLATINUM 8573C (7).
Instrumentation: k151d-pin--06 PN: numpy AVX512 still on; k151d-pin--15 PN: numpy AVX512 still on; k151d-pin--18 PN: numpy AVX512 still on.

| condition | cases | cases with one class | Intel-AMD pairs not identical | Intel-AMD pairs with different iterations | largest Intel-AMD abs dZ (raw) | lowest Intel-AMD matched abs r | same-model pairs not identical |
|---|---|---|---|---|---|---|---|
| U | 10 | 0 | 320/320 | 20 | 4.97 | 0.000598299 | 0/240 |
| P | 10 | 10 | 0/320 | 0 |    0 |       1 | 0/240 |
| PN | 8 | 8 | 0/256 | 0 |    0 |       1 | 0/192 |

- U: Intel and AMD give different models
- P: bit-identical across Intel and AMD
- PN: bit-identical across Intel and AMD

Unpinned, by vendor (distinct iteration counts and planted abs r): see pin-unpinned-by-vendor.tsv.

History (instrumentation, not a gate): k151c pin P: 96 of 96 rows match a committed row; k151b seed 10 (AMD EPYC 7763 64-Core Processor) U: 11 of 12 rows match a committed row; k151b seed 23 (INTEL(R) XEON(R) PLATINUM 8573C) U: 11 of 12 rows match a committed row; k151c pin U: 32 of 32 rows match a committed row
(k151c unpinned rows are compared on the same CPU model only, pinned rows on any; k151b rows on any CPU model, with the CPU of the k151b job named; per job in pin-history.tsv.)
