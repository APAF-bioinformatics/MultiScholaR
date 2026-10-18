# Remediation record for the 2026-10-19 implementation audit of the analysis-record series

Audit: [AUDIT-20261019-implementation.md](AUDIT-20261019-implementation.md) (post-hoc, audit only; 3 high, 8 medium, 5 low). This record maps each finding to its owner and is filled with closure evidence as the owning tickets complete. The audit file is not edited.

| Finding | Owner | Status | Closure evidence |
|---|---|---|---|
| IM-01 fixed constants dropped; `@args` dumped for IQ, P14, design | KIT-038 item 1 | open | |
| IM-02 realised RUV outcome absent; NA labelled `user_set` | KIT-038 item 2 | open | |
| IM-03 provenance vocabulary unchecked | KIT-038 item 3 | open | |
| IM-04 database identity (name, release, date, entries) | user decision (paper): not in the series | deferred | |
| IM-05 upstream search provenance and data accession | user decision (paper): not in the series | deferred | |
| IM-06 report prints `NA` software version | KIT-039 item 1 | open | |
| IM-07 SI in kit vocabulary; no parameter labels | KIT-039 item 2 | open | |
| IM-08 fixed interpretive prose in the templates | user decision (release); outside the series | deferred | |
| IM-09 build identity and numerical environment | KIT-038 item 4 | open | |
| IM-10 claim scope (same build and host, from the DIA-NN report onward) | manuscript sentence (deferred until the series is done) | deferred | |
| IM-11 method citations, `inst/CITATION` | user decision (release) | deferred | |
| IM-12 design not SDRF | later gate; none | noted | |
| IM-13 no machine-readable JSON Schema | none scheduled | noted | |
| IM-14 bundle only after "Generate Report" | none (stated deviation) | noted | |
| IM-15 series-level line-by-line review of the rendered origin record | release gate (KIT-050) | deferred | |
| IM-16 what holds | — | — | — |

Stage 2 adds LFQ/TMT records (KIT-020, 043) and Limpa records (KIT-046) on the same builder, so IM-01, 02, 03 and 09 are repaired before those tickets (KIT-038 precedes them in the index).
