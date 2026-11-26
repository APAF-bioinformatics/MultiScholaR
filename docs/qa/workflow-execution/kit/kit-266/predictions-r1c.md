# R1c declaration (head check of the new file, before R2)

Purpose: M1's path control failed (see m1/declaration.md); learn from the in-process double whether the worker branch finds
its files saved, before R2 is declared. Binding: the fix commit 7da7df16 (head). The new test file only, base none, one job,
ceiling 45 job-minutes, expected 15. Result-producing: a failure here is a finding that I correct before R2 (R2 not yet sent;
the one correction pass of R1/R2/R3/M1 is not used by it).
Prediction (from the design): every E and G row passes. Rows I regard as most at risk, by reading: E7 (Reset mechanics), E9/E9b
(snapshot tampering), E11 (export/import), E12 (c), E20 (chained order). A failure of E16, E17, E18, E19-E21 on the worker count
(application builds) would show the identities not found saved, i.e. the cause of M1's path control.
