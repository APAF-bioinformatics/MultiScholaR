# KIT-247 withdrawal of VAL-288 (the owner's VAL-331) and C05, declared before sending

- Code: b9695bf8 = `git revert` of 1d0364e0 (VAL-288 (a) and (b) removed: no formula, organism, taxon or FASTA detection,
  no limited question, no limited reset), keeping `omicsDesignCanonicalText` and its two call sites (values read as text,
  blank as missing): without it an unchanged proteomics re-import is judged "changed" and takes the full cut (R1'''
  45676491133 and R1'''b 45688025752, F288b on the tests-only base: head protein_s4_initial, QC pending). I02 does not
  depend on it. C03 and C04 fall away with the withdrawal (they were defects of the limited Saves only).
- C05 (small part, coordinator's approval): the refused confirmed import also puts back `aa_seq_tbl_final`,
  `fasta_metadata` and `uniprot_dat_cln`; the import's inputs (`recordDIAWorkflowInputs`) are recorded after the
  checkpoint. The RDS copies the sidecars write into `source_dir` before the checkpoint stay (a stated limit, a candidate
  ticket).
- Tests first: 70b43e04 (W288a label-free and lipid; W288b) and eab793f1 (C05), both on ac53976a's code.

## The one run (three KIT-247 files, no check; 1 + 1 shard)

Base eab793f1 (R/ identical to ac53976a, the current head the coordinator named; carrying the new rows so that they
show failing first), head b9695bf8. Files: test-kit247-design-resave.R, test-kit247-import-config.R,
test-kit247-resave-no-session.R. Expected about 130 job-minutes (two jobs of about 65: the design-resave file ran 50 to
56 minutes in R1'''b and the final run), 8.7 core-hours, inside the 9.1 left under 190 (180.9 used, measured).

Predictions:
- Base: W288a (prot_lfq, lipid_canonical) FAIL (the formula-only re-save asks); W288b FAILS (the organism-only import
  asks); C05 FAILS (the refused import leaves the folder's FASTA/UniProt tables on the bus and records its design file);
  every other row passes; the two small files pass.
- Head: every row of the three files passes (test-kit247-design-resave.R 188 - 9 + 11 = 190 expectations is not
  predicted exactly; every test passes).
