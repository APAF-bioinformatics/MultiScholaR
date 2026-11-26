R1c (45797066259): reported, not adjusted: the P4 head-rule expectation (`refused: ...`) passes on the old code for prot_tmt,
prot_dia and prot_dia_limpa, whose lanes stop at the Accession Cleanup step (the head is a QC state, not one of the DA's
states, so the rule already refuses). It fails by value for prot_lfq (lane until "da"), where the head is final before the cut.
The rule is therefore only discriminating on the lfq lane.
