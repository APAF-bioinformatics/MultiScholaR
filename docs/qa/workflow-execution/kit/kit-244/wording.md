# KIT-244 wording: drafts for the owner (VAL-41)

Every text below is a draft. The owner sets the final wording. A change to the refusal or to the Reset notice
updates E5/E5b or E4/E4b of `tests/testthat/test-kit244-limpa-pass-through.R` in the same commit (audit I12).

## 1. The refusal (class `prot_norm_limpa_peptide_input`), as implemented

Shown after "Error in normalization: ":

> on DIA with limpa, "Run Normalization & RUV" needs the protein table that DPC Quantification makes, and the current
> table is still a peptide table. In the Quality Control tab apply the remaining peptide steps and DPC Quantification
> (normally followed by Between-Sample Normalisation and the Peptide Count Filter), then press "Run Normalization &
> RUV" again.

Optional shorter equivalent (AUDIT-KIT-244-impl-20261125, point 2):

> on DIA with limpa, "Run Normalization & RUV" needs a protein table, and the current table is still a peptide table.
> In the Quality Control tab apply the remaining peptide steps and then DPC Quantification (usually followed by
> Between-Sample Normalisation and the Peptide Count Filter), then press "Run Normalization & RUV" again.

## 2. NEWS entry (with the audit's I01 correction), as committed

> Proteomics, DIA with limpa (chosen at the DIA-NN import; not the default route): after the correlation filter,
> pressing 'Run Normalization & RUV' again now keeps the table the first press started from, as it already did after
> the RUV covariates and as the other proteomics routes already do. Before, the correlation filter's result became the
> starting table: re-applying the filter at a lower threshold could not bring back samples a higher threshold had
> removed (the analysis record named the lower threshold), skipping the filter kept those samples removed while
> recording it as skipped, and 'Reset to Pre-Normalization' either changed nothing or, once the filter had been applied
> again with a different result or skipped, refused with a reason that did not apply. At the same or a higher threshold
> the samples kept were already the same. Pressing 'Run Normalization & RUV' before DPC Quantification now refuses with
> a message saying which Quality Control steps to apply first; before, it was accepted and the correlation filter and
> RUV covariates then failed with messages that did not say why (KIT-244).

**Open question for the owner (audit I13):** the entry sits under `### Fixes`. VAL-220 calls part A "a results change
for that sequence only", and the unreleased section also has `### Scientific changes`; KIT-225's entry (VAL-203) set
the precedent of Fixes. Fixes or Scientific changes is the owner's call.
