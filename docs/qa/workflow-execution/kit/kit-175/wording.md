# KIT-175 wording, as implemented (for the owner's edit before release; VAL-41)

Implemented as proposed in the ticket's section 12 (revision 1.6). `{…}` are counts or values, never identifiers.
Code: `uniprotAnnotationSentence`, `uniprotAnnotationNotice`, `uniprotRecordSummary`, `signalUniprotEmpty`,
`signalUniprotIncomplete` (`R/func_prot_annotation_uniprot_helpers.R`), `workflowEnrichmentAnnotationProblem` and the
refusal in `runWorkflowEnrichment` (`R/utils_workflow_enrichment.R`).

## Design-step notice: first sentence by reason

| Reason | Sentence |
|---|---|
| `request_failed` | Protein annotation was not retrieved: UniProt could not be reached ({detail}). |
| `service_status` | Protein annotation was not retrieved: UniProt answered with an error (status {status}) for {k} of {n} requests. [+ " UniProt stopped answering, so the remaining {r} requests were not sent." when r > 0] |
| `no_rows` | Protein annotation is empty: UniProt has no entry for any of the {N} accessions under organism {taxon}. Check the organism (taxon ID). |
| `no_uniprot_accessions` | Protein annotation was not requested: none of the {N} protein identifiers is a UniProt accession. |
| `no_uniprot_accessions`, N = 0 (D02) | Protein annotation could not be requested at this step: it has no protein identifiers to look up. Save the design again to retrieve the annotation. |
| incomplete | Protein annotation is incomplete: {k} of {n} requests to UniProt failed (status {status}); {r} more were not sent after UniProt stopped answering; {t} answers were cut short with a requested protein missing, so {m} of {N} proteins were not looked up. (each clause only when its count is above 0; k = 0 and r = 0: "Protein annotation is incomplete: {t} of {n} answers from UniProt were cut short with a requested protein missing, so {m} of {N} proteins were not looked up.") |
| scripts-directory copy failed | Protein annotation was retrieved but could not be saved to the project ({detail}). |
| `not_attempted` | Protein annotation was not attempted: {detail}. Details: "no protein identifier column is selected"; "the data have more than 10,000 distinct protein groups, the most the annotation step reads"; "the data have no protein identifiers"; otherwise the error's message, redacted. |
| `invalid_response` | Protein annotation was not retrieved: the answer received was not UniProt's complete table (for example a network sign-in page, or an answer cut short). |
| `id_conversion_failed` | Protein annotation was not retrieved: the Ensembl identifiers could not be converted to UniProt accessions ({detail}). |
| `processing_failed` | Protein annotation was retrieved but could not be processed ({detail}). |
| no reason known | Protein annotation is empty. |

Then: "The analysis can continue. Without annotation the results tables have no gene names, plots are labelled with
accessions, and enrichment cannot be run." (incomplete: "The analysis can continue. The proteins not looked up have no
gene names and are left out of enrichment."), then for `request_failed`, `service_status`, `invalid_response`,
incomplete and `processing_failed`: " To try again, save the design again."

## Condition messages

- "UniProt annotation is empty: {reason in words} (0 of {N} accessions annotated)"
- "UniProt annotation is incomplete: {k} of {n} requests failed, {m} of {N} accessions not looked up" ({k} here is
  `requests_failed`, unsent requests included)

## Enrichment step

- Incomplete annotation (the run proceeds): "Enrichment ran with incomplete protein annotation: {m} of {N} proteins
  have no annotation and were left out." Outcome `annotation_error`: the record summary of the incomplete case.
- An error of the step's own retrieval: the first sentence above for its reason (`creationReason`; an error without a
  recorded reason is `processing_failed`), then the refusal.
- **Added in implementation (not in section 12; for the owner):** a project annotation file that cannot be read:
  with the run refused, "Protein annotation could not be read from the project file ({detail})."; with the run
  proceeding on a retrieval, "The saved protein annotation could not be read ({detail}); it was retrieved again."
  Outcome `annotation_error`: "annotation file could not be read ({detail})".
- Refusal (D1): "Enrichment needs protein annotation and this project has none. Save the design again to retrieve it,
  then run the later steps again."

## Record `query_summary`

- complete (unchanged): "uniprotkb search, organism_id {taxon}"
- transport error: "annotation failed: UniProt could not be reached; 0 of {N} accessions annotated"
- refused: "annotation failed: UniProt refused {k} of {n} requests (status {status}) and returned no entry for the
  others; 0 of {N} accessions annotated" [+ "; {r} requests not sent"]
- unreadable or cut short, no rows: "annotation failed: the answer was not a complete UniProt table; 0 of {N}
  accessions annotated"
- no entry: "no entry found for {N} accessions under organism {taxon}" (release from the header)
- incomplete: "{k} of {n} requests failed (status {status}); {m} of {N} accessions not looked up" [+ " ({t} answers
  cut short with a requested accession missing)"] [+ "; {r} requests not sent"]
- not queried: "not queried: none of {N} identifiers is a UniProt accession"; N = 0 (D02): "not queried: this step had
  no protein identifiers (annotation is retrieved at the design step)"
- "not queried: the Ensembl identifiers could not be converted to UniProt accessions"
- "not queried: annotation was not attempted ({detail})"
- "annotation failed after retrieval: the answer could not be processed"
- "annotation retrieved but not saved to the project file ({detail})"

## Methods sentence (section 5)

"Where the UniProt annotation was incomplete, proteins without annotation were excluded from both the tested set and
the background of the enrichment analysis; their number is reported with each result."
