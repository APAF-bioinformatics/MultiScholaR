# KIT-167 feasibility record

## Step 1: capture of the UniProt references (A0; ticket revision 1.9, section 4, "Pinning")

Script: `docs/qa/workflow-execution/kit/kit-167/capture-references.R` at commit
bf56a4062893514a8b2706b000e48ddc67d80911 (branch `kit-167-reference-mapping`). Requests on `run-tests`
(commit a92716ab): `k167-capture-bf56a406-a` and `k167-capture-bf56a406-b`, pushed together on 2026-11-23.
Bounds declared before sending: ceiling 45 min per job, 6.0 core-hours for both; stop rules as section 4.

| Job | Run id | Scripts job wall (UTC) | Minutes | Script seconds | Result |
|---|---|---|---|---|---|
| a | 45127082405 | 06:52:41 to 07:01:26 | 8.8 | 70 | `ok` |
| b | 45127086068 | 06:52:44 to 07:04:41 | 12.0 | 38 | `ok` |

Core-hours (scripts jobs, 4 cores): (8.8 + 12.0) / 60 x 4 = 1.4, against the expected 2.7 and the ceiling 6.0.

### Observed

`scripts.tsv`, job a:

```
index	path	exit	seconds	status
1	docs/qa/workflow-execution/kit/kit-167/capture-references.R	0	70	ok
```

`scripts.tsv`, job b:

```
index	path	exit	seconds	status
1	docs/qa/workflow-execution/kit/kit-167/capture-references.R	0	38	ok
```

- Neither artifact holds `withheld.txt`; the scan pre-count (`scan-precount.tsv`, not committed) found no match in either job.
- `checks.tsv`: every check passed in both jobs (44 each; job a's copy in `capture/checks.tsv`).
- Comparison across the jobs, on the laptop: `diff` of the two `hashes.tsv`: no difference; `diff` of the two
  `release.tsv`: no difference.
- Release pinned: **2026_04**, read from `current_release/RELEASE.metalink` (the `reldate.txt` fallback was not
  needed); all eight UniProtKB stream replies and the three UniParc replies carry `X-UniProt-Release: 2026_04`.
  The release is the one expected, so the counts of 2.2 and the shares of 5.7 need no re-reading.
- Citation: job a's `citation.tsv` reads `confirmed` (job b's also): The UniProt Consortium, "UniProt: the
  Universal Protein Knowledgebase in 2025", Nucleic Acids Res 53:D609-D617 (2025), doi:10.1093/nar/gkae1010,
  read from https://rest.uniprot.org/help/publications (the JSON reply's `content` field).
- Entry counts equal the counts of section 2.2: K-12 4,403; MGH 78578 5,126; NCTC 8325 2,889; SF370 1,690;
  UniParc UP000001734 5,739, UP000000530 2,721, UP000002436 1,993. In each reference the FASTA's accessions
  equal both tables' accessions.

`hashes.tsv` (job a, identical in job b):

```
file	bytes	sha256	entries
UP000000265.fasta	2469712	6f307f8d5de3b81509b13a98eb2e8141c17f52ce6c6dc87f21a2dc1088cf5981	5126
UP000000265_go.tsv	1147948	8598415a91efc926a99351e4842c7894ce5d2d43554e22175a5ac84872b8e63d	5126
UP000000265_xref.tsv	232084	f07f03123a54606e33ae05cabc3d7dc49342f7d425aabdbc2b62de3e767c170f	5126
UP000000530_uniparc.fasta	883830	37f073db90c795bf67bcd0fa7b2076bfdbd6290a949846f49eed6dab456d19c5	2721
UP000000625.fasta	1891162	fd3faad09b87e46d16d1dac0f21b1f6ae3c92c36b6269e886a114c4e731a0bfb	4403
UP000000625_go.tsv	1047533	2fafcd7a700fc4b27a8f4b189875e62b9017282dfc9a663424c5ea6652030dd4	4403
UP000000625_xref.tsv	333013	2b28e804f57c6d52e4d9de5f0babb78ed1fae2b84c5a8be2be0d75be07f1d322	4403
UP000000750.fasta	735243	a777a5de83d2d065611bd514353f646d5ae4d5f1105c54b406e128b9e3f8f763	1690
UP000000750_go.tsv	310938	c828a012c39d412a7fe83960501fe5c7e01f2ab6e3de383831b2431c9f17b45f	1690
UP000000750_xref.tsv	70977	b734bdc2b8f1d657f5f44370ba4fa727c151c72386a95a8a3b57cc087f794063	1690
UP000001734_uniparc.fasta	1910385	01670c0538461898b0206c2bb4849a4dde523fca904301824b420583f264f15b	5739
UP000002436_uniparc.fasta	631889	b24026bbff1c398cd808d2db57330a8a581351669b563ae8520e872dcf52a4d6	1993
UP000008816.fasta	1227493	9177de0f9e65c32962d3f751841c96eb156843e36c450da226e62f5fe3330116	2889
UP000008816_go.tsv	578247	eead90181aebc00ff683e7eed6fae55e1807925181307761684d72f29fc4fff7	2889
UP000008816_xref.tsv	175632	75264a9ae2a555e91786c99df5fe6af7d0a7c6b550f808b657425ef497422257	2889
```

### Committed

- `inst/extdata/reference_mapping/references/`: the twelve reference files; `ATTRIBUTION.txt` in
  `inst/extdata/reference_mapping/`.
- `docs/qa/workflow-execution/kit/kit-167/references-authors/`: the three UniParc proteomes and `ATTRIBUTION.txt`.
- `docs/qa/workflow-execution/kit/kit-167/capture/`: job a's `hashes.tsv`, `release.tsv`, `checks.tsv`, `citation.tsv`.
- Both `ATTRIBUTION.txt` files were written by a shell step from `release.tsv` and `citation.tsv`.
- `.gitattributes`: `-text` for both folders, so that a Windows checkout keeps the bytes the hashes describe.
- Every committed file is under 5 MB (largest 2,469,712 bytes) and every path is under 100 bytes.

### Inferred, not observed

- The kingdom folder in which each FASTA was found is not recorded; only the file names are (`checks.tsv`, `fasta_file`).
- `release_date` is `NA`: `RELEASE.metalink` gave no date under the tags the script reads. The date is not
  needed for the pin.

### Not checked

- The `LICENSE` sentence of A0 and 8.1 is not added in this step (outside the step-1 file list given to the implementer).
- The record phase of KIT-159 (`repo://` rows hashing the committed copies) is not run.

## Gate F5: attachment point per route (code reading, ticket revision 1.10, 8.5)

Read on `kit-167-reference-mapping` after the merge d3909a30 (`v0.9-staging` at ea93d948).

### Observed (code)

- Every route has one session owner per omic, a `DIAStateAdapter` made at session start by
  `createWorkflowSessionOwner()` (`R/utils_workflow_dia_sessions.R:885-905`) with `reference_states = TRUE`
  (`:902`); the legacy-session adoption sets the same (`:1473`). DIA, DIA with limpa, LFQ and TMT share the proteomics
  owner (`omic_type = "proteomics"`); transcriptomics has its own (`omic_type = "transcriptomics"`).
- An import-time table attaches to the head through `saveWorkflowOwnedValue(table, label, manager)`
  (`R/utils_workflow_dia_state.R:2219`): a data frame becomes an owned table component and
  `snapshot$gui_references[[label]]` on the published head; it is read back with `resolveWorkflowOwnedValue()`
  (`:2272`). Input files attach by role with `recordWorkflowInputs()` (`:2984`), as the imports register theirs
  (`completeDIAInputFiles()`, `:3044`).
- `recordWorkflowExternalResource()` (`:3014`) records on every route whose owner uses reference states, except inside
  a worker (`MULTISCHOLAR_WORKFLOW_WORKER = 1`), where nothing is recorded; the mapping runs in the controller at
  import, so it records.

### Outcome

Pass on all five routes (DIA, DIA_LIMPA, LFQ, TMT, transcriptomics): the mapping table, its record sheet and the
transcriptomics gene table attach as owned tables under the labels `reference_mapping`,
`reference_mapping_record` and `reference_mapping_genes` (`R/utils_workflow_reference_mapping.R`), and the input files
by role. No route lacks an owner at import, so no route's panel waits.

### Not checked

- That the owned tables survive a session export and re-open on every route (no run; the gui-reference path is the
  one the imports already use).
