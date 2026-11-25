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

## Gate F2: size and time of one mapping (ticket revision 1.10, section 11)

Script `gate-f2.R`. Purpose: one mapping, both directions, of KPC2 against MGH 78578 at the shipped settings, threads
4; ceiling 30 minutes per job; decides the ranking key (5.2) and the determinism pins (5.1), or stops for the owner
(over 10 minutes or 8 GB for one direction).

| Run | Commit | Job minutes | Outcome |
|---|---|---|---|
| 45334031153 | 05a353d1 | 8 | instrument defect: the PRIDE API lookup returned no URL; no result |
| 45335506224 | 4f97dc12 | 25 | instrument defect (fix 1: folder listing): PXD016672's folder (2021/07) holds no FASTA; no result |
| 45338825064 | c67490b5 | 24 | ok (fix 2: also the DDA project PXD020839, where the file is) |

Two instrument runs of the three VAL-54 allows; no expected value, threshold or stop rule changed.

### Observed (run 45338825064, `scripts/1/out/f2/f2.tsv`)

- Engine: `diamond-linux64.tar.gz` 10,221,459 B, SHA-256 `33a765d1…61fbc7` (as 2.2); one member `diamond`; executable
  SHA-256 `d3d72578fa97e87cc0846a06d2f9d7dc481a2106265f9546bc91512972fcd1f1` (identical in all three runs, on AMD EPYC
  9V45, 9V74 and 7763); `diamond version 2.2.8`. Committed in `inst/extdata/diamond-pin.tsv` (Linux row).
- Strain FASTA: `https://ftp.pride.ebi.ac.uk/pride/data/archive/2021/07/PXD020839/BPA_SEP_KPC2_Contaminants.fasta`,
  15,817,761 B (the size in 2.2), SHA-256 `3139e262c1131121ec91ec5980dfc3b117263ad88fe6b687703ebd103bd288ca`; 26,286
  entries: 5,454 `>ref|` (all distinct sequences), 20,585 UniProt-style, 247 other.
- Reference: UP000000265 (5,126 entries, SHA-256 as `hashes.tsv`).
- Runner: AMD EPYC 7763, `avx512f` absent, 4 cores.
- Times and memory (threads 4; `/usr/bin/time -v`): makedb 0.1 s and 0.0 s; forward 19.1 s, peak RSS 173.0 MB, 41,399
  rows, 4,740 queries with a hit; reverse 18.5 s, peak RSS 168.4 MB, 41,403 rows, 4,798 queries with a hit; each
  output 4.0 MB.
- The raw `score` field is an integer on every row of both tables.

### Outcome

Pass: far below 10 minutes and 8 GB. Ranking key: the raw `score` (5.2). Pins kept as run: `--block-size 0.4
--index-chunks 1 --algo 0` (chosen on time and memory only; both proteomes are one block). The counted run's about 122
searches at about 20 s each are far inside its ceiling, so no job split.

### Inferred, not observed

- Section 4 and 2.2 cite PXD016672 for the KPC2 FASTA; the file is in the DDA project PXD020839 (as the research note's
  DDA row says), with the size 2.2 gives. The same may hold for the other strains' FASTAs (DDA projects PXD020545,
  PXD020791, PXD020863); for the coordinator to correct in the ticket's inputs.

### Not checked

- Run-to-run timing variation (one run); the `score` type under another strain or reference.

## Gate F1: the Windows engine (ticket revision 1.11, section 11)

Job `kit167_engine` of `distributable.yml` (dispatch input `stage = kit167-engine`, windows-2025, no build, no token,
nothing uploaded; results as `KIT167` log lines read with `gh run view <id> --log`).

| Run | Commit | Outcome |
|---|---|---|
| 45353194826 | (packaging commit before the fix) | instrument defect: the container steps could not reach the Docker engine; the executable's SHA-256 and imports were read |
| 45353725258 | (after the fix: the engine started and waited for) | complete |

### Observed (run 45353725258)

- `diamond-windows.zip` 3,366,672 B, SHA-256 `8035ece1…cd288ce1` (as 2.2); `diamond.exe` SHA-256
  `c8bd6efa915af1a09b13c7de53aec6b75e96196da6df0509646d8f83e9605eb1` (committed in `inst/extdata/diamond-pin.tsv`;
  the same value in both runs).
- `dumpbin /dependents`: KERNEL32.dll, MSVCP140.dll, MSVCP140_1.dll, VCRUNTIME140.dll, VCRUNTIME140_1.dll and the
  `api-ms-win-crt-*` set; linker version 14.51.
- Server Core ltsc2025 container (host 10.0.26100), a 100-protein search of the shipped SF370 reference against itself:
  no runtime: registry absent, `diamond version` and the search exit 0xC0000135; the oldest x64 redistributable Microsoft
  serves (14.0.24215, installed with status 0): registry `Installed = 1`, `Bld = 24215`, both still exit 0xC0000135;
  Microsoft's current redistributable (status 0): `Bld = 35211`, both exit 0, 113 hit rows.

### Outcome

The program fails for a missing DLL without the runtime, so the installer check and the refusal ship (decided, pool S6);
the refusal's status 0xC0000135 is the one `checkDiamond()` classifies. `Bld` for the installer check is set to 35211,
the oldest build observed to work (`packaging/electron/installer.nsh`); no build between 24215 and 35211 was tried.
Step (4) of F1 (the installer's registry check against both states) was run as the same 64-bit registry read in
PowerShell in each container state; the NSIS script itself was not executed in the container (it runs only inside the
installer).

### Not checked

- Status 0xC0000139 (an entry point missing) was not produced by any state tried.
- A runtime build between 24215 and 35211; Windows client editions.

## Gates F3 and F6 (run 45352823145, `pass2-inputs.R` at e5a258d2; tables in `gates/`)

### F3: the authors' accessions to UniParc

Each `Results<strain>.xlsx` (PRIDE DIA projects PXD016840, PXD016672, PXD016504, PXD016913) was read; the accession
columns found were `sseqid` and `Uniprot` of the sheets `Comp 1` and `AllProteinsNorm`. Each distinct accession was
looked up as `https://rest.uniprot.org/uniprotkb/<accession>.json` (release 2026_04) and its UniParc id read from the
entry (active) or the inactive entry. Observed: B36 2,175 accessions, KPC2 1,927, BPH2900 1,242, HKU419 1,179; every
one resolved to a UPI (share 1.000 each; 3,566 active and 2,957 inactive entries over the four). Outcome: every strain
evaluable in V5 (at least 90 % resolved). Product committed as `references-authors/authors-accession-upi.tsv` (strain,
accession, UniProtKB status, UPI, release; UniProt CC BY 4.0), for V5 as a `repo://` input.
Inferred, not observed: the counts are of distinct accessions, so they are below the ticket's protein counts (B36
2,229 proteins with an authors' hit), as several proteins share a hit.

### F6: the reference order

For each shipped reference, the INSDC assembly UniProt names for the proteome, its NCBI GFF, and the share of entries
whose `reference_locus` (5.7) joins a protein-coding locus by `locus_tag` or `old_locus_tag`:

| Reference | Assembly | Entries | With a locus | Joined | Share | Pass (>= 90 %) |
|---|---|---|---|---|---|---|
| UP000000625 | GCA_000005845.2 | 4,403 | 4,402 | 4,289 | 0.974 | yes |
| UP000008816 | GCA_000013425.1 | 2,889 | 2,852 | 2,851 | 0.987 | yes |
| UP000000265 | GCA_000016305.1 | 5,126 | 5,126 | 5,126 | 1.000 | yes |
| UP000000750 | GCA_000006785.2 | 1,690 | 1,688 | 1,687 | 0.998 | yes |

Outcome: F6 passes on all four, so V9 carries its criterion on every strain (D1). The GFF SHA-256 values are in
`gates/F6-reference-order.tsv` for the four `inputs.tsv` rows.

### Not checked

- F3: whether the authors' `sseqid` and `Uniprot` columns ever disagree for one protein (both were taken).

## Pass 2 inputs: where the search FASTAs are (ruling (a)), gate F4, the GEO cut, the demo tables (run 45386584323)

`pass2-inputs.R`; runs 45352823145 (F3, F6 produced; A and CUT failed in the script's mechanics), 45372182412
(instrument run 2: timed out in A's per-file HEAD requests, no result), 45386584323 (instrument run 3: complete).

### Observed

- Projects listed (`gates/A-project-files.tsv`): the DIA projects PXD016840, PXD016672, PXD016504, PXD016913 hold the
  `Results<strain>.xlsx`; the DDA projects PXD020839, PXD020791, PXD020863 hold the `BPA_SEP_<strain>_Contaminants.fasta`
  files; PXD020545 (B36) holds `Ecoli_B36_MQtxt.zip` and no FASTA. Downloaded (`gates/A-strain-fasta.tsv`): KPC2
  15,817,761 B `3139e262c1131121ec91ec5980dfc3b117263ad88fe6b687703ebd103bd288ca`; BPH2900 14,666,884 B
  `0ab6a83cd0f15bb201c807df28a85088fd0860fbea1dfdcc6b1e81dc4b8b9abe`; HKU419 14,375,644 B
  `a21bd43f56e151da34ae8b36247a9f8358b78ef956ad2aa8271fcc3345c21ab2` (sizes as 2.2).
- F4: the zip central directory of `Ecoli_B36_MQtxt.zip`, read by an HTTP range request, lists 19 entries, all
  `txt/*.txt` (`gates/F4-b36-zip-entries.tsv`), no FASTA; B36's protein FASTA built with NCBI efetch from the 4,922
  distinct WP_ accessions of its GEO GFF, all returned (1,950,215 B, `8faeed71d51b4b19c9fefc56e474b35c4b16e7a644748e6f3fe18cdeecffd776`;
  `inputs/B36-efetch.fasta`, the versions in `inputs/B36-efetch-versions.tsv`). The four GEO annotations
  (`gates/F4-geo-annotation.tsv`, `gates/F4-geo-attributes.tsv`): GFF3 with `protein_id` on most CDS rows (B36 4,995 of
  5,194); the four GTF files carry no `protein_id` on any CDS row, so the parser refuses them (`reference_annotation_no_protein_id`,
  5.5) and the GFF3 is the file to use.
- CUT: 54 B36 loci, 113 feature lines, expectations by awk (`tests/testdata/reference-mapping/`); no GTF cut.
- DEMO, the real-engine lane (`gates/DEMO-summary.tsv`), DIAMOND v2.2.8 from the pin, threads 4:

| Strain | Reference | Seconds | Proteins | one_to_one | many_to_one | one_way | tie | below_threshold | no_hit | Median identity |
|---|---|---|---|---|---|---|---|---|---|---|
| B36 | UP000000625 | 34.9 | 4,922 | 3,664 | 269 | 22 | 6 | 230 | 731 | 99.3 |
| KPC2 | UP000000265 | 43.2 | 5,454 | 4,456 | 106 | 4 | 4 | 170 | 714 | 100.0 |
| BPH2900 | UP000008816 | 21.1 | 2,711 | 2,447 | 34 | 6 | 0 | 88 | 136 | 99.5 |
| HKU419 | UP000000750 | 13.9 | 1,861 | 1,574 | 35 | 4 | 1 | 55 | 192 | 100.0 |

  KPC2 mapped again at threads 1: the same mapping SHA-256 and decision-table SHA-256 as at threads 4.

### Outcome

F4 gives B36's FASTA by its second route (every WP_ accession of the GEO GFF returned by efetch) and the GTF attributes. Not checked here: that every DDA/DIA-measured B36 WP_ is among them (F4's stop rule; V4 counts it). The demo
tables are committed as `inst/extdata/reference_mapping/demo/<strain>.tsv.gz` (gzip without a time stamp). These are
descriptive outputs of the shipped configuration (V6's territory), not validation results.

### Not checked

- That the B36 efetch FASTA equals the FASTA the authors searched (that file was not deposited); it is RefSeq's current
  version of each accession.
