# KIT-079 outcome: partitioned omics steps, and the composite's cheaper read

Status 2026-11-03: **completed**. Stage 3b, ticket 8 of 10. The intensity filter, log2 and duplicate resolution have partitioned variants, exact against `whole`, on both omics; `whole` stays the default until KIT-080. By user decision the composite QC figure also reads its PNGs as native rasters, which took the slice's peak below the old in-process build's. No memory claim is made here; that is KIT-080's.

- **Commits** on `performance-v2`: code and tests 4ced76c5, records ca848dc5, bindings eb7d8a27, evidence after them.
- **Base:** dc970263 (KIT-078 closed).
- **Install:** `kit_079` = `core/browser-smoke.VxhycO` (4ced76c5).
- **Ticket:** `tickets/kit/KIT-079.md`, revision 1.2.

## What changed

1. **The three kit gaps,** each failing first ([gaps-fail-first.txt](gaps-fail-first.txt)):
   - `openWorkflowTableGroups` reads more than 64 columns, in blocks;
   - `saveWorkflowPartitionedTable` writes one table in batches;
   - a list slot's column selection can differ by element.
2. **The variants.** `runProjected` hands a `partitioned` request to the declaration's own function, which hydrates no assay. The kernel on the object with zero-row assays supplies the recorded arguments and the design; each assay's table is then written from its pages or groups.
   - **Intensity filter:** per assay, the percentile over every intensity once, then the kernel's helper per page.
   - **Log2:** the unchanged kernel per page.
   - **Duplicates:** an assay with no duplicate comes back unchanged in source order. One with any is resolved by batches of whole groups in key order with every numeric column double, which is what the kernel does to the whole assay. That order is proved by the awkward-assay test, with mixed-case and missing ids.
   - **ITSD:** `whole` only. Its kernel computes the cross-sample mean inside; recorded.
3. **The composite's read (user decision).** The two generators' default reader is `png::readPNG(native = TRUE)`, the resolution unchanged. The figure's store identity comes from the per-assay files' bytes, so it does not move.

## Results

| Gate | Evidence | Result |
|---|---|---|
| Variants exact, fixtures | `test-workflow-partitions.R` | Exact on committed state, value, recorded arguments, dimensions and provenance: 300-sample wide lipid assay; uneven metab assays (12,700 and 28 rows, so more than one page and batch); integer columns, mixed-case and missing ids, all-missing and empty assays on both omics. |
| Variants exact, screen corpora | [lipid](partitioned-exact-079l.txt), [metab](partitioned-exact-079m.txt) | Each variant on a copy of the finished member against the state it committed: assays, design, arguments, dimensions, value. All exact. |
| The worker passes the variant | same test file | A qualified request runs `partitioned` through `workflowWorkerMain`; `whole` stays every default; ITSD has no variant. |
| Composite, real files | [composite-reader.txt](composite-reader.txt) | PNG byte-identical on both omics, PDF page byte-identical at 100 dpi. Peak in a fresh process 2.02 → 0.87 and 2.03 → 0.87 GiB. |
| Journeys (default `whole`, on the worker) | [079l](kit079-kit-journey-079l-stage.json), [079m](kit079-kit-journey-079m-stage.json) | Accepted exact, 15 dispatches each, gaps under 0.1 s. |
| Member records | [079l](member-record-replay-079l.json), [079m](member-record-replay-079m.json) | Replay exact, 10 states; tampered input refused. |
| Oracles and proteomics records | [exactness-gates.txt](exactness-gates.txt) | `lipid79`, `metab79` identical to KIT-078's. LFQ, TMT, DIA_LIMPA identical to kit_065. |
| KIT-028 replay | [replay-check-kit079.json](replay-check-kit079.json) | Exact, 16 states. |
| Proteomics on the chain | [gate-closure.txt](gate-closure.txt) | `kit-gui-079` observed: 13 of 13 dispatches, 0 mismatches against `kit-gui-078`. |
| Figures | [079l](figure-compare-079l.txt), [079m](figure-compare-079m.txt) | 99 of 99 and 91 of 91 stored, rooted assets; the tool's `RESULT differs` is the 12 PCA products only. |
| Lifecycle members | [lipid](kit079-kit-lifecycle-lipid-stage.json), [metab](kit079-kit-lifecycle-metab-stage.json) | Observed. |
| Tests and check | [affected-tests.txt](affected-tests.txt), [check.log](check.log) | 174 files, 15,460 pass; the 7 failures were my pins and one assertion, updated and passing alone. Check at baseline. No full suite (next at KIT-080). |

## Memory and time (reported; no claim)

| GiB anonymous, lifecycle | kit_077 in process | kit_078 worker | **kit_079** |
|---|---|---|---|
| Lipid whole slice | 2.50 | 2.90 | **1.86** |
| Metab whole slice | 2.52 | 2.91 | **1.80** |
| Lipid total peak | 2.62 | 2.99 | **1.96** |
| Metab total peak | 2.63 | 3.00 | **1.89** |
| Normalisation (lipid / metab) | 2.50 / 2.52 | 2.90 / 2.91 | **1.78 / 1.76** |
| DA (lipid / metab) | 1.88 / 1.89 | 1.85 / 1.81 | **1.86 / 1.80** |
| Intensity filter (lipid and metab) | 0.98 | 1.24 | **1.24** |
| Duplicates (lipid and metab) | 1.00 | 1.27 | **1.27** |

- **The peak moved.** It is now the DA step, whose worker's publication sits beside the application; normalisation fell below it. On these corpora the whole slice is below the in-process build by 0.64 to 0.72 GiB.
- **These corpora are small.** About 200 to 300 features per assay, so the peaks are mostly loaded packages and figures rather than matrices. Real data can behave differently, and that is what KIT-080 measures.
- **For KIT-080's rules:** on the screen corpus `whole_slice_anon` would no longer refuse, but its per-step guard (≤ 1.02) would: the intensity and duplicate steps rose from about 0.98 and 1.00 to 1.24 and 1.27 GiB, the worker beside the application. Protocol v1.6 is KIT-081's to freeze on the user's go; this is the question to settle before then. Nothing was changed for it here.
- **Time.** Journeys 148 and 146 s, against kit_077's 277 and 265 s; normalisation 56 s against 117 s; the pause before correlation 6 s against 50 s. These members ran on mains in the `low-power` profile ([power-state.txt](power-state.txt)), and kit_077's power state was not recorded, so this is reported, not claimed. Why the pause shrank was not investigated.
- **Disk.** The lifecycle project is 110 and 100 MiB, against 196 and 188 at KIT-078 and 114 and 111 at KIT-077. The composite's `.rds` went from about 41 to 11 MB.

## Deviations

1. **Two mechanisms in one ticket:** the user's decision, revision 1.2.
2. **`whole` and `partitioned` share an identity.** The variant is not in the declared parameters, so one variant's saved product is the other's hit: the contract is that they are exact.
3. **Batch sizes.** The row streamer's page (10,000 rows) and the cursor's delivery bound are fixed, so "every tested batch size" is the page and batch counts the fixtures reach. No size parameter was added for a test.

## Limits

- **The partitioned variants have not run in a journey or been measured.** The default is `whole`, and KIT-081 makes their operation members runnable. Their benefit may be nil at real sizes; the intensity filter's pre-pass holds every intensity of the assay, about the matrix's size in doubles, and exact quantiles need them.
- **Screen corpora are small** (about 200 to 300 rows per assay), one page. The larger exactness cases are the synthetic fixtures.
- **Proteomics reads its composite the old way** (`R/mod_prot_norm_qc_support_helpers.R:225`). Left alone as scoped; the same change would apply there.
- **Earlier projects:** a composite already stored stays stored, and its bytes are the same.
- **Not inspected by eye:** the composite. PNG bytes and a PDF page were compared instead.

## Rollback

`git revert` of 4ced76c5. The default is unchanged, so no stored state depends on a variant unless KIT-080 changes it.

## Handoff

KIT-081: three partitioned steps per omic to make runnable as operation members, and the per-step guard question above. KIT-080: the defaults, on real data.
