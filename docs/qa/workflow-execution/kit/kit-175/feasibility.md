# KIT-175 feasibility: gate F3 (stopped for the owner)

Ticket revision 1.4. Branch `kit-175-uniprot-annotation`, from `v0.9-staging` 271cd3c4 (a ledger-only commit after
b052b824: `git diff --stat b052b824..271cd3c4` touches `docs/methods-decision-ledger.md` only, so the ticket's line
numbers hold).

## F3: what was run

- Run 45118859939 (`run-tests kit175-f3-gate-4546b50b`), one scripts job on a GitHub-hosted runner (ubuntu-22.04,
  R 4.6.1), network open. Script `f3-uniprot-pages.R` at commit 4546b50b; job wall time 13.9 minutes, about 0.92
  core-hours (ticket R0b: expected 1, ceiling 4).
- Lists: one request each to `https://rest.uniprot.org/uniprotkb/stream` for `proteome:UP000005640 AND reviewed:true`
  (20,414 entries) and `proteome:UP000000625` (4,403 entries), fields `accession,organism_id`; sorted by accession on
  the runner; the first 2,000 (human) and 1,000 (E. coli) taken. All 3,000 match the UniProt accession pattern; the
  organism ids of the taken entries are 9606 and 83333 only (`f3-lists.tsv`). The accession lists are not kept.
- Annotation requests: exactly as `directUniprotDownload` sends them at the base commit (`httr::GET` on
  `/uniprotkb/search`, query `(<25 accessions joined by OR>) AND organism_id:<taxon>`, `format=tsv`, the same ten
  fields, no `size`, time-out 30 s, 1 s apart, one attempt): 80 human requests under 9606, 40 E. coli under 83333.
  Per request the script wrote counts, header presence and the release only (`f3-requests.tsv`, `f3-summary.txt`).

## F3: what was observed

| Measure | Human (80 requests) | E. coli (40 requests) |
|---|---|---|
| status 200 | 80 | 40 |
| `X-Total-Results` present | 80 | 40 |
| requests whose `X-Total-Results` exceeds the rows returned | 76 | 39 |
| hits beyond the rows returned, summed | 9,956 | 568 |
| largest `X-Total-Results` of one request | 559 | 68 |
| `Link` header present (all `rel="next"`) | 76 | 39 |
| rows returned per request | 25 (every request) | 25 (every request) |
| requested accessions not returned | 0 | 0 |
| returned entries outside the requested accessions | 0 | 0 |
| `Retry-After` seen | none (no refusal occurred) | none |
| release header | 2026_04 | 2026_04 |

## Case reached: (ii), stop for the owner

The ticket's stop rule: "Any request with more hits than rows: stop for the owner". 115 of 120 requests had more hits
than rows, so the work stopped here, before any test or product code was committed. Case (iii) was not reached (the
header is present on every TSV answer); case (iv): no refusal occurred, so no `Retry-After` value was seen.

What the numbers show (observed): every one of the 3,000 requested accessions came back on the first page, the first
page held only requested entries, and the extra hits are further entries that the query matches beyond the 25
requested ones, served on later pages (`Link: rel="next"`).

Inferred, not verified (later pages were not fetched): the extra hits are entries that the bare accession terms
match in other fields (for example a secondary accession or a cross-reference that names a requested accession). On
these lists the first page was never short of a requested accession; the default page size (25) equals the chunk
size, and UniProt appears to rank the exact accession matches first.

Consequence for the ticket as written (item 1a): a 200 answer whose `X-Total-Results` exceeds its rows counts as a
failed request and its rows are dropped (C02). On these real lists that rule would mark 96 % of requests failed and
drop their rows although every requested accession was returned: nearly every real retrieval would end
"incomplete" with most of its annotation discarded. Implementing item 1a as written would turn a correct
annotation into a near-empty one.

Questions for the owner (the ticket's stop rule says a query change is a separate ticket, D7):

1. Truncation signal: keep `X-Total-Results > rows` (as written, unusable on real lists), or replace it, for
   example by "a requested accession is absent from the answer and the answer has `Link: rel=\"next\"`", which
   detects a cut that loses a requested entry and did not fire on any of these 120 requests.
2. Query: field-qualified terms (`accession:<acc>`) or a page size of 500 would make the hit count equal the
   requested entries; both change the query, which the ticket puts out of scope (D7).

## Not run

- F1 (R0): not run; it follows the gate.
- The re-run of 2.3's greps on the implementation base: the base differs from b052b824 by a ledger row only; not
  re-run.

# F1 (R0): the consumers without annotation (revision 1.6)

Run 45281585187 (`run-tests k175-r1-def420a1`), scripts job, script `f1-consumers.R` at def420a1 (the tests-only
commit; `R/` equal to the merged base af04046f). 12.3 minutes, about 0.82 core-hours (R0: expected 1, ceiling 4).
Network closed (proxy at a closed port), test mode off, no service called (`clusterProfiler::enricher` wrapped to
record its arguments; `gprofiler2::gost` replaced by a recorder returning no result). Fixture: KIT-173's committed
E. coli project (70 proteins; `tests/testthat/helper-go-propagation.R`). Table: `r1/f1-consumers.tsv`.

## Observed

| Case | Annotation rows | Outcome | Service calls | Identifiers sent per call | Background |
|---|---|---|---|---|---|
| clusterProfiler, complete | 70 | completed | 2 | 25; 10 | 70; 70 |
| clusterProfiler, empty | 0 | error: column `go_id_go_biological_process` not found (unclassed `rlang_error`) | 0 | | |
| clusterProfiler, one chunk of 25 missing | 45 | completed | 2 | 25; 10 | 70; 70 |
| g:Profiler (gost recorded), complete | 70 | completed | 10 | 0 each | 0 each |
| g:Profiler, empty | 0 | completed | 10 | 25 or 10 | 70 |
| g:Profiler, one chunk missing | 45 | completed | 10 | 0 each | 0 each |
| organism filter on, empty (clusterProfiler) | 0 | the 8-second notice "Multi-species filtering enabled but no organism mapping available. Proceeding without filtering.", then the column error | 0 | | |
| organism filter on, complete | 70 | completed, the same notice | 2 | 25; 10 | 70; 70 |
| integration features, complete / empty / one chunk missing | 70 / 0 / 45 | NA gene names 0 / 70 / 25; display name = feature id 0 / 70 / 25 | | | |

DA export: 70 rows in every case; the `gene_name` column blank on 70 of 70 rows with the complete and the partial
table, on 0 of 70 with the empty table. Summary report: not run (Not checked).

## Inferred (not verified)

- The fixture's own protein identifier table carries `gene_names` (`GENE1` ... `GENE70`). With no annotation the DA
  table's `gene_name` comes from that column, so g:Profiler was queried on the fixture's own names; with an annotation
  joined, the two `gene_names` columns collide (`.x`, `.y`) and `gene_name` is blank, so g:Profiler was asked about
  nothing. On data whose identifier table carries no gene names, the empty case sends nothing (section 2.2). The
  collision with an annotation present is a separate defect of the fixture or the DA join, outside this ticket.
- clusterProfiler's background was the 70 DA proteins in every case, also with 25 unannotated: the universe passed is
  the DA background; DOSE then intersects it with the annotated identifiers (section 2.2; not observed here).

## For the owner (D1 premise)

Section 2.2 says the gene-set service is asked about nothing without annotation. F1 shows that this holds only when
the data carry no gene names of their own: with a protein identifier table that names genes, the g:Profiler route ran
without UniProt annotation. The refusal of item 5a (D1, VAL-30) is implemented as decided and also refuses that case.

## Live observation by the implementation auditor (2026-11-24; implementation audit I01)

Release 2026_04, read by the auditor, not on a runner:
1. `(accession:P05067-4 OR accession:P04637) AND organism_id:9606`: status 200, `x-total-results: 10`, no `link`;
   rows `P05067-4`, `P04637`, `P04637-2` to `-9`. The canonical `P05067` is not returned.
2. `(accession:P05067 OR accession:P04637) AND organism_id:9606`: 2 rows, `P04637` and `P05067`, total 2 (F4's shape).
3. The bare query as base sent it, `(P05067-4 OR P04637) AND organism_id:9606`: total 510 and a next link; page one
   starts `P04637`, `P05067` (the isoform term brought back its canonical entry).
F3 and F4 measured canonical accessions only. Consequence: since the correction pass, `directUniprotDownload` removes
isoform suffixes before the query.
