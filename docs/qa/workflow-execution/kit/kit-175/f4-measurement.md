# KIT-175 measurement F4 (before the owner's decision D8; commits to neither option)

Run 45121855739 (`run-tests kit175-f4-measure-e27c9d69`), script `f4-query-forms.R` at e27c9d69, one scripts job of
20.7 minutes (script 572 s), about 1.38 core-hours. With F3 (0.92), R0b has used 2.30 of its ceiling of 4. Release
2026_04 on every answer. Both forms sent as `directUniprotDownload` sends today (httr, the same fields, no `size`,
30 s, one attempt, 1 s apart); only the terms differ: bare `(A OR B ...) AND organism_id:T` against field-qualified
`(accession:A OR accession:B ...) AND organism_id:T`. Files: `f4a-compare.tsv`, `f4a-forms.tsv`, `f4b-probes.tsv`,
`f4-summary.txt` (counts only).

## Observed

(a) F3's 120 requests (2,000 reviewed human, 1,000 E. coli K-12, current primary accessions):

| | bare | qualified |
|---|---|---|
| status 200 | 120 | 120 |
| `Entry` sets identical between the forms | 120 of 120 | |
| next link | 115 | 0 |
| `X-Total-Results` above rows | 115 | 0 |
| requested accessions returned | 3,000 of 3,000 | 3,000 of 3,000 |

(b) Absent kinds, each request = 20 current reviewed human primaries + 5 probes, under organism 9606. Samples taken
on the runner: 200 of the 1,175,771 secondary-primary pairs of UniProt's `sec_ac.txt` whose primary is a reviewed
human entry (seeded, 175); 200 of the 1,914 deleted accessions of `delac_sp.txt` (seeded); the first 100 E. coli
K-12 accessions by accession.

| Kind, form | probes found as requested | found under the merged primary | not found | fillers returned | next link (requests) | entries outside the request on page one |
|---|---|---|---|---|---|---|
| secondary, bare | 0 | 1 | 199 | 800 of 800 | 40 of 40 | 198 |
| secondary, qualified | 0 | 0 | 200 | 800 of 800 | 0 of 40 | 0 |
| deleted, bare | 0 | 0 | 200 | 800 of 800 | 40 of 40 | 200 |
| deleted, qualified | 0 | 0 | 200 | 800 of 800 | 0 of 40 | 0 |
| E. coli under 9606, bare | 0 | 0 | 100 | 400 of 400 | 20 of 20 | 100 |
| E. coli under 9606, qualified | 0 | 0 | 100 | 400 of 400 | 0 of 20 | 0 |

Under the bare form every page one held 25 rows: the 20 requested primaries and 5 entries that were not requested
(unrelated hits fill the page up to 25); under the qualified form page one held exactly the 20 primaries.

## Inferred, not verified

- `accession:` does not match secondary accessions: a merged accession is not found under the qualified form.
- Under the bare form the merged primary is probably on a later page (1 of 200 reached page one); later pages were
  not fetched.

## Not anticipated by the register's F4 text

1. The qualified form does not find merged (secondary) accessions at all (0 of 200, no next link), where the
   register expected them to come back under the primary.
2. The bare form does not find them on page one either (1 of 200): with absences in a request, unrelated entries
   take their places on page one, so a requested entry (the merged primary) is lost to a later page. This is the
   one real truncation case F3 could not show.
3. The bare form puts entries that were not requested into the returned rows (198, 200 and 100 above: one per
   absent probe), so today's annotation table can hold unrelated proteins whenever a request has absences.
