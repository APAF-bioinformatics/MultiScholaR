# KIT-175 outcome: stopped at gate F3 (case ii), for the owner

Run: 45118859939 at commit 4546b50b (script only; no package code changed on this branch). Evidence and the
questions for the owner: `feasibility.md`; per-request counts `f3-requests.tsv`; lists `f3-lists.tsv`.

Observed: on 120 real requests (2,000 human reviewed, 1,000 E. coli K-12 accessions, sent exactly as the base code
sends them), `X-Total-Results` exceeded the rows returned on 115, while every requested accession was returned on
the first page and no returned entry was outside the request.

Inferred: item 1a's truncation rule as written would mark almost every real request failed and drop its rows.

| Criterion | State |
|---|---|
| A1 | F3 counts and case reached: recorded (case ii); F1 not run |
| A2 to A6 | not reached (no test or product commit; the gate stops the work before code) |

## Not checked

- What the extra hits are (later pages were not fetched).
- F1, every test (T1 to T8), the check, the full suite.
