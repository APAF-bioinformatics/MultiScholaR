# KIT-209 dated section (2026-11-24): the candidate, gate B, the base's library, measured times, limits and layouts

Committed alone before M1 (ticket 9.1, "Order of commits", item 3). Declaration f6b21504 and its addendum 1 (55fbe19c).

## The release candidate

- **RC = `91e3cdcb5e3a6fa34712a1575ecbfb855ed02cf2`**, the staging head on 2026-11-24 when the coordinator said to
  proceed (after KIT-207's merge, 2ead4fec). G0 ran on c03cb02a; between the two, 30 files under `R/`, `inst/` changed
  (KIT-182, KIT-203, KIT-207). A later change under `R/`, `src/`, `inst/` or `DESCRIPTION` before the release commit
  re-measures by the ticket's 9.4 map.
- **Base = `b468a885c28ea7dfe228a0fdaad947825bc4839e`.**
- **RUV panel:** `MultiScholaR.ruv_design_panel` unset, `omicGuiRuvDesignEnabled()` FALSE at G0 (recorded per journey;
  the counted journeys record it at the RC).

## Gate B and the options taken

- B1 to B6 hold for the base on every route (`g0.md`). No option (a) to (d) is taken.
- **A1:** U4 and U5: both builds through the session export; DA and the report base only, labelled "base only; the
  candidate refuses this design (one biological unit per group)".
- **A2:** U11 and U12: RUV skipped on both builds. The base also refuses the bound manual RUV (probe run 45383701044:
  "RUV-III cannot estimate 3 factors from 0 control proteins observed in every sample" on `lfq_300`, "from 2" on
  `tmt_300`; the same counts as the candidate's in G0p).
- **The base's library (A4):** both builds on the candidate's internal dependency versions. Every base figure carries
  "base on the candidate's internal dependency versions (`fastduckdb` 0.2.x, `fast_limpa_dpc` 0ad8c527)"; the two
  upgrades are outside what the comparison shows.

## Measured times and the journey limits (ticket 9.3: 3 × the measured process wall)

| Unit | Measured s (source) | Limit s |
|---|---|---|
| U1 | 125.3 (G0a-2, the larger build) | 376 |
| U2 | 117.2 (G0a-2) | 352 |
| U3 | 77.9 (G0a-2, base) | 234 |
| U4 | 287.4 (G0a-2, base) | 863 |
| U5 | 254.3 (G0a-2, base) | 763 |
| S1, S2 | 67.4, 67.0 (G0a-2) | 203, 202 |
| U6 | 1,555.2 (g0p3, RC) | 4,666 |
| U7 | 1,344.3 (G0p, c03cb02a) | 4,033 |
| U8 | 295.4 (g0p3, RC) | 887 |
| U9 | 2,469.2 (G0p) | 7,408 |
| U10 | 2,566.2 (G0p) | 7,699 |
| U11 | 740.8 through `normalisation_ruv` (the base's probe, the larger of the two builds' measured times to that step; the skipped-RUV journey is otherwise unmeasured) | 2,223 |
| U12 | 738.5 through `normalisation_ruv` (likewise) | 2,216 |
| D1 to D7, X10 | RC only, reported | the job's deadline |

Each limit is cut to the job's remaining time less 12 minutes. A journey stopped at its limit is reported with its class.

## Layouts (measured basis; setup measured 7 to 11 min; the base's size times assumed equal to the candidate's)

| Job | Content | Expected min | Ceiling min |
|---|---|---|---|
| M1 | U1 to U5, 3 pairs each (B C, C B, B C); S1, S2, 1 pair | 101 | 240 |
| J1 | U6, 2 pairs; U8, 2 pairs; D1; D3 | 166 | 240 |
| J2 | U7, 2 pairs; D2 | 128 | 240 |
| J3 | U9, 2 pairs; D4 | 230 | 330 |
| J4 | U10, 2 pairs; D5 | 238 | 330 |
| J5 | U11 (RUV skipped), 2 pairs | about 65 (4 × at most 12.4 + setup and corpus) | 240 |
| J6 | U12 (RUV skipped), 2 pairs | about 65 | 240 |
| D6 | LFQ shipped default (automatic sweep) | 110 (estimate; not probed) | 240 |
| D7 | TMT shipped default (automatic sweep) | 256 (estimate; not probed) | 340 |
| X10 | U10 with `limpa` | 193 (estimate) | 330 |
| M3 | I-M (3 per build), I-L (2 per build), RC and C on KIT-150's lane | 220 (estimate) | 180 / 330 |

Sent in two waves so that at most 8 of these jobs are queued at once: M1, J1 to J6, D6; then D7, X10, M3. P3 (a third
pair) by the ticket's 7.2 rule after M2.

## Corpora (A5)

Both builds read one generated copy per job. The runner's `screen` `report.tsv` and `trx_screen`/`trx_60k_200`
`truth.json` differ from the bound sha256 identically across jobs (deterministic on the runner); `full`'s is compared
across J3, J4 and X10.
