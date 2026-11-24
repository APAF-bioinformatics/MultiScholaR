# KIT-209 outcome so far: G0 (the instrument, gate B, the integration lifecycle, the time probes)

Ticket revision 1.3 (VAL-151: a memory discovery, no gate; VAL-145: in process on GitHub runners). The counted runs
(M1 to M3, D, X10, P3) are not sent: they wait for the release candidate and the dated section. Declaration:
`declaration.md` (f6b21504, committed alone before any run). Runs: `runs.tsv`. G0 detail: `g0.md`.

**State: G0 done.** G0a holds (after two instrument runs under VAL-54), gate B holds for the base on every route,
G0i holds, G0p measured U7, U9 and U10 and found U11 and U12 refused by the package under the bound manual RUV.
G0a' is not run (it needs the coordinator to extend the lane's token handling). Three items go to the coordinator
before M1 and M2 (last section).

## Observed

Labels: OBSERVED = in the returned records named; INFERRED = reasoning from them, not measured.

### Screen units, one run per build (G0a-2, run 45356238103; G0a is not a counted run)

Journey peak = the largest step `VmHWM`; retained heap = the largest R heap after `gc(full = TRUE)` at a step boundary;
wall = baseline to the last step. Per step: `g0a-steps.tsv`.

| Unit | Base: peak MiB / retained MiB / wall s | Candidate: peak / retained / wall | Note |
|---|---|---|---|
| U1 lipid_screen | 1,274.5 (normalisation) / 534.6 / 114.6 | 1,279.1 (normalisation) / 536.4 / 118.9 | |
| U2 metab_screen | 1,273.9 / 529.2 / 110.6 | 1,261.2 / 531.0 / 110.8 | |
| U3 trx_screen | 1,823.3 (ruviii_export) / 937.6 / 72.8 | 1,635.5 (ruviii_export) / 928.6 / 67.3 | candidate lower at normalisation: 1,809.9 against 1,599.6 MiB |
| U4 DIA screen | 1,350.6 (report) / 659.6 / 262.2 | 1,281.6 (rollup) / 337.2 / 242.1, 17 of 19 steps | candidate refuses DA (KIT-128); rollup 1,177.9 against 1,281.6 MiB |
| U5 DIA_LIMPA screen | 1,315.4 (report) / 657.3 / 238.1 | 1,183.2 (rollup) / 335.7 / 218.1, 12 of 14 steps | as U4 |
| S1 lfq_fixture | 969.9 (report) / 628.3 / 62.3 | 969.0 / 627.5 / 60.8 | 5 proteins; fixed overhead only |
| S2 tmt_fixture | 1,015.1 (report) / 628.4 / 61.6 | 968.6 / 627.6 / 61.0 | as S1 |

Base LFQ/TMT figures include the four protein QC steps the base's GUI bypassed by default (`g0.md`, B5). The base runs on
the candidate's internal dependency versions (`fastduckdb` 0.2.x, `fast_limpa_dpc` 0ad8c527).

### Size units, candidate only, one run (G0p, run 45366618275)

| Unit | Outcome | Process wall s | Journey peak MiB (step) | Retained heap max MiB |
|---|---|---|---|---|
| U7 metab_large_missing | completed | 1,344 | 4,232.7 (report); normalisation 3,656.4 | 1,325.7 |
| U9 DIA full | completed | 2,469 | 4,072.1 (rollup); other steps 1,076 to 3,056 | 712.6 |
| U10 DIA_LIMPA full, `fastlimpadpc` | completed | 2,566 | 3,821.6 (rollup); dpc_quant 2,752.3 | 693.8 |
| U11 lfq_300, manual RUV | package error at normalisation_ruv | 718 | 2,035.0 (accession) | 321.6 |
| U12 tmt_300, manual RUV | package error at normalisation_ruv | 728 | 1,990.6 (accession) | 325.8 |

The longest steps: U10 peptide_accession 1,004 s; U9 imputation 639 s; U7 normalisation 715 s; U11/U12 accession 518
and 524 s.

### Integration (G0i, run 45351244978), one run each at S

Stage peaks equal within 21 MiB between the head (`RC`) and build C; training 197.6 s on `RC` against 48.2 s on C.

### Run-to-run spread seen so far (same build; G0 only; G0b was removed by revision 1.3)

| Same-build pair | Journey peak MiB | Retained heap MiB | Wall s |
|---|---|---|---|
| U1 candidate, two runs in one job (G0a-2) | 1,277.9 / 1,279.1 | 536.3 / 536.4 | 119.0 / 118.9 |
| U3 base, G0a-1 against G0a-2 (two jobs) | 1,825.7 / 1,823.3 | 937.2 / 937.6 | 72.7 / 72.8 |
| U3 candidate, two jobs | 1,639.4 / 1,635.5 | 928.3 / 928.6 | 68.9 / 67.3 |
| U4 base, two jobs | 1,336.4 / 1,350.6 | 659.4 / 659.6 | 259.1 / 262.2 |
| U5 base, two jobs | 1,291.3 / 1,315.4 | 657.0 / 657.3 | 233.9 / 238.1 |
| reset control, two jobs | 1,148.2 / 1,147.9 | 198.1 / 198.4 | 2.9 / 2.9 |

OBSERVED: peak ranges 0.3 to 24.1 MiB, retained heap 0.1 to 0.9 MiB, wall within 4 s, at screen size. Nothing is known
yet of the spread at size (one run per unit).

### Corpora

Every generated file matches its bound sha256 except `screen` and `full` `report.tsv` and `trx_screen` `truth.json`
(`g0a-corpora.tsv`, `g0p*-corpora.tsv`): the DIA reports generated on the runner are not byte-identical to the bound
laptop files (INFERRED cause: writer formatting on another library; not checked).

## Against the ticket's estimates (section 9.2, laptop × 1.5)

| Unit | Estimate min | Measured (candidate, runner) min |
|---|---|---|
| U7 | 47 (pre-fix upper value) | 22.4 |
| U9 | 46 | 41.2 |
| U10 | 46 | 42.8 |
| U11, U12 (manual RUV) | 29, 36 | not completed (refused after 12.0 and 12.1 min) |
| screen units | 1 to 5 | 1.0 to 4.8 (process wall) |
| job setup | 35 | 7.1 to 11.3 (snapshot library, both builds, pandoc) |

## DRAFT of the dated section (not committed as the dated section; it needs the RC and the coordinator's decisions)

**Journey time limits (9.3):** 3 × the measured process wall, the larger of the two builds' at screen (G0a-2), the
staging head's at size (G0p), each cut to the job's remaining time less 12 minutes:

| Unit | Measured s | Limit s |
|---|---|---|
| U1 | 125.3 | 376 |
| U2 | 117.2 | 352 |
| U3 | 77.9 | 234 |
| U4 | 287.4 (base) | 863 |
| U5 | 254.3 (base) | 763 |
| S1, S2 | 67.4, 67.0 | 203, 202 |
| U7 | 1,344.3 | 4,033 |
| U9 | 2,469.2 | 7,408 |
| U10 | 2,566.2 | 7,699 |
| U6, U8 | none: not in G0p's list | no limit can be set (a probe is needed) |
| U11, U12 | none: refused | no limit until the binding is decided |

**Layouts and budget restated from measured times** (assuming the base's size times equal the candidate's, which is
unmeasured until M2; setup 11 min, corpus generation 5 to 13 min, measured):

| Job | Content | Expected min (measured basis) | Ticket expected | Ceiling min |
|---|---|---|---|---|
| M1 | U1 to U5 × 3 per build, S1, S2 × 1 | 101 (journeys 87 + setup 11 + corpora 3) | 115 | 240 |
| J2 | U7 × 2 per build (90) + D2 (22) | 128 | 270 | 240 |
| J3 | U9 × 2 per build (165) + D4 (41) | 230 | 265 | 330 |
| J4 | U10 × 2 per build (171) + D5 (43) | 238 | 265 | 330 |
| J1 | U6, U8, D1, D3 | unmeasured (ticket 150) | 150 | 240 |
| J5, J6, D6, D7 | U11, U12 and their shipped defaults | blocked: binding | 151, 179, 110, 256 | |
| X10 | `limpa` on U10 | unmeasured (ticket 193) | 193 | 330 |

M1 + J2 + J3 + J4 on measured times: 697 min = 46.5 core-hours, against the ticket's 915 min = 61.0 for the same rows.
J2 fits one job (R13).

## Runner time used (G0)

| Row | Ticket expected / ceiling min | Used min |
|---|---|---|
| G0a | 75 / 150 | 32.6 (G0a-1) |
| G0 instrument re-runs (VAL-54) | 100 / 540 | 62.6 (G0a-2 52.3, cancelled 0.8, G0a-3 9.5) |
| G0p | 274 / 480 | 155.9 |
| G0i | 40 / 60 | 18.0 |
| G0a' | 60 / 120 | 0 (not run) |
| **Total** | 549 / 1,350 | **269.1 min = 17.9 core-hours** (ticket G0 expected 36.6, ceiling 90.0) |

SIM: not run (removed by revision 1.3). Pushing the ticket branch also started the repository's Module CI on it (not
counted here).

## For the coordinator

1. **U4 and U5 at screen size on the candidate:** KIT-128 refuses DA when every group is one biological unit, which the
   `screen` corpus's design with the bound technical-replicate assignment produces; the base runs it. Decide before M1:
   measure U4/U5 on the candidate without DA and report, or change the screen design binding on both builds, or another
   route.
2. **U11 and U12 (LFQ, TMT, manual RUV):** the bound manual RUV (percentage 5, the tab's k NULL taken as 3) is refused on
   `lfq_300` and `tmt_300` by the candidate ("cannot estimate 3 factors from 0 [2] control proteins observed in every
   sample"). Whether the base refuses is unknown. A binding decision (k, percentage, or the automatic sweep, or a
   complete-case filter) is needed before M2, and goes to the owner with R01's limit.
3. **U6 and U8 have no measured time** (G0p's list did not include them), so 9.3 gives them no limit: a probe (about 30
   min of runner time on the candidate) or a ruling.
4. **G0a'** (the base on its own pins) needs the lane's token handling extended; until then every base figure carries
   the note "base on the candidate's internal dependency versions".
5. **Integration:** the head's MOFA training took about 4 times build C's at S in one run each (equal peaks); for the
   owner as context, to be seen in M3.
6. Unrelated, left alone: the lipid between-sample kernel returns unnormalised data with a warning on an unknown
   method instead of refusing (`R/func_lipid_s4_normalization_methods.R:751-760`, a tryCatch around the switch; seen in
   G0a-2).

## Not checked

- The base at size (no G0 probe binds it); the base's response to U11/U12's manual RUV.
- The run-to-run spread at size.
- Any GUI-path fingerprint record for DIA `screen` on either build (none exists on this image).
- Why the DIA reports' bytes differ from the bound files.
- U6, U8, X10, D1 to D7 times on the runner.

## For the change report

No behaviour changed; measurement only.
