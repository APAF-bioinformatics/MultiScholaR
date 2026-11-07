# KIT-101 outcome: the GUI workover's design, contract and behaviour baseline, plus 28 defects found on the way (2026-11-08)

**Done.**
- **Design.** The proteomics verdict, the shared layer's boundary, the omic configuration contract, the core-to-module back-references and the family work map are written in [design.md](design.md).
- **Behaviour baseline.** Ten entry-point lanes are frozen and deterministic: five proteomics routes, lipid, three metab and trx. They pass on the final code, also under the strict-skip sentinel.
- **Test infrastructure.** The sentinel, the per-file baseline runner, the offline UniProt cache, a trx fixture and richer proteomics fixtures.
- **Findings.** The work surfaced 28 defects. The user directed them fixed, queued or parked, and none is left unhandled ([findings.md](findings.md)). Each fix had a test that failed before it.

**Final build:** `kit_101` = `core/browser-smoke.rvol2q` (source `23558bbe`). Two later commits:
- `b91a2b13` (F-28), a one-line change to the UniProt failure path, not installed. It is covered by its unit test, `R CMD check` on `0fc510d9`, and the e2e lanes, which load the source.
- `0fc510d9`, which changes tests only.

## Proof

| Check | Result |
|---|---|
| Entry-point suite, 10 lanes, frozen from `ddf42372` | 356 passed, 0 failed, also under `MULTISCHOLAR_STRICT_SKIPS=1` ([entry-suite.md](entry-suite.md)). On the final code: lipid 35, metab 105, proteomics 187, trx 29, all passing. |
| Screen journeys on `kit_101` | `101a` lipid, `101b` metab and `101c` trx are accepted and exact. `101m` DIA_LIMPA, `101l` LFQ and `101t` TMT were recorded as rejected, then re-verified offline as accepted and exact (0 mismatches; see the section below). |
| e2e browser lanes ([e2e-final.txt](e2e-final.txt), [fixes8/e2e.txt](fixes8/e2e.txt)) | lipid 88/0; metab 190/0; trx 14/0; proteomics DIA 42/0, DIA_LIMPA 12/0 plus 1 skip by design (F-25, KIT-104), enrichment 98/0, TMT/LFQ 153/0. The proteomics lanes make no network calls. |
| Full suite, `load_all`, one process per file ([full-suite-final.tsv](full-suite-final.tsv)) | 381 files, 34,033 passed. Three items: 1 stale pin from F-07, since fixed (58/0), and the 2 known diagnostics. The skips match the baseline, plus Q-23 (F-04). |
| `R CMD check` ([check.log](check.log)) | 3 WARNINGs, 3 NOTEs, the same sections as the baseline. The code-problems NOTE lost F-09's "no visible binding for 'parameters'". |
| Skip baseline ([test-baseline.tsv](test-baseline.tsv), [test-baseline.md](test-baseline.md)) | Recorded in both load modes. The rule for KIT-102 onwards is in `quarantined-tests.md`. |

## The three proteomics journeys and protocol v1.7

**`101m` (DIA_LIMPA):**
- **The rejection:** an oracle mismatch on `design/components/peptide_data`. The `limpa65_screen` record still held its v1.6 fingerprint. v1.7 (KIT-092) re-hashes a record that nests tables, and protocol 9.3 requires a re-fingerprint before any new comparison.
- **Re-fingerprint:** done from its store, with KIT-092's method ([refingerprint.R](refingerprint.R), [refingerprint.txt](refingerprint.txt)). The v1.6 recomputation equals the record, one `result_table` component was re-hashed, and 7 products were renamed.
- **Result:** re-verified offline as accepted and exact ([reverify.txt](reverify.txt)).

**`101l` and `101t` (LFQ and TMT fixtures):**
- **The rejection:** the only mismatch was the identity product `normalization_science`, which KIT-092 also listed as changed by v1.7.
- **No re-fingerprint possible:** the kit-065 stores no longer exist, so those oracles are listed as unavailable.
- **Replacement oracles:** in-process records of the same fixtures on `kit_101`, bound as `lfq101_fixture` and `tmt101_fixture`. Every state field of every step equals kit-065's record; only the `<type>:<sha12>` result products differ ([fixture-oracles-cross-build.txt](fixture-oracles-cross-build.txt)).
- **Result:** both re-verified offline against them as accepted and exact.
- **Claim limit:** for LFQ and TMT, the cross-build check rests on that state-field equality, not on a v1.7 re-fingerprint.

## Scientific changes, each approved by the user

- **F-18:** MaxQuant and FragPipe read a 0 intensity ("not quantified") as NA, and F-26 does the same for PD TMT.
  - The synthetic `lfq_300` oracles were recorded with zeros kept, and now carry `recorded_before`.
  - A realistic LFQ corpus (`lfq_300_mnar`: 8.01% missing, intensity-dependent, 1,508 proteins complete) has oracles recorded on `kit_101` (RUV picks 3%).
  - The KIT-054 real-data oracles are re-run at the real-data gate (P-04).
- **F-19:** the DIA-NN fixtures now carry a properly normalised `Precursor.Normalised`. On DIA_LIMPA this takes the lane from 23 false positives to 4.
- **F-23:** a normalisation chosen on DIA_LIMPA now reaches the matrix DA fits. The default `none` is unchanged and byte-identical.

## Queued and parked

- **Queued:** F-16 to KIT-107; F-25 and the two display-only F-24 sites to KIT-104.
- **New ticket:** the statistical caveats (q-value method, RUV correction mode, DIA_LIMPA normalisation point, RUV control selection) are KIT-114.
- **Parked:** P-01 to P-06 in [../gui-workover-parked-decisions.md](../gui-workover-parked-decisions.md).

## Process notes

- **Correction budget:** the ticket ran several correction rounds, each authorised by the user ("fix all findings", "fix all", "of course fix it", "Approve", and the standing authority of 2026-11-07).
- **Freeze base:** the entry suite's freeze base predates F-27, F-28 and the cross-omic pin. The final-code run shows they change no capture.
- **Rule breach:** one diagnosis agent used a read-only `pgrep -f` once in a wait loop, against the host rule. It stopped and waited on files instead.

## Rollback

`git revert` of the KIT-101 commits. Tests, fixtures and docs revert cleanly. The R fixes are independent, one per finding.
