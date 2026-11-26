# Writes script-diff.md: every changed line range of KIT-281's scripts against KIT-211's on staging (difflib, line by line), with the ticket
# item or the declared reason of each. Usage, from the package root: python3 -I docs/qa/workflow-execution/kit/kit-281/make-script-diff.py
import difflib, hashlib

K211, K281 = 'docs/qa/workflow-execution/kit/kit-211/', 'docs/qa/workflow-execution/kit/kit-281/'


def lines(p):
    return open(p, encoding='utf-8').read().split('\n')


def sha(p):
    return hashlib.sha256(open(p, 'rb').read()).hexdigest()


def ops(a, b):
    return difflib.SequenceMatcher(None, a, b, autojunk=False).get_opcodes()


def rng(x, y):
    return '%d-%d' % (x + 1, y) if y - x > 1 else ('%d' % (x + 1) if y - x == 1 else '(none; inserted after %d)' % x)


SIM = {  # old start line -> (class, what)
    1: ('ticket :1-22', 'header comments (ticket id, usage, environment, the G3 and bound-tree notes)'),
    5: ('ticket :1-22', 'header comments'), 11: ('ticket :1-22', 'header comments'), 14: ('ticket :1-22', 'header comments'),
    20: ('ticket :1-22', 'header comments'),
    281: ('comment', 'the comment above CELLS (Av and Ag keep their entries and the generator lines; KIT-281 runs neither)'),
    289: ('ticket :289', '`CELL_MASKS` = A to D at MASKS[1:3] (the twelve vd settings)'),
    293: ('ticket :296 (comment :293-295)', 'ARMS comment and `ARMS <- list(main = c("S", "B", "B1", "N"))`'),
    299: ('ticket :299', '`SUBSETS` gains "g_added"'),
    368: ('ticket :368-402', '`proteinTable`: the one new argument `kept_b`, LAST, named at both call sites (no positional shift)'),
    397: ('ticket :368-402', '`proteinTable`: the new column `in_b = h_ids %in% kept_b`'),
    426: ('ticket :418-429', '`armRows`: for arm B1 only, one more part `aggregateSubset(tab, !tab$in_b, "g_added", c("class", "stratum", "ell"))`'),
    460: ('ticket :459-467', '`kept$A` removed'), 462: ('ticket :459-467', '`kept$Bp` removed'),
    464: ('ticket :459-467', '`kept$B1 <- bFilterKept(masked, sim$group, MIN_REPS, 1)` unconditional (it was inside `if (ag)`)'),
    483: ('ticket :482-484', 'the `proteinTable(...)` calls pass `kept_b = kept$B` (the arms call, and the C0 call that cannot run)'),
    497: ('NOT ON THE TICKET\'S LIST (ticket section 7 requires it): ADAPTED', 'E3: KIT-211 tested arm A\'s DA input (no longer an arm: it would be a constant TRUE and the lifecycle summary would record a check that was not run); here `identical(input$B1, masked[kept$B1, , drop = FALSE]) && setequal(das$B1$id, kept$B1)`'),
    510: ('NOT ON THE TICKET\'S LIST: ADAPTED', 'the metrics column `e3_pass` carries E3 on arm B1 (was arm A)'),
    516: ('ticket :516-527', 'lifecycle: E2\'s 1-group variant unconditional; the `e9_before` line removed'),
    520: ('ticket :516-527', 'lifecycle: `package_b1`/`e2g1` unconditional; the `helper`/`e9` lines removed'),
    525: ('ticket :516-527', 'lifecycle: the `e9_pass`, `e9_warnings` columns removed'),
    571: ('new code (ticket "Added" (c); K281-P05; section 8; section 3 commit rule)', '`partialGuard()` (scoped partial-match guard), `G3_TABLE`, `COUNTED_FIRST`, `COUNTED_CAP`, `KIT280_MARKER`'),
    573: ('ticket :573', 'usage message names sim-k281.R'),
    588: ('new code (guard repair, coordinator\'s check)', 'the last-net `options(error = )` writer: records `instrument-stop.txt` (the scorer refuses it) before the non-zero exit'),
    591: ('new code (task item 2; ticket section 8)', 'the refusal of a counted job whose range lies outside rule G3\'s table (g3-rule.tsv)'),
    600: ('ticket :600; "Added" (a), (b)', 'console line names KIT-281; `progress.tsv` header and `provenance-sim-start.json` written before the first task'),
    604: ('new code (ticket section 3 commit rule)', 'the assertion on the replay marker table: `threshold_values_in_filter` must be absent'),
    609: ('ticket "Added" (b)', 'one line per completed task (cell, r, seconds) appended to progress.tsv'),
    625: ('ticket :625', 'console line names KIT-281'),
    629: ('ticket :629-647; "Added" (c)', 'the lifecycle gathering runs under `partialGuard`'),
    634: ('ticket :634-645', 'the lifecycle summary has four checks E2, E3, E7, E10 (E9 row dropped); its E2, E3 details and F1 text restated; `[[ ]]` reads'),
    645: ('ticket :634-645', '`e9-warnings.tsv` dropped'),
    647: ('ticket :634-645', 'end of the guarded lifecycle block'),
    649: ('ticket "Added" (c)', '`assemble(...)` runs under `partialGuard`'),
    660: ('ticket "Added" (c)', 'the determinism comparison (own table handling only; not the re-run of the tasks) runs under `partialGuard`'),
    670: ('ticket "Added" (c)', 'end of the guarded determinism comparison'),
    677: ('ticket :678; "Added" (a)', 'the final provenance says KIT-281 and the script path; adds status, smoke, counted, marker names, G3 table hash, the guard\'s scope; `seconds_s` read under `partialGuard`'),
    694: ('ticket :694', 'console line names KIT-281'),
}

out = ['# script-diff.md: the changed lines of KIT-281\'s scripts against KIT-211\'s (generated by make-script-diff.py)', '',
       'KIT-211\'s files are those committed on `v0.9-staging` (`docs/qa/workflow-execution/kit/kit-211/`; SHA-256 `51c2d467...` sim, `aefc66ba...` score,',
       '`a387ac1e...` l1, asserted by `e0-k281.py`). Line numbers: "old" = KIT-211\'s file, "new" = KIT-281\'s file. Every changed block of the simulation',
       'carries a `KIT-281` comment or is in the header. **Verbatim parts are located by their markers, not by line numbers** (`e0-k281.py`).', '']
for name, key in (('sim', 'sim-k'), ):
    a, b = lines(K211 + 'sim-k211.R'), lines(K281 + 'sim-k281.R')
    out += ['## `sim-k281.R` against `sim-k211.R` (complete)', '', 'SHA-256 of `sim-k281.R`: `%s`' % sha(K281 + 'sim-k281.R'), '',
            '| old lines | new lines | ticket item or status | what |', '|---|---|---|---|']
    for tag, i1, i2, j1, j2 in ops(a, b):
        if tag == 'equal':
            continue
        cls, what = SIM.get(i1 + 1, ('UNANNOTATED', '?'))
        if i1 + 1 == 20 and False:
            pass
        out.append('| %s | %s | %s | %s |' % (rng(i1, i2), rng(j1, j2), cls, what))
    out += ['', 'Everything not in the table is byte-identical to `sim-k211.R`, including: the verbatim block (KIT-211 :24-239, :241-279: KIT-197\'s lines plus',
            'KIT-211\'s two guarded lines), `OLD_RULE` (:291-292, Av and Ag entries harmless), `bFilterKept`, `packageBKept`, `lfcDeviation`,',
            '`aggregateSubset`, `runTask`\'s DA calls, `writeTask`, `assemble`, `readPart`, `partPaths`. `helperShareKept`, `shareKept`, `masked_a` and the',
            '`ag` branches remain as unused definitions (declared harmless). `e0-k281.py` asserts the region and these named functions.', '']
    ticket_items = ('Not in the ticket\'s list but made: the E3 adaptation (:497-498 and :510), the counted-job refusal (:591), the marker assertion (:604; the ticket\'s',
                    'commit rule, implemented as the replay-marker assertion), `provenance-sim-start.json` (the start record), and the `smoke` flag (`K281_SMOKE=1`; the',
                    'ticket requires the scorer to stop on a smoke flag, which needs a producer).')
    out += list(ticket_items) + ['']
a, b = lines(K211 + 'score-k211.R'), lines(K281 + 'score-k281.R')
out += ['## `score-k281.R` against `score-k211.R`', '', 'SHA-256 of `score-k281.R`: `%s`' % sha(K281 + 'score-k281.R'), '',
        'Identical ranges (3 lines or more; KIT-211 line to KIT-281 line). Everything else in KIT-281\'s file is adapted or new (list below).', '',
        '| old lines | new lines | content |', '|---|---|---|']
CONTENT = {42: 'constants STRATA, VD (unchanged lines among the adapted constants), `k211Abort`, the verbatim blocks (num/writeTsv, jackSe, quant)',
           71: 'readInput\'s first lines (adapted function; E0 does not assert them)', 78: 'readInput\'s accumulators',
           91: '`armAtoms` up to `n0t_kept_both` (the new atoms are inserted before `calls <- ...`) ', 169: '`armAtoms` end', 188: '`armStats` (new entries inserted at the head of `out`)',
           220: '`computeStats` head', 234: '`computeStats` (the Phi block)', 246: '`computeStats` end', 252: '`bootstrapWeights`, `resampleStats`, `summariseStats`, `q1`, `limitsOf` (E0: identical by name)',
           327: '`fmt`, `mdTable` (E0: identical by name)', 433: 'falsifier\'s interval arithmetic lines'}
for tag, i1, i2, j1, j2 in ops(a, b):
    if tag == 'equal' and i2 - i1 >= 3:
        out.append('| %s | %s | %s |' % (rng(i1, i2), rng(j1, j2), CONTENT.get(i1 + 1, '')))
out += ['', 'Adapted (ticket section 6 list) or new, by KIT-281 line:', '',
        '- constants (new 19-58): `CELLS_ALL`, `CELL_MASKS`, `ARMS_MAIN`, `CANDIDATES`, `POOL`, `SEED` 280000001, `ALPHA_PRIMARY` 0.025, `ALPHA_SETTING` 0.05/12,',
        '  `ALPHA_F_POOL` 0.05, `THIN`, `AG_STATS`, new `RHO_MARGIN`, `FDP_LIMIT`, `GUARD_MARGIN`, `POWER_MARGIN`, `PC_MARGIN`; `ARMS_AG`, `ALPHA_PC`, `ALPHA_F_SETTING`,',
        '  `MARGIN`, `MCAR`, `CELLS_DECIDING`\'s Av/Ag entries not carried; new `DESIGN_STATS`, `G3_TABLE`, `ROW_*`/`METRIC_*` column vocabularies; `options(warnPartialMatchDollar = TRUE)`.',
        '- input (85-216): `readInput` is ADAPTED (lesson 1: `colClasses = "character"`, `na.strings` empty, explicit numeric and logical conversion, vocabulary checks, the',
        '  declared range, the smoke flag, instrument-stop.txt) with new helpers `readText`, `toNum`, `toLogical`, `needColumns`, `needLabels`, `declaredRange`.',
        '- `armAtoms`: new atoms `n0_gadd`, `p05_gadd`, `n0_gadd_l0/l1`, `p05_gadd_l0/l1`, `bhV_whole`, `bhT_whole` (new 296-305).',
        '- `armStats`: new `alpha_gadd`, `alpha_gadd_l0`, `alpha_gadd_l1`, `alpha_whole`, `fdr_bh`, and `gadd_per_r` (a count per replicate-mask, reported).',
        '- `computeStats` (new 360-388): :232-233 and :240-245 removed; new `rhog.B1`, `rhog_l0.B1`, `rhog_l1.B1`, `dfdrbh.B1`, `dpowbh.B1`, `rhoratio.B1`; `Phi.S`, `Phi.B`, `Phi.B1`',
        '  over `N.rho_l` by the existing Phi loop (the `bandratio` line is removed with C0).',
        '- NEW decision section (the scoring, 451 to the end) replacing KIT-211 :336-593: checks table, the pools, the deciding row and guards G1a, G1b, G2, F5\', the positive',
        '  control, per-design, per-setting and by-flag tables, the outcome rule, `decision-table.md`, `provenance-score.json`, `progress-score.tsv`, the driver with the',
        '  instrument-stop handling (`stopRecord`). Not carried: `primaryState`, `costState`, `classOf`, `outcomeOf`, the A, B\', C0, Av and Ag paths, the 24-setting pool,',
        '  the mcar pool, the per-setting falsifier intervals.',
        '- After the script audit (S02, S03, S08): a non-finite ESTIMATE of rho_g or a guard is an instrument stop (`not_formed`); the q-method count is reported, not a stop',
        '  (`fallback_by`, `qmethod-fallback.tsv`, `n_qmethod_fallback` and `notes` in `outcome.tsv`); every `addCheck`, `stopRecord` and `progress` call names its arguments.', '']
a, b = lines(K211 + 'l1-k211.R'), lines(K281 + 'l1-k281.R')
out += ['## `l1-k281.R` against `l1-k211.R`', '', 'SHA-256 of `l1-k281.R`: `%s`' % sha(K281 + 'l1-k281.R'), '',
        'Verbatim (E0, markers): `run` (l1-k197.R:34-37) and `sha` (l1-k197.R:38-38). Adapted: the header; the partial-match guard (global); provenance and progress before the',
        'first unit; G0 as a check; E2, E3, E7, E10 (E9 dropped); E8; **E1\'** in place of E1 (S and B against KIT-211\'s committed rows); **E4** is the scorer self-test (`selftest/`)',
        'in place of KIT-211\'s e4 tables; E5 (scorer twice, mechanics env); E6 (three masks); G2 (75 minutes); **new G3** (the count, rule G3 as `g3Rule`, `g3-count.tsv`,',
        '`g3-rule.tsv`); after the script audit: E11 (counted-range refusal through the last net, RC1-5), E1\'\'s missing-unit comparison (S04), G2 from set-up plus script with derived timeouts (S05, `g2-plan.tsv`), named arguments in every `addCheck`, `run` and `failureNote` call (S08), PATH for the self-test runner. Unchanged ranges (3 lines or more):', '', '| old lines | new lines |', '|---|---|']
for tag, i1, i2, j1, j2 in ops(a, b):
    if tag == 'equal' and i2 - i1 >= 3:
        out.append('| %s | %s |' % (rng(i1, i2), rng(j1, j2)))
open(K281 + 'script-diff.md', 'w').write('\n'.join(out) + '\n')
print('script-diff.md written')
