### Tests (head): 466 files; passed 47855, failed 928, errors 13, skipped 582; ok 431, failed 35, timeout 0, crashed 0

### Tests (base): 466 files; passed 47206, failed 924, errors 13, skipped 578; ok 430, failed 35, timeout 0, crashed 0, absent 1

### Base versus head: same 455, new-failure 1, message-changed 3, count-changed 6, fixed 0, only-on-head 1, compare-error 0
**Not clean under the rule**: see the rows below and compare-failures.txt

| file | base p/f/e/s | head p/f/e/s | verdict | tests | note |
|---|---|---|---|---|---|
| test-gui-entry-proteomics.R | 189/4/0/0 failed | 185/8/0/0 failed | new-failure | new-failure: KIT-101 entry lane prot_dia_limpa: DIA-NN, the DIA_LIMPA route; counts: KIT-101 entry lane prot_dia_limpa: DIA-NN, the DIA_LIMPA route (passed 38 -> 34, skipped 0 -> 0) |  |
| test-gui-entry-transcriptomics.R | 16/13/0/0 failed | 16/13/0/0 failed | message-changed | message-changed: KIT-101 entry lane trx_screen: the KIT-040 count corpus, 12 samples, 300 genes |  |
| test-integration-enrichment-mofa.R | 0/0/0/2 ok | 0/0/0/3 ok | count-changed | counts: KIT-184 K13: with MOFA2 trained on corrected values, camera with the removed factors, recorded and repeated offline (passed NA -> 0, skipped NA -> 1) |  |
| test-integration-enrichment-tests.R | 180/0/0/2 ok | 181/0/0/2 ok | count-changed | counts: KIT-183 T11: the tab: the test selector, its default, notes, counts, caution, seed label and key (passed 25 -> 26, skipped 0 -> 0) |  |
| test-integration-gui.R | 395/0/0/0 ok | 397/0/0/0 ok | count-changed | counts: KIT-138: the correction for unwanted variation each input names is stated in words (passed 6 -> 8, skipped 0 -> 0) |  |
| test-integration-input-lipid-metab.R | 132/0/0/0 ok | 208/0/0/0 ok | count-changed | counts: KIT-135: a new RUV design writes the RUV state's corrected table and names the design (M10); design none is normalised (passed 21 -> 40, skipped 0 -> 0); counts: KIT-135: lipidomics with today's RUV writes the corrected matrix per assay as omic_assay views, beside the RDS (passed 47 -> 65, skipped 0 -> 0); counts: KIT-135: metabolomics with today's RUV writes the corrected matrix per assay |  |
| test-integration-input-transcriptomics.R | 57/0/0/0 ok | 78/0/0/0 ok | count-changed | counts: KIT-136: a new design names its RUV record fields; an export from an earlier RUV step is refused until run again (passed 10 -> 19, skipped 0 -> 0); counts: KIT-136: with RUV applied the input is the accepted RUVIII export (INT-1), corrected, beside today's TSV (passed 28 -> 38, skipped 0 -> 0); counts: KIT-136: with RUV off (design none) the TMM logCPM is written, normalised, without a RUV |  |
| test-integration-removed-factors.R | 0/0/0/0 absent | 493/0/0/3 ok | only-on-head |  |  |
| test-omic-enrich-module.R | 1051/19/0/0 failed | 1051/19/0/0 failed | same |  | identical message without the actual digest value: (KIT-113 W1) each W1 (b) session's committed record, results, product, provenance and analysis record step |
| test-omic-summary-module.R | 1675/28/0/0 failed | 1675/28/0/0 failed | message-changed | message-changed: (b) trx_main: the outputs, notifications, progress, reactive values, calls, global and files of each step are today's; message-changed: (b) trx_managed: the outputs, notifications, progress, reactive values, calls, global and files of each step are today's |  |
| test-omics-summary-core.R | 676/25/0/2 failed | 676/25/0/2 failed | message-changed | message-changed: (W3, core side) trx_main: each entry gives W1 (b)'s notifications, progress, arguments, global and files; message-changed: (W3, core side) trx_managed: each entry gives W1 (b)'s notifications, progress, arguments, global and files |  |
| test-workflow-integration-input.R | 181/0/0/0 ok | 241/0/0/0 ok | count-changed | counts: KIT-061 (D14): DIA_LIMPA in covariates mode writes the RUV-adjusted product, not the state's uncorrected table (passed 17 -> 27, skipped 0 -> 0); counts: KIT-061 (K061-1): after a classic re-run with RUV skipped the head's normalised matrix is written, not the old correction (passed 13 -> 15, skipped 0 -> 0); counts: KIT-061: a classic route with RUV writes the corrected matrix after the c |  |

### R CMD check, base versus head
```
base: ERRORs 0, WARNINGs 3, NOTEs 3 (check exit 0)
head: ERRORs 0, WARNINGs 3, NOTEs 3 (check exit 0)
findings identical on both: 6; only on head: 0; only on base: 0

== only on head (new)
none

== only on base (gone)
none
```

**Run outcome: success** (every head file ran to completion; failing expectations, if any, are in the tables)
