# DIA fixture and measurement protocol work

DIA-002, version 1.0, fixture/protocol preparation complete. Runtime package
source is unchanged. No G0 or scientific/performance acceptance is implied.

## Freeze predecessors

These exact version-0.1 texts and hashes are preserved before recording approval
and fixture dispositions. Subsequent metadata changes do not rewrite preflight
receipts, regenerate inputs or upgrade the failed capacity case.

```json
[
  {
    "path": "tests/testdata/workflow-checkpoints/scenarios.json",
    "sha256": "54f64237b6d4ae554bfa26b199a5988844d507279ea2a60f2fb8029b10d539ea",
    "content": "{\n  \"schema\": \"multischolar.workflow_checkpoint_recipes.v1\",\n  \"version\": \"0.1-preflight\",\n  \"baseline\": \"c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017\",\n  \"routes\": {\n    \"iq\": \"first\",\n    \"limpa\": \"recipes only until DIA-017 real execution; identity doubles not science\"\n  },\n  \"small_default\": \"small\",\n  \"max_chunk_rows\": 200000,\n  \"recipes\": [\n    {\n      \"id\": \"small\",\n      \"n_proteins\": 64,\n      \"n_samples\": 12,\n      \"peptides_per_protein\": 3,\n      \"forms_per_peptide\": 2,\n      \"history_depth\": 6,\n      \"annotation_chars\": 48,\n      \"seed\": 902001,\n      \"purpose\": \"ordinary correctness; not representative performance\",\n      \"precursor_rows\": 4608,\n      \"protein_numeric_bytes\": 6144,\n      \"preflight\": \"pending\",\n      \"frozen_for_timing\": false\n    },\n    {\n      \"id\": \"representative\",\n      \"n_proteins\": 10000,\n      \"n_samples\": 48,\n      \"peptides_per_protein\": 4,\n      \"forms_per_peptide\": 2,\n      \"history_depth\": 14,\n      \"annotation_chars\": 96,\n      \"seed\": 902002,\n      \"purpose\": \"plausible deep DIA precursor workload; shape proxy, not empirical representativeness\",\n      \"precursor_rows\": 3840000,\n      \"protein_numeric_bytes\": 3840000,\n      \"preflight\": \"pending\",\n      \"frozen_for_timing\": false\n    },\n    {\n      \"id\": \"tall\",\n      \"n_proteins\": 12000,\n      \"n_samples\": 24,\n      \"peptides_per_protein\": 6,\n      \"forms_per_peptide\": 2,\n      \"history_depth\": 14,\n      \"annotation_chars\": 96,\n      \"seed\": 902003,\n      \"purpose\": \"tall precursor stress\",\n      \"precursor_rows\": 3456000,\n      \"protein_numeric_bytes\": 2304000,\n      \"preflight\": \"pending\",\n      \"frozen_for_timing\": false\n    },\n    {\n      \"id\": \"wide\",\n      \"n_proteins\": 8000,\n      \"n_samples\": 256,\n      \"peptides_per_protein\": 2,\n      \"forms_per_peptide\": 2,\n      \"history_depth\": 14,\n      \"annotation_chars\": 96,\n      \"seed\": 902004,\n      \"purpose\": \"many-sample correlation/PCA and lookup stress\",\n      \"precursor_rows\": 8192000,\n      \"protein_numeric_bytes\": 16384000,\n      \"preflight\": \"pending\",\n      \"frozen_for_timing\": false\n    },\n    {\n      \"id\": \"annotation_heavy\",\n      \"n_proteins\": 10000,\n      \"n_samples\": 48,\n      \"peptides_per_protein\": 4,\n      \"forms_per_peptide\": 2,\n      \"history_depth\": 14,\n      \"annotation_chars\": 2048,\n      \"seed\": 902005,\n      \"purpose\": \"high-entropy annotation/component-sharing stress\",\n      \"precursor_rows\": 3840000,\n      \"protein_numeric_bytes\": 3840000,\n      \"preflight\": \"pending\",\n      \"frozen_for_timing\": false\n    },\n    {\n      \"id\": \"cycles\",\n      \"n_proteins\": 10000,\n      \"n_samples\": 48,\n      \"peptides_per_protein\": 4,\n      \"forms_per_peptide\": 2,\n      \"history_depth\": 14,\n      \"annotation_chars\": 96,\n      \"seed\": 902006,\n      \"cycles\": 10,\n      \"purpose\": \"ten repeated normalize/revert/style/export cycles, not independent replicates\",\n      \"precursor_rows\": 3840000,\n      \"protein_numeric_bytes\": 3840000,\n      \"preflight\": \"pending\",\n      \"frozen_for_timing\": false\n    },\n    {\n      \"id\": \"capacity\",\n      \"n_proteins\": 20000,\n      \"n_samples\": 512,\n      \"peptides_per_protein\": 4,\n      \"forms_per_peptide\": 2,\n      \"history_depth\": 24,\n      \"annotation_chars\": 256,\n      \"seed\": 902007,\n      \"purpose\": \"capacity-only stress; not an empirical proteome/population or finite-speedup control\",\n      \"precursor_rows\": 81920000,\n      \"protein_numeric_bytes\": 81920000,\n      \"preflight\": \"pending\",\n      \"frozen_for_timing\": false\n    }\n  ],\n  \"scenario_matrix\": [\n    {\n      \"id\": \"C01\",\n      \"scenario\": \"Save, close/reopen, hydrate each supported S4 class\",\n      \"expected\": \"Exact object and metadata round-trip; usable S4 dispatch\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C02\",\n      \"scenario\": \"Find by method/design/assay/source factors\",\n      \"expected\": \"Correct match/no-match/ambiguity; zero payload reads\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C03\",\n      \"scenario\": \"Small versus large objects; zero optional cache\",\n      \"expected\": \"Correct uncached behavior; no silent cap increase\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C04\",\n      \"scenario\": \"Unchanged/shared metadata versus changed quantities\",\n      \"expected\": \"Memory accounting distinguishes shared allocations; changed science never hits an old product\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C05\",\n      \"scenario\": \"Row/column subsets, reordering, duplicate selectors and patches\",\n      \"expected\": \"Defined selector semantics, exact IDs/order and safe cache invalidation\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C06\",\n      \"scenario\": \"Factors/ordered factors, names, row names, annotations and nested args\",\n      \"expected\": \"Preservation under the codec; explicit unsupported-value errors\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C07\",\n      \"scenario\": \"NA, NaN, infinities, zero, all-missing and constant features\",\n      \"expected\": \"Preserve values and existing unavailable-panel/rejection behavior\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C08\",\n      \"scenario\": \"Forward history, immediate revert, deep revert, reset and reapply\",\n      \"expected\": \"Same state/history semantics; only requested state hydrated\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C09\",\n      \"scenario\": \"Revert then create a new branch of history\",\n      \"expected\": \"Old required history remains recoverable; no stale active refs\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C10\",\n      \"scenario\": \"Cache hit, eviction, reload and changed style\",\n      \"expected\": \"No refit on valid hits; correct invalidation; byte ceiling respected\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C11\",\n      \"scenario\": \"Original owners released versus deliberately retained aliases\",\n      \"expected\": \"Memory savings occur only when all relevant references are released\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C12\",\n      \"scenario\": \"Many checkpoints and metadata matches\",\n      \"expected\": \"Bounded metadata lookup and resident cache; no full-history hydrate\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C13\",\n      \"scenario\": \"Failed save, corruption, truncation, missing payload or invalid codec\",\n      \"expected\": \"Classed failure; previous active state/history remains usable\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C14\",\n      \"scenario\": \"Process exit during save or before acknowledgement\",\n      \"expected\": \"Persistent recovery and idempotent retry in a fresh process\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C15\",\n      \"scenario\": \"Cancellation and project replacement\",\n      \"expected\": \"Stale results cannot activate/publish into the replacement project\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C16\",\n      \"scenario\": \"Two stores/projects and conflicting writers\",\n      \"expected\": \"No cross-project reuse; reject or serialize conflicts explicitly\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C17\",\n      \"scenario\": \"Permission failure and disk-full injection\",\n      \"expected\": \"No falsely accepted checkpoint, lost previous output or orphaned authority\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C18\",\n      \"scenario\": \"Complete legacy export and timestamped/latest outputs\",\n      \"expected\": \"Full historical contents and consumer compatibility; export peak included\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C19\",\n      \"scenario\": \"Ten normalize/revert/style/export cycles in one process\",\n      \"expected\": \"Resources/cache bounds checked each cycle; retain every RSS-guard failure\",\n      \"recipes\": [\n        \"cycles\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C20\",\n      \"scenario\": \"Restart with no resident cache\",\n      \"expected\": \"Cold restore cost and compatibility measured independently of cache hits\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C21\",\n      \"scenario\": \"Same DA request, fit/result reuse and changed contrast\",\n      \"expected\": \"No fitting on a fully compatible result hit; correct invalidation/recomputation at each statistical layer\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C22\",\n      \"scenario\": \"Changed design, weights, sample order or tested feature population\",\n      \"expected\": \"Reject incompatible fits/results even when dimensions or filenames match\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C23\",\n      \"scenario\": \"Plot style-only change, scientific change and asset hit\",\n      \"expected\": \"Separate zero-fit, zero-source-hydration and zero-render expectations according to the actual hit level\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C24\",\n      \"scenario\": \"DA result pages, volcano projections and heatmap selections\",\n      \"expected\": \"Correct ordered values and full-population statistics; actual bounded reads for capable codecs\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C25\",\n      \"scenario\": \"One-off access versus frequent reuse; hit/miss mixtures\",\n      \"expected\": \"Include initial production/persistence, lookup, cold read, invalidation and eviction costs\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    },\n    {\n      \"id\": \"C26\",\n      \"scenario\": \"Combined consumer sequence\",\n      \"expected\": \"Saving for one consumer cannot shift more cost or peak memory into another consumer without reporting it\",\n      \"recipes\": [\n        \"small\",\n        \"representative\",\n        \"tall\",\n        \"wide\",\n        \"annotation_heavy\",\n        \"capacity\"\n      ],\n      \"execution\": \"not_run\",\n      \"fixture_variants\": [\n        \"independently_changed\",\n        \"shared_metadata_only\",\n        \"retained_alias\"\n      ],\n      \"owner\": \"DIA-002\"\n    }\n  ]\n}\n"
  },
  {
    "path": "tests/testdata/workflow-checkpoints/protocol.json",
    "sha256": "48a9ec2434744673e9e1f6007c653f429f89826f8072ce864df14913d75cda5a",
    "content": "{\n  \"schema\": \"multischolar.dia_measurement_protocol.v1\",\n  \"version\": \"0.1-proposed\",\n  \"source_control\": \"c7c12851ea4d9e96f91df0e1f9e7b91c2eb51017\",\n  \"status\": \"proposed_not_frozen\",\n  \"approval\": {\n    \"status\": \"requested\",\n    \"approved_by\": null,\n    \"question\": \"G2 30% lower retained heap; <=10% elapsed/CPU overhead; 30s max save or cold revert. G4 unchanged.\"\n  },\n  \"independence_unit\": \"fresh installed-package journey process tree and project for each independent replicate\",\n  \"pairing_block\": \"same recipe seed/settings/dependencies/host/boundary; original control and treatment each fresh, serial\",\n  \"repetitions\": 7,\n  \"member_order\": \"alternate control/candidate then candidate/control by block index; no prewarming outside boundary\",\n  \"paired_summary\": \"geometric mean of within-block candidate/control ratios (mean log ratios); no pooling shapes or routes\",\n  \"uncertainty_method\": {\n    \"method\": \"paired log-ratio Student-t one-sided upper confidence bound\",\n    \"confidence\": 0.95,\n    \"df\": 6,\n    \"assumptions\": \"independent blocks; approximately normal block log ratios; retain raw pairs; obvious violations hold instead of switching estimator\"\n  },\n  \"decision_rule\": \"all required point and upper-bound ratios must meet their thresholds, with science/files/limits valid; point miss fails, uncertainty-only miss holds\",\n  \"borderline_rule\": \"point threshold met but upper confidence bound exceeds threshold -> held; no extra favorable reruns\",\n  \"non_regression_rule\": \"candidate/matched immediate accepted comparator <=1.00 for both point and upper-bound elapsed and combined peak; no margin\",\n  \"combined_peak_definition\": {\n    \"metric\": \"fresh service cgroup v2 memory.peak\",\n    \"units\": \"bytes\",\n    \"membership\": \"GUI process, reusable worker if adopted, driver and descendants under one service\",\n    \"interval\": \"process startup through all required outputs/acknowledgement and teardown\",\n    \"shared_memory\": \"charged once by cgroup accounting; do not add RSS\",\n    \"file_cache\": \"include cache and kernel charges reported by same cgroup; separate host page-cache condition recorded\",\n    \"driver\": \"included identically for paired journeys\",\n    \"exclusions\": \"fixture generation and independent output oracle excluded from timed service; their own preflight accounted separately\"\n  },\n  \"g2\": {\n    \"retained_heap_ratio_max\": 0.7,\n    \"sequence_elapsed_ratio_max\": 1.1,\n    \"sequence_cpu_ratio_max\": 1.1,\n    \"save_max_seconds\": 30,\n    \"cold_revert_max_seconds\": 30,\n    \"retained_heap_measurement\": \"unique live R allocations across all owners after standardized GC, excluding freed objects; never sum object.size across aliases\",\n    \"scope\": \"qualified history workloads; capacity-only failures excluded from finite-ratio acceptance\"\n  },\n  \"g4\": {\n    \"elapsed_ratio_max\": 0.9,\n    \"combined_peak_ratio_max\": 0.9,\n    \"immediate_comparator_ratio_max\": 1,\n    \"required\": \"complete independent IQ then real Limpa journey; startup, full-history legacy RDS, copies, automatic outputs and included restore boundary\"\n  },\n  \"service_limits\": {\n    \"memory_max_bytes\": 8589934592,\n    \"memory_swap_max_bytes\": 0,\n    \"runtime_max_seconds\": 900,\n    \"optional_cache_max_bytes\": 67108864\n  },\n  \"run_budget\": {\n    \"qualified_blocks_per_build_recipe_route\": 7,\n    \"diagnostic_profiles_per_case\": 1,\n    \"fixture_preflights_per_recipe\": 1\n  },\n  \"retry_rule\": \"one replacement only for independently documented infrastructure-invalid block, preserving failed receipt; no retry of valid unfavorable outcomes\",\n  \"correction_budget\": {\n    \"design_corrections_per_failed_pilot\": 1,\n    \"requires\": \"measured attribution, scoped owner, new build ID; rerun only invalidated evidence\"\n  },\n  \"invalid_baseline_rule\": \"OOM/timeout/missing science or outputs cannot be a finite denominator; capacity-only observations labeled separately\",\n  \"missing_metric_rule\": \"required absent/nonfinite metrics -> held; never zero, impute or omit\",\n  \"boundaries\": {\n    \"baseline\": \"DIA-003 binds original browser actions and required-file inventory before timing\",\n    \"histories\": \"include initial save and realistic forward/revert/reset/reapply; include full export separately and in composed sequence\",\n    \"consumers\": \"resident vs full RDS vs small saved product vs selective read; creation cost and misses included\",\n    \"cold_state\": \"fresh process, empty optional cache, disjoint project; filesystem cache condition explicitly paired, not asserted cold from restart alone\"\n  },\n  \"routes\": {\n    \"iq\": {\n      \"status\": \"awaiting fixture preflight, approval and G0\"\n    },\n    \"limpa\": {\n      \"status\": \"held until DIA-017\",\n      \"inherits_iq_evidence\": false\n    }\n  }\n}\n"
  }
]
```


## Current boundaries

The five ticket-owned files define test helpers, tests, recipes, a proposed
measurement protocol and this evidence record. Builders use real package S4
classes. Protein allocation fixtures are explicitly synthetic, not claimed IQ
results. EList controls test identity only; real Limpa science remains DIA-017.
The helper allocates nothing at load time. Default tests create small data only.

Precursor intensities have feature effects, a treatment shift on every tenth
protein, independent continuous noise, MAR and intensity-dependent missingness.
Raw and precursor-normalized linear quantities are separate columns. Peptide
sequences contain only amino-acid letters; modified forms and charge states have
distinct precursor IDs. Per-protein seeds make chunk partitioning reproducible.
Every recipe has balanced biological groups and paired technical replicates.

The checked-in recipe dimensions and seeds are now frozen, with the recorded
preflight dispositions and user-approved intermediate allowance. Tall/wide/
annotation/capacity stress cases do not establish empirical
representativeness; later scaling claims need a real representative workload.
The capacity case is not a finite speedup denominator if its control fails.

## Protocol approval

User approved the G2 allowance on 2026-09-12: at least 30% lower retained unique R
heap, at most 10% total elapsed/CPU overhead for the history sequence, and maximum
30-second save/cold-revert latency for qualified cases. This is an intermediate
tradeoff only. Final G4 remains at least 10% lower complete elapsed and combined
peak memory, and no accepted-comparator regression. The user's exact response
was “Yes that allowance is fine”. This approves the intermediate history tradeoff,
not a slower GUI, a changed scientific method or an increased resource limit.
Protocol version 1.0 records that scope. G0 still requires an installed, qualified
control; this preparation ticket does not itself authorize candidate promotion.

The predeclared measurement procedure uses seven independent paired blocks per
recipe/build/route. Aggregate
within-block log ratios, with a predeclared one-sided 95% Student-t upper bound;
point failures fail, uncertainty-only failures hold. No pooling across routes or
shapes, no optional favorable reruns, and no treating ten within-session cycles
as independent replicates. DIA-003 owns the evaluator and adversarial receipt tests.

## Resource preflight

Dedicated external root approved/created:
/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5

The repository is on btrfs; /tmp is tmpfs and cannot support disk-offload claims.
Each preflight must run alone in a fresh systemd service with MemoryMax=8G,
MemorySwapMax=0 and RuntimeMaxSec=900. A no-op service probe verified those
settings could execute: success, zero swap. The host has existing swap usage;
this does not mean the fixture service may swap. No host-wide settings changed.

The materializer writes one report in chunks of at most 200000 rows, then opens
one complete peptide object, releases it, and builds one protein object. It never
retains an oracle or all history generations. It refuses to reuse an existing
case directory. Partial reports and receipts remain after failure. Current
preflight receipts measure feasibility only, not candidate/control performance.
Large outputs stay external; do not commit or automatically regenerate them in
ordinary tests. Retain this one preflight set for the upcoming control/harness;
no automatic repeated run or deletion of failed evidence is permitted.

## Evidence so far

- Initial exact fixture test command: 66 passed, one error. Actual precursor
  rollup dropped Q.Value, but the existing S4 validity still requires that field.
- After explicit baseline-error characterization: 74 passed, no warnings/skips.
  This did not repair the intermediate S4 or qualify it for valid checkpoint or
  scientific acceptance. Raw peptide and synthetic protein fixtures pass actual
  validity; rolled/imputed original intermediates remain a G0 correctness risk.
- Added actual IQ consumption: 81 assertions passed, no warnings/skips. IQ read
  2304 peptide quantities for 64 proteins and 12 samples and produced a valid
  real ProteinQuantitativeData result. This is not complete GUI/scientific parity.
- Preliminary representative preflight v01 passed: 3840000 precursor rows;
  556028391-byte external report; peptide object.size 445383144 bytes and protein
  object.size 8506184 bytes. Service output reported 40.398 seconds, 41.068 CPU
  seconds, approximately 2.4 GiB peak and zero swap. These are feasibility costs,
  not performance gains; object.size is not unique retained heap or combined peak.
- The completed transient service was collected, so a later systemctl query
  returned no live counters. The fixture materializer now records actual cgroup
  limits/peak and source/recipe hashes inside the process, and writes a running
  receipt before large work so OOM/worker loss cannot erase the attempt. The first
  receipt remains preliminary evidence and is not promoted to a frozen gate.
- Tall preflight v01 also passed: 3456000 precursor rows, 66.118 seconds service
  runtime and zero swap. It binds helper SHA-256
  29ba69bc649f5f273e1bf02a9ac3ccd7c3eb6cddae7b21a8df42ef8ad04ab47e and recipe
  SHA-256 54f64237b6d4ae554bfa26b199a5988844d507279ea2a60f2fb8029b10d539ea.
  Actual in-process limits were memory.max=8589934592 and memory.swap.max=0.
  The full numeric peak/footprint receipt is retained at tall/preflight.json in
  the external root. This is a resource preflight, not science/performance acceptance.
- Commands: Rscript --vanilla -e 'devtools::test(filter = "workflow-checkpoint-fixtures")'.
- The first resource service probe hit the sandbox user-bus restriction; the
  escalated retry succeeded. No limits were weakened.
- air is not installed. New helpers use four-space indentation and explicit
  namespaces; no repository-wide formatting or dependency installation occurred.

## Remaining downstream qualification

All six large recipes have been attempted; capacity failed with OOM and is not
an eligible finite control. The representative receipt remains explicitly
preliminary. Allocation/selector/patch controls, scenario construction details,
focused regressions and approval are recorded here. DIA-003 must bind the fresh
installed build and complete browser/output boundary; DIA-004 must qualify IQ.
The known original rollup validity defect and other DIA-001 baseline risks are
not fixed or accepted by this closeout. Real Limpa remains separately held.

The freeze changed metadata/dispositions only: a verifier compared every original
recipe input and scenario, G2/G4 number and service limit against the exact saved
version-0.1 documents above. The post-freeze fixture suite passed 99 assertions,
zero failures/warnings, with the documented service-only skip. No large case was
rerun, no failed receipt was replaced, and no performance gain is claimed.

Frozen manifest SHA-256: cfdf7da2fc5a6123da269f85011550d8ffd2e8a611de071e0cd1a52a6f35d152.
Frozen protocol SHA-256: 517edccb0906165698f3c1f165dfe7d961a17910d93363a60e843fa466c62379.

## Coverage and service-path verification, 2026-09-12

The first covr::file_coverage run measured 68.30189% for the new helper. Most
uncovered lines belonged to its intentionally guarded physical preflight path.
An explicitly enabled small service test then exercised real disk writes,
cgroup reads, duplicate-output refusal and a deliberately injected construction
failure with a retained failed receipt. Coverage increased to 96.60377%, with
112 successful assertions in that service run. These are fixture-helper coverage
results, not package-wide coverage or scientific/performance acceptance.

The small coverage service used the same 8 GiB/zero-swap/900-second constraints
and completed successfully in 19.576 seconds. Its reports/receipts remain under
coverage-small-v01 in the approved external root. Ordinary tests skip this
explicit-service test when DIA_FIXTURE_PREFLIGHT_ROOT is unset, so default tests
do not create external directories or rerun large recipes. Both coverage results
are retained; the earlier deficiency is not erased.

Command shape: pkgload::load_all(...); library(testthat);
covr::file_coverage("tests/testthat/helper-workflow-checkpoint-fixtures.R",
"tests/testthat/test-workflow-checkpoint-fixtures.R"). The successful service set
DIA_FIXTURE_PREFLIGHT_ROOT to the fresh coverage-small-v01 directory before
running that exact test file; no estimator was mocked. The injected failure
replaces only precursor construction to test failure-receipt retention.

## Scenario fixture construction contract

These are fixture inputs/oracles for the existing C01-C26 owners, not claims that
the future saved-entity operations already exist. The scenario manifest selects
the dimensions; these instructions specify which concrete variants to construct.
The raw peptide and synthetic protein constructors return real valid S4 classes.
Original rolled/imputed peptide intermediates are a separate known-invalid
baseline characterization, not substitutes for valid codec fixtures.

| Cases | Concrete fixture setup and expected distinction |
| --- | --- |
| C01 | Raw peptide, actual small IQ protein output and synthetic protein; preserve complete S4 slots/config/design exactly through RDS; retain invalid rolled-state rejection separately |
| C02 | Same quantity generation with normalization method none/cyclicloess/quantile/scale, changed design-factor levels, typed integer/double/text/missing factors; include two otherwise identical metadata matches to test ambiguity |
| C03 | Small and every declared large shape with optional cache zero and 64 MiB; capacity protein payload exceeds cache size while still constructible independently; cache limits are not object-format limits |
| C04 | Six independently generated small protein states versus one payload shared by six metadata revisions; allocation addresses must differ/share respectively; large callers iterate generation indices without accumulating an oracle |
| C05 | Rows 3,1,3 and columns ID, fourth sample, second sample; independent row/sample reordering; change row 9/sample 2 to 123.456 in a copy; original quantities/annotations/design remain identical |
| C06 | Ordered acquisition factor including unused middle level; unused batch factor level; nested named integers/factors, automatic and subset-derived row names, table units attribute |
| C07 | NA, NaN, +Inf, -Inf, zero and signed zero in distinct cells, an all-missing row and a constant row; preserve classifications exactly; do not require estimators to accept rejected inputs |
| C08-C09 | Save each declared generation, immediate/deep revert, reapply existing name, then branch; preserve inactive state records and audit lineage rather than trimming stored contents to active history |
| C10 | Exact input/style hit, changed style only, changed quantity and changed design; distinguish resident eviction from deletion of persistent checkpoint/product dependencies |
| C11 | Keep an explicit external list of original/shared S4 aliases while dropping manager references; compare with a fresh process that releases those owners; no reclaimed-memory claim while aliases persist |
| C12 | Iterate full declared history_depth and metadata revisions; repeat typed metadata lookup without hydrating payloads (later utility test), not repeated serialization of one object disguised as independent history |
| C13 | Save a valid small object; independently corrupt/truncate/remove its owned payload or change codec tag in test-owned artifacts; expect rejection and unchanged prior active state |
| C14 | Use a subprocess and the same valid small object; terminate at declared before/after-write and before-acknowledgement hooks; reopen in another process, never infer recovery from the old R6 object |
| C15-C16 | Two independent project roots with identical recipe/seed/dimensions but distinct incarnations; replace one project while a request is outstanding and attempt conflicting writers |
| C17 | Explicit file-write refusal/disk-full injection on small test-owned outputs; keep committed science/previous complete outputs and classify failure; real capacity OOM remains separate retained evidence |
| C18 | Complete states list including inactive branches and both current S4 aliases, all timestamped/latest/TSV/figure obligations; full-history RDS cannot become reference-only |
| C19 | cycles recipe and exactly ten normalize/revert/style/export actions in one session; the input preflight does not execute or qualify these actions |
| C20 | Close all parent-process objects/cache; fresh-process reopen and external legacy restore from immutable files; changed/corrupt external source must still be validated |
| C21-C22 | Actual IQ output for DA; same full result request, changed contrast, changed factor levels/order, technical blocks, model-selected samples and feature population; Limpa same E with different standard.error or n.observations and aligned support-order changes |
| C23 | Distinct scientific-product, style-only and exact rendered-asset requests; preserve each effective PCA/density population and distinguish zero-fit from zero-render expectations |
| C24 | Full contrast statistics before table pages/volcano filters; duplicate/reordered heatmap row selectors and design-ordered sample selections; no page-specific FDR |
| C25-C26 | One cold production/persistence followed by zero, one and repeated consumers; explicit hit/miss/style/quantity-change sequences; sum exclusive phases only and retain full export/startup/publication costs |

Failure hooks, cache operations and consumer execution above belong to the
existing downstream tickets, not this fixture helper. Numeric method outputs
must be captured from the actual qualified control by the scientific/measurement
owners; synthetic allocation values are never labeled fitted IQ/Limpa results.


## Rollback procedure

No storage-backed session or runtime adapter exists. Preserve the original
resident path and all earlier failed receipts. Supersede only these fixture
recipes/helpers/protocol if a case is rejected; retain provenance and do not
reset the branch or alter scientific methods to make a fixture pass.

## Retained initial failing diagnostic

The later green characterization does not supersede this original validity failure.

```text
Rscript --vanilla -e 'devtools::test(filter = "workflow-checkpoint-fixtures")'
ℹ Testing MultiScholaR

-------------------------------------------------------------------
Welcome to MultiScholaR!

IMPORTANT: Please run MultiScholaR::loadDependencies() 
           to ensure all required packages are installed and loaded.
-------------------------------------------------------------------
✔ | F W  S  OK | Context

⠏ |          0 | workflow-checkpoint-fixtures                                   
⠸ |         14 | workflow-checkpoint-fixtures                                   
⠸ |         24 | workflow-checkpoint-fixtures                                   
⠴ |         56 | workflow-checkpoint-fixtures                                   
⠦ | 1       66 | workflow-checkpoint-fixtures                                   
✖ | 1       66 | workflow-checkpoint-fixtures [1.4s]
────────────────────────────────────────────────────────────────────────────────
Error ('test-workflow-checkpoint-fixtures.R:129:5'): fixtures exercise real precursor aggregation and peptide imputation
Error in `validityMethod(object)`: Q value column must be in the peptide data table
Backtrace:
    ▆
 1. ├─testthat::expect_true(methods::validObject(rolled)) at test-workflow-checkpoint-fixtures.R:129:5
 2. │ └─testthat::quasi_label(enquo(object), label)
 3. │   └─rlang::eval_bare(expr, quo_get_env(quo))
 4. └─methods::validObject(rolled)
 5.   ├─methods (local) anyStrings(validityMethod(object))
 6.   │ └─base::isTRUE(x)
 7.   └─MultiScholaR (local) validityMethod(object)
────────────────────────────────────────────────────────────────────────────────

══ Results ═════════════════════════════════════════════════════════════════════
Duration: 1.4 s

── Failed tests ────────────────────────────────────────────────────────────────
Error ('test-workflow-checkpoint-fixtures.R:129:5'): fixtures exercise real precursor aggregation and peptide imputation
Error in `validityMethod(object)`: Q value column must be in the peptide data table
Backtrace:
    ▆
 1. ├─testthat::expect_true(methods::validObject(rolled)) at test-workflow-checkpoint-fixtures.R:129:5
 2. │ └─testthat::quasi_label(enquo(object), label)
 3. │   └─rlang::eval_bare(expr, quo_get_env(quo))
 4. └─methods::validObject(rolled)
 5.   ├─methods (local) anyStrings(validityMethod(object))
 6.   │ └─base::isTRUE(x)
 7.   └─MultiScholaR (local) validityMethod(object)

[ FAIL 1 | WARN 0 | SKIP 0 | PASS 66 ]

```

## Preflight command provenance

The representative and tall runs used the following explicit commands. The first
has only its preliminary external receipt and service summary; the tall receipt
also binds the generator/recipe hashes and actual cgroup limits.

```bash
systemd-run --user --unit=multischolar-dia002-representative-v01 --wait --pipe --property=MemoryMax=8G --property=MemorySwapMax=0 --property=RuntimeMaxSec=900 --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 --setenv=VROOM_THREADS=1 --working-directory=/home/doktersmol/Projects/APAF/MultiScholaR /usr/bin/Rscript --vanilla -e 'pkgload::load_all(quiet=TRUE, helpers=FALSE, attach_testthat=FALSE); sys.source("tests/testthat/helper-workflow-checkpoint-fixtures.R", envir=environment()); preflightWorkflowCheckpointRecipe("representative", "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5", getwd())'

systemd-run --user --unit=multischolar-dia002-tall-v01 --wait --pipe --property=MemoryMax=8G --property=MemorySwapMax=0 --property=RuntimeMaxSec=900 --setenv=OMP_NUM_THREADS=1 --setenv=OPENBLAS_NUM_THREADS=1 --setenv=VROOM_THREADS=1 --working-directory=/home/doktersmol/Projects/APAF/MultiScholaR /usr/bin/Rscript --vanilla -e 'pkgload::load_all(quiet=TRUE, helpers=FALSE, attach_testthat=FALSE); sys.source("tests/testthat/helper-workflow-checkpoint-fixtures.R", envir=environment()); preflightWorkflowCheckpointRecipe("tall", "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5", getwd())'
```

Earlier focused rerun on the helper: 81 passed, 0 failed, 0 warnings,
0 skips. DIA-002 remains in progress; this is not an installed or performance gate.

## Continuation evidence, 2026-09-12

## Bound generated inputs

All six reports were fingerprinted once after their processes terminated. The
capacity file has 81920001 lines including its header, confirming generation
finished before the unsuccessful whole-object feasibility check. That does not
turn the OOM outcome into a pass.

```text
2414125be92a11fbe4342fdf354e0214915e30592bec66f3e9d6a984f9290841  /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/representative/report.tsv
41ace829401f18b9c700bd4b0f4223f8eb08b35f98a2691546906ebc221b9e3e  /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/tall/report.tsv
bc673ae8e691259e1b4bd68d69154969f15c7841c00387ff9f7ddff278f3664e  /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/wide/report.tsv
4a464eb9ffa37f8f335991eb5640a1fdfe060e806ec7f9818faa385feafd8b9e  /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/annotation_heavy/report.tsv
b9c19e6410d017d8072ff2383c67061280ffab3c47d38f7170d5ff2d8f7f4581  /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/cycles/report.tsv
93c5d6e63234ea755b69bf7546c3571a88fb62883754b2abd67c03c3968654a2  /home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/capacity/report.tsv
```

The exact-order IQ assertions also pass: final fixture suite 99 passed, zero
failures/warnings/skips. Hashing and line-count checks are separate input
verification operations, not timed candidate members.

Final focused regression command:
Rscript --vanilla -e 'devtools::test(filter = "workflow-checkpoint-fixtures|prot-(02s-identification-evidence|02b-qc-peptide-module-contracts)", reporter = "summary")'.
It completed without failures or warnings; the explicit service-only test skipped
in this ordinary process and was separately executed successfully in the guarded
coverage service. Current source/authority hash checks and git diff --check pass.
A first read-only verification command had a stray here-document marker after
its successful assertions and exited with NameError; the corrected command
passed. No source or evidence was modified by that failed command.


The interrupted wide process was recovered through its existing handle and found
completed successfully. It was not restarted. Annotation-heavy and cycle-input
preflights then ran as separate services under the same limits. The cycle input
preflight is not an execution of the ten-cycle acceptance scenario.

| Recipe | Rows | Report bytes | In-process peak bytes | Swap bytes |
| --- | ---: | ---: | ---: | ---: |
| annotation_heavy | 3840000 | 556026987 | 2772549632 | 0 |
| cycles | 3840000 | 556029924 | 2621620224 | 0 |
| representative | 3840000 | 556028391 | not captured (preliminary) | service summary: 0 |
| tall | 3456000 | 500607058 | 2562809856 | 0 |
| wide | 8192000 | 1185864748 | 4147425280 | 0 |

Small tests now pass 99 assertions with no failures, warnings or skips. New
checks use rlang allocation addresses to distinguish actual shared payloads
from independent generations, and verify copy-safe patches and typed factors.
No new runtime code or dependency was introduced. The capacity preflight is
still active; its initial running receipt is not a completed result.

Capacity monitor handoff: unit multischolar-dia002-capacity-v01.service,
unified exec session 36789, MainPID 1034408, started 2026-09-12 07:23:24 AEST.
At 07:28:20 AEST it was confirmed active/running, with MemorySwapCurrent=0,
MemoryMax=8589934592, MemorySwapMax=0 and approximately 5.05 GB current charge.
The explicit capacity recipe has 81920000 rows, not 40960000; its output was
still being written (6662033368 bytes at that observation). Await terminal
process evidence through that handle/unit. Do not restart from the running
receipt or an observation timeout. The service's 900-second deadline remains.

Completed receipt bindings (not comparative acceptance):

### Capacity terminal result (supersedes the live status above, not its evidence)

On 2026-09-12 at 07:32:15 AEST, the service terminated with Result=oom-kill,
ExecMainCode=2 and ExecMainStatus=9. Terminal service runtime was 531.524 seconds,
CPU time 565.074 seconds, MemoryPeak=8589934592 bytes and MemorySwapPeak=0.
MemoryMax=8589934592, MemorySwapMax=0 and RuntimeMaxUSec=15min remained enforced.
The process was killed, so its initial running JSON could not be finalized by
on.exit. That original receipt remains intact; this externally observed terminal
record prevents it being mistaken for a live or successful run. No retry occurred.

The retained report is 11862665965 bytes. This 81920000-row recipe is a failed
capacity preflight, not a qualified in-memory control or a finite timing
denominator. Its declaration and failure remain in the series. The five smaller
large-input cases passed resource feasibility only. No candidate performance
measurement, reduced scientific population or increased resource limit follows
from this failure.

```json
[
  {
    "recipe": "annotation_heavy",
    "rows": 3840000,
    "report_bytes": 556026987,
    "peak_bytes": 2772549632,
    "swap_current_bytes": 0,
    "receipt": "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/annotation_heavy/preflight.json",
    "receipt_sha256": "2fad26cfec6e413f67d3641abf0d2fd1bee23cd11515cbc5ffbe9a682cd7532a",
    "helper_sha256": "29ba69bc649f5f273e1bf02a9ac3ccd7c3eb6cddae7b21a8df42ef8ad04ab47e",
    "scenarios_sha256": "54f64237b6d4ae554bfa26b199a5988844d507279ea2a60f2fb8029b10d539ea"
  },
  {
    "recipe": "cycles",
    "rows": 3840000,
    "report_bytes": 556029924,
    "peak_bytes": 2621620224,
    "swap_current_bytes": 0,
    "receipt": "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/cycles/preflight.json",
    "receipt_sha256": "9baca2942b0bbd37ed99f3a1d5d9d318db1ee4afa029e64798319d4e09256bce",
    "helper_sha256": "29ba69bc649f5f273e1bf02a9ac3ccd7c3eb6cddae7b21a8df42ef8ad04ab47e",
    "scenarios_sha256": "54f64237b6d4ae554bfa26b199a5988844d507279ea2a60f2fb8029b10d539ea"
  },
  {
    "recipe": "representative",
    "rows": 3840000,
    "report_bytes": 556028391,
    "peak_bytes": null,
    "swap_current_bytes": null,
    "receipt": "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/representative/preflight.json",
    "receipt_sha256": "49988d9a6624066424e0bf39ae434daf8937eadb4536012054af2446f7dff435",
    "helper_sha256": null,
    "scenarios_sha256": null
  },
  {
    "recipe": "tall",
    "rows": 3456000,
    "report_bytes": 500607058,
    "peak_bytes": 2562809856,
    "swap_current_bytes": 0,
    "receipt": "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/tall/preflight.json",
    "receipt_sha256": "642ad0d366524e574071cb4737a052d186ab5bec7f294f4474e29418546e120a",
    "helper_sha256": "29ba69bc649f5f273e1bf02a9ac3ccd7c3eb6cddae7b21a8df42ef8ad04ab47e",
    "scenarios_sha256": "54f64237b6d4ae554bfa26b199a5988844d507279ea2a60f2fb8029b10d539ea"
  },
  {
    "recipe": "wide",
    "rows": 8192000,
    "report_bytes": 1185864748,
    "peak_bytes": 4147425280,
    "swap_current_bytes": 0,
    "receipt": "/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/wide/preflight.json",
    "receipt_sha256": "42d2728b1f30ee03cecaee3f9e4ddca29cf667329249725ca0cb30425506b329",
    "helper_sha256": "29ba69bc649f5f273e1bf02a9ac3ccd7c3eb6cddae7b21a8df42ef8ad04ab47e",
    "scenarios_sha256": "54f64237b6d4ae554bfa26b199a5988844d507279ea2a60f2fb8029b10d539ea"
  }
]
```
