omicEnrichW113Frozen[["gprofiler"]] <- list(
    history = c("initial", "correlation_filtered", "enrichment_completed"),
    state = c(
        `config$selected_contrast` = "B_vs_A", `config$analysis_method` = "gprofiler2",
        `config$organism_supported` = "TRUE", `config$up_cutoff` = "0", `config$down_cutoff` = "0",
        `config$q_cutoff` = "0.05", `config$organism_taxid` = "9606",
        `config$pathway_dir` = "<ROOT>/results/pathway_enrichment",
        `config$outcomes` = "<digest> 83c639e0603586f8",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A",
        previous = "correlation_filtered", products = ""),
    results = c(
        class = "EnrichmentResults", enrichment_data = "ab5dd21b32dbb33d",
        enrichment_summaries = "6122a4961d0ad898", contrasts = "950ddbb7d830df8c", plots = "e7dada16ba17cceb"),
    products = c(
        type = "dia-enrichment-results", status = "available", value = "e3f94d9e3b685fcb",
        identity = "ad3ef17b26bf1067",
        roots = "checkpoint:source|result_table:da_table:groupB-groupA|result_table:other|result_table:other|result_table:other"),
    provenance = c(
        operation = "enrichment", variant = "whole",
        parameter_names = "taxon_id|up_cutoff|down_cutoff|q_cutoff|exclude_iea|correction_method|protein_id_column|contrast_names|results",
        rest = "032782518e99aba7"),
    record = c(
        state = "enrichment_completed", operation = "enrichment", variant = "whole",
        step_label = "enrichment_completed",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A", products = "0",
        `parameters$taxon_id` = "provenance ; input ; 9606",
        `parameters$up_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$down_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$q_cutoff` = "provenance ; no_shipped_default ; 0.050000000000000003",
        `parameters$exclude_iea` = "provenance ; no_shipped_default ; FALSE",
        `parameters$correction_method` = "provenance ; no_shipped_default ; gSCS",
        `parameters$protein_id_column` = "provenance ; input ; gene_name",
        `parameters$contrast_names` = "provenance ; input ; groupB-groupA",
        `parameters$results` = "provenance ; input ; groupB-groupA=<id>"))

omicEnrichW113Frozen[["clusterprofiler"]] <- list(
    history = c("initial", "correlation_filtered", "enrichment_completed"),
    state = c(
        `config$selected_contrast` = "B_vs_A", `config$analysis_method` = "clusterprofiler",
        `config$organism_supported` = "FALSE", `config$up_cutoff` = "0", `config$down_cutoff` = "0",
        `config$q_cutoff` = "0.05", `config$organism_taxid` = "999999",
        `config$pathway_dir` = "<ROOT>/results/pathway_enrichment",
        `config$outcomes` = "<digest> 78f2d3a35ecc6ef7",
        description = "Enrichment analysis completed using clusterprofiler for contrast: B_vs_A",
        previous = "correlation_filtered", products = ""),
    results = c(
        class = "EnrichmentResults", enrichment_data = "f3960c635154ee90",
        enrichment_summaries = "5d6736b697306a94", contrasts = "950ddbb7d830df8c", plots = "e7dada16ba17cceb"),
    products = c(
        type = "dia-enrichment-results", status = "available", value = "b85cc85cec7a640b",
        identity = "ad3ef17b26bf1067",
        roots = "checkpoint:source|result_table:da_table:groupB-groupA|result_table:other|result_table:other|result_table:other"),
    provenance = c(
        operation = "enrichment", variant = "whole",
        parameter_names = "taxon_id|up_cutoff|down_cutoff|q_cutoff|exclude_iea|correction_method|protein_id_column|contrast_names|results|go_propagation|significance_rule",
        rest = "cd8fd89d6f6482f1"),
    record = c(
        state = "enrichment_completed", operation = "enrichment", variant = "whole",
        step_label = "enrichment_completed",
        description = "Enrichment analysis completed using clusterprofiler for contrast: B_vs_A",
        products = "0", `parameters$taxon_id` = "provenance ; input ; 999999",
        `parameters$up_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$down_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$q_cutoff` = "provenance ; no_shipped_default ; 0.050000000000000003",
        `parameters$exclude_iea` = "provenance ; no_shipped_default ; FALSE",
        `parameters$correction_method` = "provenance ; no_shipped_default ; BH",
        `parameters$protein_id_column` = "provenance ; input ; uniprot_acc",
        `parameters$contrast_names` = "provenance ; input ; groupB-groupA",
        `parameters$results` = "provenance ; input ; groupB-groupA=<id>",
        `parameters$go_propagation` = "provenance ; no_shipped_default ; is_a_part_of_regulates",
        `parameters$significance_rule` = "provenance ; no_shipped_default ; bh_whole_family"))

omicEnrichW113Frozen[["kernel_gprofiler"]] <- list(
    history = c("initial", "correlation_filtered", "enrichment_completed"),
    state = c(
        `config$selected_contrast` = "B_vs_A", `config$analysis_method` = "gprofiler2",
        `config$organism_supported` = "TRUE", `config$up_cutoff` = "0", `config$down_cutoff` = "0",
        `config$q_cutoff` = "0.05", `config$organism_taxid` = "9606",
        `config$pathway_dir` = "<ROOT>/results/pathway_enrichment",
        `config$outcomes` = "<digest> 83c639e0603586f8",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A",
        previous = "correlation_filtered", products = ""),
    results = c(
        class = "EnrichmentResults", enrichment_data = "962ee51f2c925016",
        enrichment_summaries = "18ee1ae7f00ccf86", contrasts = "950ddbb7d830df8c", plots = "28df520b992b4a4b"),
    products = c(
        type = "dia-enrichment-results", status = "available", value = "387eecf6099e563f",
        identity = "ad3ef17b26bf1067",
        roots = "checkpoint:source|result_table:da_table:groupB-groupA|result_table:other|result_table:other|result_table:other"),
    provenance = c(
        operation = "enrichment", variant = "whole",
        parameter_names = "taxon_id|up_cutoff|down_cutoff|q_cutoff|exclude_iea|correction_method|protein_id_column|contrast_names|results",
        rest = "032782518e99aba7"),
    record = c(
        state = "enrichment_completed", operation = "enrichment", variant = "whole",
        step_label = "enrichment_completed",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A", products = "0",
        `parameters$taxon_id` = "provenance ; input ; 9606",
        `parameters$up_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$down_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$q_cutoff` = "provenance ; no_shipped_default ; 0.050000000000000003",
        `parameters$exclude_iea` = "provenance ; no_shipped_default ; FALSE",
        `parameters$correction_method` = "provenance ; no_shipped_default ; gSCS",
        `parameters$protein_id_column` = "provenance ; input ; gene_name",
        `parameters$contrast_names` = "provenance ; input ; groupB-groupA",
        `parameters$results` = "provenance ; input ; groupB-groupA=<id>"))

omicEnrichW113Frozen[["kernel_clusterprofiler"]] <- list(
    history = c("initial", "correlation_filtered", "enrichment_completed"),
    state = c(
        `config$selected_contrast` = "B_vs_A", `config$analysis_method` = "clusterprofiler",
        `config$organism_supported` = "FALSE", `config$up_cutoff` = "0", `config$down_cutoff` = "0",
        `config$q_cutoff` = "0.05", `config$organism_taxid` = "999999",
        `config$pathway_dir` = "<ROOT>/results/pathway_enrichment",
        `config$outcomes` = "<digest> 78f2d3a35ecc6ef7",
        description = "Enrichment analysis completed using clusterprofiler for contrast: B_vs_A",
        previous = "correlation_filtered", products = ""),
    results = c(
        class = "EnrichmentResults", enrichment_data = "609cb7b1bc247d2f",
        enrichment_summaries = "28df520b992b4a4b", contrasts = "950ddbb7d830df8c", plots = "28df520b992b4a4b"),
    products = c(
        type = "dia-enrichment-results", status = "available", value = "1186145dadf9295c",
        identity = "ad3ef17b26bf1067",
        roots = "checkpoint:source|result_table:da_table:groupB-groupA|result_table:other|result_table:other|result_table:other"),
    provenance = c(
        operation = "enrichment", variant = "whole",
        parameter_names = "taxon_id|up_cutoff|down_cutoff|q_cutoff|exclude_iea|correction_method|protein_id_column|contrast_names|results|go_propagation|significance_rule",
        rest = "cd8fd89d6f6482f1"),
    record = c(
        state = "enrichment_completed", operation = "enrichment", variant = "whole",
        step_label = "enrichment_completed",
        description = "Enrichment analysis completed using clusterprofiler for contrast: B_vs_A",
        products = "0", `parameters$taxon_id` = "provenance ; input ; 999999",
        `parameters$up_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$down_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$q_cutoff` = "provenance ; no_shipped_default ; 0.050000000000000003",
        `parameters$exclude_iea` = "provenance ; no_shipped_default ; FALSE",
        `parameters$correction_method` = "provenance ; no_shipped_default ; BH",
        `parameters$protein_id_column` = "provenance ; input ; uniprot_acc",
        `parameters$contrast_names` = "provenance ; input ; groupB-groupA",
        `parameters$results` = "provenance ; input ; groupB-groupA=<id>",
        `parameters$go_propagation` = "provenance ; no_shipped_default ; is_a_part_of_regulates",
        `parameters$significance_rule` = "provenance ; no_shipped_default ; bh_whole_family"))

omicEnrichW113Frozen[["organism_filter/no_mapping"]] <- list(
    history = c("initial", "correlation_filtered", "enrichment_completed"),
    state = c(
        `config$selected_contrast` = "B_vs_A", `config$analysis_method` = "gprofiler2",
        `config$organism_supported` = "TRUE", `config$up_cutoff` = "0", `config$down_cutoff` = "0",
        `config$q_cutoff` = "0.05", `config$organism_taxid` = "9606",
        `config$pathway_dir` = "<ROOT>/results/pathway_enrichment",
        `config$outcomes` = "<digest> a3fad64ece50efbd",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A",
        previous = "correlation_filtered", products = ""),
    results = c(
        class = "EnrichmentResults", enrichment_data = "ab5dd21b32dbb33d",
        enrichment_summaries = "6122a4961d0ad898", contrasts = "950ddbb7d830df8c", plots = "e7dada16ba17cceb"),
    products = c(
        type = "dia-enrichment-results", status = "available", value = "64dbb43f7ac98b4f",
        identity = "ad3ef17b26bf1067",
        roots = "checkpoint:source|result_table:da_table:groupB-groupA|result_table:other|result_table:other|result_table:other"),
    provenance = c(
        operation = "enrichment", variant = "whole",
        parameter_names = "taxon_id|up_cutoff|down_cutoff|q_cutoff|exclude_iea|correction_method|protein_id_column|contrast_names|results|organism_filter_enabled|organism_filter_taxon",
        rest = "db53d45b43df9549"),
    record = c(
        state = "enrichment_completed", operation = "enrichment", variant = "whole",
        step_label = "enrichment_completed",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A", products = "0",
        `parameters$taxon_id` = "provenance ; input ; 9606",
        `parameters$up_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$down_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$q_cutoff` = "provenance ; no_shipped_default ; 0.050000000000000003",
        `parameters$exclude_iea` = "provenance ; no_shipped_default ; FALSE",
        `parameters$correction_method` = "provenance ; no_shipped_default ; gSCS",
        `parameters$protein_id_column` = "provenance ; input ; gene_name",
        `parameters$contrast_names` = "provenance ; input ; groupB-groupA",
        `parameters$results` = "provenance ; input ; groupB-groupA=<id>",
        `parameters$organism_filter_enabled` = "provenance ; no_shipped_default ; TRUE",
        `parameters$organism_filter_taxon` = "provenance ; input ; 9606"))

omicEnrichW113Frozen[["organism_filter/mapped"]] <- list(
    history = c("initial", "correlation_filtered", "enrichment_completed"),
    state = c(
        `config$selected_contrast` = "B_vs_A", `config$analysis_method` = "gprofiler2",
        `config$organism_supported` = "TRUE", `config$up_cutoff` = "0", `config$down_cutoff` = "0",
        `config$q_cutoff` = "0.05", `config$organism_taxid` = "9606",
        `config$pathway_dir` = "<ROOT>/results/pathway_enrichment",
        `config$outcomes` = "<digest> e219d467ecf658c4",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A",
        previous = "correlation_filtered", products = ""),
    results = c(
        class = "EnrichmentResults", enrichment_data = "ab5dd21b32dbb33d",
        enrichment_summaries = "6122a4961d0ad898", contrasts = "950ddbb7d830df8c", plots = "e7dada16ba17cceb"),
    products = c(
        type = "dia-enrichment-results", status = "available", value = "64dbb43f7ac98b4f",
        identity = "ad3ef17b26bf1067",
        roots = "checkpoint:source|result_table:da_table:groupB-groupA|result_table:other|result_table:other|result_table:other"),
    provenance = c(
        operation = "enrichment", variant = "whole",
        parameter_names = "taxon_id|up_cutoff|down_cutoff|q_cutoff|exclude_iea|correction_method|protein_id_column|contrast_names|results|organism_filter_enabled|organism_filter_taxon",
        rest = "db53d45b43df9549"),
    record = c(
        state = "enrichment_completed", operation = "enrichment", variant = "whole",
        step_label = "enrichment_completed",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A", products = "0",
        `parameters$taxon_id` = "provenance ; input ; 9606",
        `parameters$up_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$down_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$q_cutoff` = "provenance ; no_shipped_default ; 0.050000000000000003",
        `parameters$exclude_iea` = "provenance ; no_shipped_default ; FALSE",
        `parameters$correction_method` = "provenance ; no_shipped_default ; gSCS",
        `parameters$protein_id_column` = "provenance ; input ; gene_name",
        `parameters$contrast_names` = "provenance ; input ; groupB-groupA",
        `parameters$results` = "provenance ; input ; groupB-groupA=<id>",
        `parameters$organism_filter_enabled` = "provenance ; no_shipped_default ; TRUE",
        `parameters$organism_filter_taxon` = "provenance ; input ; 9606"))

omicEnrichW113Frozen[["refusals/pre_da"]] <- list(
    history = c("initial", "correlation_filtered"),
    state = "<none>",
    results = "<none>",
    products = "<none>",
    provenance = "<none>",
    record = "<none>")

omicEnrichW113Frozen[["refusals/q_cutoff"]] <- list(
    history = c("initial", "correlation_filtered"),
    state = "<none>",
    results = "<none>",
    products = "<none>",
    provenance = "<none>",
    record = "<none>")

omicEnrichW113Frozen[["refusals/table_removed"]] <- list(
    history = c("initial", "correlation_filtered"),
    state = "<none>",
    results = "<none>",
    products = "<none>",
    provenance = "<none>",
    record = "<none>")

omicEnrichW113Frozen[["restore"]] <- list(
    history = c("initial", "correlation_filtered", "enrichment_completed"),
    state = c(
        `config$selected_contrast` = "B_vs_A", `config$analysis_method` = "gprofiler2",
        `config$organism_supported` = "TRUE", `config$up_cutoff` = "0", `config$down_cutoff` = "0",
        `config$q_cutoff` = "0.05", `config$organism_taxid` = "9606",
        `config$pathway_dir` = "<ROOT>/results/pathway_enrichment",
        `config$outcomes` = "<digest> 83c639e0603586f8",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A",
        previous = "correlation_filtered", products = ""),
    results = c(
        class = "EnrichmentResults", enrichment_data = "ab5dd21b32dbb33d",
        enrichment_summaries = "6122a4961d0ad898", contrasts = "950ddbb7d830df8c", plots = "e7dada16ba17cceb"),
    products = c(
        type = "dia-enrichment-results", status = "available", value = "e3f94d9e3b685fcb",
        identity = "ad3ef17b26bf1067",
        roots = "checkpoint:source|result_table:da_table:groupB-groupA|result_table:other|result_table:other|result_table:other"),
    provenance = c(
        operation = "enrichment", variant = "whole",
        parameter_names = "taxon_id|up_cutoff|down_cutoff|q_cutoff|exclude_iea|correction_method|protein_id_column|contrast_names|results",
        rest = "032782518e99aba7"),
    record = c(
        state = "enrichment_completed", operation = "enrichment", variant = "whole",
        step_label = "enrichment_completed",
        description = "Enrichment analysis completed using gprofiler2 for contrast: B_vs_A", products = "0",
        `parameters$taxon_id` = "provenance ; input ; 9606",
        `parameters$up_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$down_cutoff` = "provenance ; no_shipped_default ; 0",
        `parameters$q_cutoff` = "provenance ; no_shipped_default ; 0.050000000000000003",
        `parameters$exclude_iea` = "provenance ; no_shipped_default ; FALSE",
        `parameters$correction_method` = "provenance ; no_shipped_default ; gSCS",
        `parameters$protein_id_column` = "provenance ; input ; gene_name",
        `parameters$contrast_names` = "provenance ; input ; groupB-groupA",
        `parameters$results` = "provenance ; input ; groupB-groupA=<id>"))

omicEnrichW113Frozen[["prefill"]] <- list(
    history = c("initial", "correlation_filtered"),
    state = "<none>",
    results = "<none>",
    products = "<none>",
    provenance = "<none>",
    record = "<none>")

