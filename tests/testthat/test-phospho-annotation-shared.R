# fidelity-coverage-compare: shared
library(testthat)

localNamespaceBinding <- function(env, name, value, .local_envir = parent.frame()) {
  had_binding <- exists(name, envir = env, inherits = FALSE)
  old_value <- if (had_binding) get(name, envir = env, inherits = FALSE) else NULL
  was_locked <- had_binding && bindingIsLocked(name, env)

  if (was_locked) {
    unlockBinding(name, env)
  }
  assign(name, value, envir = env)
  if (was_locked) {
    lockBinding(name, env)
  }

  withr::defer({
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (had_binding) {
      assign(name, old_value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (was_locked && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }, envir = .local_envir)
}

localNamespaceBindings <- function(env, bindings, .local_envir = parent.frame()) {
  for (name in names(bindings)) {
    localNamespaceBinding(
      env = env,
      name = name,
      value = bindings[[name]],
      .local_envir = .local_envir
    )
  }
}

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

test_that("phosphosite helper primitives preserve ranking, position, and motif formatting", {
  local_mocked_bindings(
    future_map2 = purrr::map2,
    future_pmap_chr = purrr::pmap_chr,
    .package = "furrr"
  )

  add_columns <- makeFunctionWithOverrides(
    addColumnsToEvidenceTbl,
    list(str_replace_all = stringr::str_replace_all)
  )
  get_max_prob <- makeFunctionWithOverrides(
    getMaxProb,
    list(
      str_match_all = stringr::str_match_all,
      keep = purrr::keep
    )
  )
  get_max_prob_future <- makeFunctionWithOverrides(
    getMaxProbFutureMap,
    list(getMaxProb = get_max_prob)
  )
  get_best_position <- makeFunctionWithOverrides(
    getBestPosition,
    list(
      str_detect = stringr::str_detect,
      str_match_all = stringr::str_match_all,
      str_replace_all = stringr::str_replace_all,
      str_locate_all = stringr::str_locate_all,
      getMaxProb = get_max_prob
    )
  )
  get_best_position_future <- makeFunctionWithOverrides(
    getBestPositionFutureMap,
    list(getBestPosition = get_best_position)
  )
  get_pos_string <- makeFunctionWithOverrides(
    getPosString,
    list(
      cross2 = purrr::cross2,
      map_dbl = purrr::map_dbl
    )
  )
  get_x_mer_string <- makeFunctionWithOverrides(
    getXMerString,
    list(
      str_length = stringr::str_length,
      str_sub = stringr::str_sub
    )
  )
  get_x_mers_list <- makeFunctionWithOverrides(
    getXMersList,
    list(
      cross2 = purrr::cross2,
      map_dbl = purrr::map_dbl,
      getXMerString = get_x_mer_string
    )
  )
  format_position <- makeFunctionWithOverrides(
    formatPhosphositePosition,
    list(
      str_split = stringr::str_split,
      str_replace_all = stringr::str_replace_all
    )
  )

  evidence_tbl <- tibble::tibble(
    phospho_sty_probabilities = c("S(0.90)T(0.80)Y", "A(0.30)B")
  )

  cleaned <- add_columns(evidence_tbl, phospho_sty_probabilities)
  expect_equal(cleaned$evidence_id, c(0L, 1L))
  expect_equal(cleaned$cleaned_peptide, c("STY", "AB"))

  expect_equal(get_max_prob("S(0.90)T(0.80)Y", 2), c(0.9, 0.8))
  expect_null(get_max_prob("ABC", 1))
  expect_equal(get_max_prob_future(c("S(0.90)T(0.80)Y", "A(0.30)B"), c(2, 1)), list(c(0.9, 0.8), 0.3))

  expect_equal(get_best_position("S(0.90)T(0.80)Y", 2), c(1L, 2L))
  expect_error(get_best_position("Sp(0.90)TY", 1), "little 'p'")
  expect_equal(get_best_position_future(c("S(0.90)T(0.80)Y"), c(2)), list(c(1L, 2L)))

  expect_identical(suppressWarnings(get_pos_string(3, c(1, 2))), "3;4")
  expect_identical(suppressWarnings(get_pos_string(c(3, 8), c(1, 2))), "(3;4)|(8;9)")

  expect_identical(get_x_mer_string("AASTYGG", "P1", 1, padding_length = 2), "__AAS")
  expect_identical(get_x_mer_string("AASTYGG", "P1", 4, padding_length = 2), "ASTYG")
  expect_identical(get_x_mers_list("AASTYGG", "P1", 3, c(1, 2), padding_length = 1), "AST;STY")

  expect_identical(format_position("GENE1", "(3;4)|7", "S;T"), "GENE1;S3,T4")
})

test_that("phosphosite filtering and reshaping helpers preserve evidence-to-site tables", {
  local_mocked_bindings(
    future_map2 = purrr::map2,
    future_pmap_chr = purrr::pmap_chr,
    .package = "furrr"
  )

  get_max_prob <- makeFunctionWithOverrides(
    getMaxProb,
    list(
      str_match_all = stringr::str_match_all,
      keep = purrr::keep
    )
  )
  get_max_prob_future <- makeFunctionWithOverrides(
    getMaxProbFutureMap,
    list(getMaxProb = get_max_prob)
  )
  get_best_position <- makeFunctionWithOverrides(
    getBestPosition,
    list(
      str_detect = stringr::str_detect,
      str_match_all = stringr::str_match_all,
      str_replace_all = stringr::str_replace_all,
      str_locate_all = stringr::str_locate_all,
      getMaxProb = get_max_prob
    )
  )
  get_best_position_future <- makeFunctionWithOverrides(
    getBestPositionFutureMap,
    list(getBestPosition = get_best_position)
  )
  get_pos_string <- makeFunctionWithOverrides(
    getPosString,
    list(
      cross2 = purrr::cross2,
      map_dbl = purrr::map_dbl
    )
  )
  get_x_mer_string <- makeFunctionWithOverrides(
    getXMerString,
    list(
      str_length = stringr::str_length,
      str_sub = stringr::str_sub
    )
  )
  get_x_mers_list <- makeFunctionWithOverrides(
    getXMersList,
    list(
      cross2 = purrr::cross2,
      map_dbl = purrr::map_dbl,
      getXMerString = get_x_mer_string
    )
  )
  filter_peptides <- makeFunctionWithOverrides(
    filterPeptideAndExtractProbabilities,
    list(
      str_detect = stringr::str_detect,
      getMaxProbFutureMap = get_max_prob_future,
      map_lgl = purrr::map_lgl,
      map2_lgl = purrr::map2_lgl,
      getBestPositionFutureMap = get_best_position_future
    )
  )
  add_peptide_bounds <- makeFunctionWithOverrides(
    addPeptideStartAndEnd,
    list(
      str_locate_all = stringr::str_locate_all,
      map_chr = purrr::map_chr
    )
  )
  add_position_strings <- makeFunctionWithOverrides(
    addPhosphositesPositionsString,
    list(
      map_chr = purrr::map_chr,
      map2_chr = purrr::map2_chr,
      cross2 = purrr::cross2,
      map_dbl = purrr::map_dbl,
      getPosString = get_pos_string
    )
  )
  add_xmers <- makeFunctionWithOverrides(
    addXMerStrings,
    list(getXMersList = get_x_mers_list)
  )
  filter_by_score <- makeFunctionWithOverrides(
    filterByScoreAndGetSimilarPeptides,
    list(
      map2_lgl = purrr::map2_lgl,
      str_split = stringr::str_split,
      map_int = purrr::map_int,
      str_replace_all = stringr::str_replace_all,
      map_chr = purrr::map_chr
    )
  )
  pivot_longer_sites <- makeFunctionWithOverrides(
    allPhosphositesPivotLonger,
    list(
      as_name = rlang::as_name,
      enquo = rlang::enquo,
      str_replace = stringr::str_replace,
      map_int = purrr::map_int
    )
  )
  group_paralogs <- makeFunctionWithOverrides(
    groupParalogPeptides,
    list(
      as_name = rlang::as_name,
      enquo = rlang::enquo
    )
  )
  pivot_wider_sites <- makeFunctionWithOverrides(
    allPhosphositesPivotWider,
    list(
      as_name = rlang::as_name,
      enquo = rlang::enquo
    )
  )

  evidence_tbl_cleaned <- tibble::tibble(
    evidence_id = c(0L, 1L, 2L, 3L),
    cleaned_peptide = c("STY", "STY", "STY", "AB"),
    phospho_sty_probabilities = c(
      "S(0.90)T(0.80)Y",
      "S(0.85)T(0.60)Y",
      "S(0.95)T(0.92)Y",
      "A(0.30)B"
    ),
    phospho_sty = c(2L, 2L, 2L, 1L),
    leading_proteins = c("P1", "P1", "CON__P3", "P4"),
    experiment = c("exp1", "exp1", "exp2", "exp1"),
    `reporter intensity corrected_1` = c(10, 12, 8, 0),
    `reporter intensity corrected_2` = c(11, 13, 9, 0)
  )
  accession_gene_name_tbl <- tibble::tibble(
    evidence_id = c(0L, 1L),
    uniprot_acc = c("P1", "P1"),
    gene_name = c("GENE1", "GENE1")
  )
  aa_seq_tbl <- tibble::tibble(
    uniprot_acc = "P1",
    seq = "AASTYGG"
  )

  removed <- removePeptidesWithoutAbundances(evidence_tbl_cleaned, "reporter intensity corrected")
  expect_equal(removed$evidence_id, c(0L, 1L, 2L))

  sites_probability_tbl <- filter_peptides(
    evidence_tbl_cleaned = evidence_tbl_cleaned,
    accession_gene_name_tbl = accession_gene_name_tbl,
    col_pattern = "reporter intensity corrected",
    accession_col = leading_proteins,
    phospho_site_prob_col = phospho_sty_probabilities,
    num_phospho_site_col = phospho_sty
  )
  expect_equal(sites_probability_tbl$evidence_id, c(0L, 1L))
  expect_equal(sites_probability_tbl$best_phos_prob[[1]], c(0.9, 0.8))
  expect_equal(sites_probability_tbl$best_phos_pos[[1]], c(1L, 2L))

  peptide_start_and_end <- add_peptide_bounds(sites_probability_tbl, aa_seq_tbl)
  expect_identical(peptide_start_and_end$pep_start[[1]], "3")
  expect_identical(peptide_start_and_end$pep_end[[1]], "5")

  phosphosite_tbl <- add_position_strings(peptide_start_and_end)
  expect_identical(phosphosite_tbl$protein_site_positions[[1]], "3;4")

  xmer_tbl <- add_xmers(phosphosite_tbl, padding_length = 1)
  expect_identical(xmer_tbl$phos_15mer_seq[[1]], "AST;STY")

  filtered_sites <- filter_by_score(
    xmer_tbl,
    site_prob_threshold = 0.75,
    secondary_site_prob_threshold = 0.5,
    num_phospho_site_col = phospho_sty
  )
  filtered_sites <- dplyr::mutate(filtered_sites, sequence = seq)
  expect_equal(filtered_sites$evidence_id, c(0L, 1L))

  all_sites_long <- pivot_longer_sites(
    filtered_sites,
    additional_cols = "experiment",
    col_pattern = "reporter intensity corrected",
    phospho_site_prob_col = phospho_sty_probabilities,
    num_phospho_site_col = phospho_sty
  )
  expect_equal(nrow(all_sites_long), 4L)
  expect_equal(sort(unique(all_sites_long$replicate)), c(1L, 2L))

  paralog_sites_long <- group_paralogs(
    all_sites_long,
    additional_cols = "experiment",
    phospho_site_prob_col = phospho_sty_probabilities,
    num_phospho_site_col = phospho_sty
  )
  expect_true(all(c("gene_names", "uniprot_acc", "phos_15mer_seq") %in% names(paralog_sites_long)))

  all_sites_wide <- pivot_wider_sites(
    paralog_sites_long,
    additional_cols = "experiment",
    phospho_site_prob_col = phospho_sty_probabilities,
    num_phospho_site_col = phospho_sty
  )
  expect_gt(ncol(all_sites_wide), 6L)

  summarised_long <- uniquePhosphositesSummariseLongList(
    paralog_sites_long,
    additional_cols = "experiment"
  )
  expect_equal(sort(names(summarised_long)), c("mean", "median", "sum"))
  expect_equal(length(summarised_long), 3L)

  summarised_wide <- uniquePhosphositesSummariseWideList(
    summarised_long,
    additional_cols = "experiment"
  )
  expect_equal(length(summarised_wide), 3L)
  expect_true(all(c("sites_id", "maxquant_row_ids") %in% names(summarised_wide[[1]])))
})

test_that("multisite orchestration and site-id ranking preserve staged helper wiring", {
  call_log <- character()
  process_multisite <- makeFunctionWithOverrides(
    processMultisiteEvidence,
    list(
      evidence_col_to_use = rlang::sym("phospho_sty_probabilities"),
      parseFastaFile = function(fasta_file) {
        call_log <<- c(call_log, "parseFastaFile")
        tibble::tibble(uniprot_acc = "P1", seq = "AASTYGG")
      },
      addColumnsToEvidenceTbl = function(evidence_tbl, phospho_site_prob_col) {
        call_log <<- c(call_log, "addColumnsToEvidenceTbl")
        tibble::tibble(
          evidence_id = 1L,
          phospho_sty_probabilities = "S(0.90)T(0.80)Y",
          leading_proteins = "P1",
          experiment = "exp1",
          corrected_1 = 10
        )
      },
      chooseBestAccession = function(...) {
        call_log <<- c(call_log, "chooseBestAccession")
        tibble::tibble(evidence_id = 1L, uniprot_acc = "P1", gene_name = "GENE1")
      },
      removePeptidesWithoutAbundances = function(...) {
        call_log <<- c(call_log, "removePeptidesWithoutAbundances")
        tibble::tibble(evidence_id = 1L)
      },
      filterPeptideAndExtractProbabilities = function(...) {
        call_log <<- c(call_log, "filterPeptideAndExtractProbabilities")
        tibble::tibble(evidence_id = 1L, uniprot_acc = "P1", gene_name = "GENE1")
      },
      addPeptideStartAndEnd = function(...) {
        call_log <<- c(call_log, "addPeptideStartAndEnd")
        tibble::tibble(evidence_id = 1L)
      },
      addPhosphositesPositionsString = function(...) {
        call_log <<- c(call_log, "addPhosphositesPositionsString")
        tibble::tibble(evidence_id = 1L)
      },
      addXMerStrings = function(...) {
        call_log <<- c(call_log, "addXMerStrings")
        tibble::tibble(evidence_id = 1L)
      },
      filterByScoreAndGetSimilarPeptides = function(...) {
        call_log <<- c(call_log, "filterByScoreAndGetSimilarPeptides")
        tibble::tibble(evidence_id = 1L)
      },
      allPhosphositesPivotLonger = function(...) {
        call_log <<- c(call_log, "allPhosphositesPivotLonger")
        tibble::tibble(evidence_id = 1L)
      },
      groupParalogPeptides = function(...) {
        call_log <<- c(call_log, "groupParalogPeptides")
        tibble::tibble(evidence_id = 1L)
      },
      allPhosphositesPivotWider = function(...) {
        call_log <<- c(call_log, "allPhosphositesPivotWider")
        tibble::tibble(sites_id = "P1!GENE1!S3!AST")
      },
      uniquePhosphositesSummariseLongList = function(...) {
        call_log <<- c(call_log, "uniquePhosphositesSummariseLongList")
        list(mean = tibble::tibble(sites_id = "P1!GENE1!S3!AST"))
      },
      uniquePhosphositesSummariseWideList = function(...) {
        call_log <<- c(call_log, "uniquePhosphositesSummariseWideList")
        list(mean = tibble::tibble(sites_id = "P1!GENE1!S3!AST"))
      }
    )
  )

  results <- process_multisite(
    fasta_file = "demo.fasta",
    evidence_tbl = tibble::tibble(
      leading_proteins = "P1",
      experiment = "exp1",
      phospho_sty_probabilities = "S(0.90)T(0.80)Y"
    ),
    accession_col = leading_proteins,
    group_id = experiment,
    additional_cols = "experiment",
    col_pattern = "corrected"
  )

  expect_equal(
    call_log,
    c(
      "parseFastaFile",
      "addColumnsToEvidenceTbl",
      "chooseBestAccession",
      "removePeptidesWithoutAbundances",
      "filterPeptideAndExtractProbabilities",
      "addPeptideStartAndEnd",
      "addPhosphositesPositionsString",
      "addXMerStrings",
      "filterByScoreAndGetSimilarPeptides",
      "allPhosphositesPivotLonger",
      "groupParalogPeptides",
      "allPhosphositesPivotWider",
      "uniquePhosphositesSummariseLongList",
      "uniquePhosphositesSummariseWideList"
    )
  )
  expect_true(all(c("summarised_wide_list", "summarised_long_list", "all_phos_sites_wide", "all_phos_sites_long") %in% names(results)))

  rank_fun <- makeFunctionWithOverrides(
    getUniprotAccRankFromSitesId,
    list(str_split = stringr::str_split)
  )

  rank_tbl <- rank_fun(
    tibble::tibble(
      uniprot_acc = c("P1:iso1", "P2:iso2"),
      sites_id = c("P1:P3!GENE1:GENE3!S3!AST", "P4:P2!GENE4:GENE2!T5!STY")
    ),
    uniprot_acc,
    sites_id
  )
  expect_equal(rank_tbl$gene_list_position, c(1L, 2L))
})
