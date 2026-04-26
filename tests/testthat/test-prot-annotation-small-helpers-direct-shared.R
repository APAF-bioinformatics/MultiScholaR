# fidelity-coverage-compare: shared
library(testthat)

makeFunctionWithOverrides <- function(fun, replacements) {
  fun_override <- fun
  environment(fun_override) <- list2env(replacements, parent = environment(fun))
  fun_override
}

localNamespaceBinding <- function(name, value, package = "MultiScholaR", env = parent.frame()) {
  ns <- asNamespace(package)
  had_value <- exists(name, envir = ns, inherits = FALSE)
  old_value <- if (had_value) get(name, envir = ns, inherits = FALSE) else NULL
  was_locked <- had_value && bindingIsLocked(name, ns)

  if (was_locked) {
    unlockBinding(name, ns)
  }
  assign(name, value, envir = ns)
  if (was_locked) {
    lockBinding(name, ns)
  }

  withr::defer({
    if (exists(name, envir = ns, inherits = FALSE) && bindingIsLocked(name, ns)) {
      unlockBinding(name, ns)
    }
    if (had_value) {
      assign(name, old_value, envir = ns)
    } else if (exists(name, envir = ns, inherits = FALSE)) {
      rm(list = name, envir = ns)
    }
    if (was_locked) {
      lockBinding(name, ns)
    }
  }, envir = env)
}

getFastaFields <- get(
  "getFastaFields",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
parseFastaObject <- get(
  "parseFastaObject",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
parseFastaFile <- get(
  "parseFastaFile",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
chooseBestPhosphositeAccession <- get(
  "chooseBestPhosphositeAccession",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
.prepareUniprotBatchInput <- get(
  ".prepareUniprotBatchInput",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
batchQueryEvidenceHelper <- get(
  "batchQueryEvidenceHelper",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
batchQueryEvidence <- get(
  "batchQueryEvidence",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
batchQueryEvidenceHelperGeneId <- get(
  "batchQueryEvidenceHelperGeneId",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
batchQueryEvidenceGeneId <- get(
  "batchQueryEvidenceGeneId",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
appender_shiny <- get(
  "appender_shiny",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
setup_shiny_logger <- get(
  "setup_shiny_logger",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)
.capture_checkpoint <- get(
  ".capture_checkpoint",
  envir = asNamespace("MultiScholaR"),
  inherits = FALSE
)

test_that("protein annotation helper utilities preserve FASTA parsing, batch preparation, and checkpoint logging", {
  expect_identical(
    getFastaFields("sp|P12345|ID OS=Homo sapiens OX=9606", "OS"),
    "Homo sapiens"
  )
  expect_true(is.na(getFastaFields("sp|P12345|ID OX=9606", "GN")))

  aa_seq <- c("MAAAA", "MBBBB")
  names(aa_seq) <- c(
    "sp|P12345-2|ID1 OS=Homo sapiens OX=9606 GN=GENE1 PE=1 SV=2",
    "tr|Q99999|ID2 OS=Mus musculus OX=10090 GN=GENE2 PE=3 SV=1"
  )

  parsed_fasta <- parseFastaObject(aa_seq)
  expect_identical(parsed_fasta$cleaned_acc, c("P12345", "Q99999"))
  expect_identical(as.character(parsed_fasta$is_isoform), c("Isoform", "Canonical"))
  expect_identical(as.character(parsed_fasta$status), c("reviewed", "unreviewed"))

  parse_fasta_file <- makeFunctionWithOverrides(
    parseFastaFile,
    list(
      read.fasta = function(file, seqtype, whole.header, as.string) {
        aa_seq
      },
      parseFastaObject = function(aa_seq) parsed_fasta,
      str_match = stringr::str_match,
      map_chr = purrr::map_chr,
      str_length = stringr::str_length
    )
  )
  fasta_tbl <- parse_fasta_file("fake.fasta")
  expect_identical(fasta_tbl$seq, unname(aa_seq))
  expect_identical(fasta_tbl$seq_length, c(5L, 5L))

  phospho_tbl <- chooseBestPhosphositeAccession(
    input_tbl = data.frame(
      group = c("g1", "g2"),
      accessions = c("P12345-2;Q99999", "Q99999"),
      cleaned_peptide = c("AAAA", "BBBB"),
      stringsAsFactors = FALSE
    ),
    acc_detail_tab = data.frame(
      uniprot_acc = c("P12345-2", "Q99999"),
      cleaned_acc = c("P12345", "Q99999"),
      gene_name = c("GENE1", "GENE2"),
      protein_evidence = c(1, 3),
      status = c("reviewed", "unreviewed"),
      is_isoform = c("Isoform", "Canonical"),
      isoform_num = c(2, 0),
      seq_length = c(500L, 300L),
      seq = c("MAAAAK", "MBBBBK"),
      annotation_score = c(10, 5),
      stringsAsFactors = FALSE
    ),
    accessions_column = accessions,
    group_id = group
  )
  expect_true(all(c("group", "gene_name", "uniprot_acc") %in% colnames(phospho_tbl)))
  expect_identical(phospho_tbl$uniprot_acc, c("P12345-2", "Q99999"))

  prepare_batch <- makeFunctionWithOverrides(
    .prepareUniprotBatchInput,
    list(
      getUniprotRegexPatterns = function() list(entry = "^[PQ][0-9]+(?:-[0-9]+)?$"),
      cleanIsoformNumber = function(x) sub("-[0-9]+$", "", x)
    )
  )

  prepared <- prepare_batch(
    input_tbl = data.frame(ids = c("P12345-2;Q99999", "BAD", "P12345-3"), stringsAsFactors = FALSE),
    id_column = ids,
    batch_size = 1
  )
  expect_identical(prepared$ids, c("P12345", "Q99999"))
  expect_identical(prepared$round, c(1, 2))

  prepared_gene <- prepare_batch(
    input_tbl = data.frame(ids = c("P12345:Q99999", "Q11111"), stringsAsFactors = FALSE),
    id_column = ids,
    delim = ":",
    batch_size = 2,
    clean_isoform = FALSE
  )
  expect_identical(prepared_gene$round, c(1, 1, 2))

  helper_batch <- makeFunctionWithOverrides(
    batchQueryEvidenceHelper,
    list(
      .prepareUniprotBatchInput = prepare_batch
    )
  )
  helper_tbl <- helper_batch(
    uniprot_acc_tbl = data.frame(accession = c("P12345-2;Q99999"), stringsAsFactors = FALSE),
    uniprot_acc_column = accession
  )
  expect_identical(helper_tbl$accession, c("P12345", "Q99999"))

  batch_query <- makeFunctionWithOverrides(
    batchQueryEvidence,
    list(
      batchQueryEvidenceHelper = function(uniprot_acc_tbl, uniprot_acc_column) {
        data.frame(accession = c("P12345", "Q99999"), round = c(1, 2), stringsAsFactors = FALSE)
      },
      subsetQuery = function(data, subset, accessions_col_name, uniprot_handle, uniprot_columns, uniprot_keytype) {
        data.frame(round = subset, accession = data$accession[data$round == subset], stringsAsFactors = FALSE)
      },
      partial = purrr::partial
    )
  )
  query_tbl <- batch_query(
    uniprot_acc_tbl = data.frame(accession = c("P12345;Q99999"), stringsAsFactors = FALSE),
    uniprot_acc_column = accession,
    uniprot_handle = "fake"
  )
  expect_identical(query_tbl$accession, c("P12345", "Q99999"))

  helper_gene <- makeFunctionWithOverrides(
    batchQueryEvidenceHelperGeneId,
    list(
      .prepareUniprotBatchInput = prepare_batch
    )
  )
  helper_gene_tbl <- helper_gene(
    input_tbl = data.frame(gene_id = c("P12345:Q99999"), stringsAsFactors = FALSE),
    gene_id_column = gene_id
  )
  expect_identical(helper_gene_tbl$gene_id, c("P12345", "Q99999"))

  batch_query_gene <- makeFunctionWithOverrides(
    batchQueryEvidenceGeneId,
    list(
      batchQueryEvidenceHelperGeneId = function(input_tbl, gene_id_column, delim = ":") {
        data.frame(gene_id = c("P12345", "Q99999"), round = c(1, 2), stringsAsFactors = FALSE)
      },
      subsetQuery = function(data, subset, accessions_col_name, uniprot_handle, uniprot_columns, uniprot_keytype) {
        data.frame(round = subset, gene_id = data$gene_id[data$round == subset], stringsAsFactors = FALSE)
      },
      partial = purrr::partial
    )
  )
  query_gene_tbl <- batch_query_gene(
    input_tbl = data.frame(gene_id = c("P12345:Q99999"), stringsAsFactors = FALSE),
    gene_id_column = gene_id,
    uniprot_handle = "fake"
  )
  expect_identical(query_gene_tbl$gene_id, c("P12345", "Q99999"))

  shiny_log_store <- local({
    current <- paste(rep("old", 1001), collapse = "\n")
    function(value) {
      if (missing(value)) {
        current
      } else {
        current <<- value
        invisible(current)
      }
    }
  })
  localNamespaceBinding("log_messages", shiny_log_store)
  appender_shiny(c("new entry", "second line"))
  updated_logs <- shiny_log_store()
  expect_true(startsWith(updated_logs, "new entry\nsecond line"))
  expect_lte(length(strsplit(updated_logs, "\n", fixed = TRUE)[[1]]), 1000)

  expect_no_error(setup_shiny_logger())

  checkpoint_dir <- tempfile("checkpoint-dir-")
  withr::defer(unlink(checkpoint_dir, recursive = TRUE, force = TRUE))
  options(
    multischolar.capture_test_checkpoints = TRUE,
    multischolar.checkpoint_dir = checkpoint_dir,
    multischolar.checkpoint_dataset = "test",
    multischolar.checkpoint_omics_layer = "proteomics"
  )
  withr::defer(options(
    multischolar.capture_test_checkpoints = NULL,
    multischolar.checkpoint_dir = NULL,
    multischolar.checkpoint_dataset = NULL,
    multischolar.checkpoint_omics_layer = NULL
  ))

  .capture_checkpoint(
    data = list(value = 1),
    checkpoint_id = "cp01",
    label = "small-helper"
  )
  expect_true(file.exists(file.path(checkpoint_dir, "cp01_small-helper.rds")))
})
