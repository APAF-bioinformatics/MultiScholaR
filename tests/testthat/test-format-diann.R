library(testthat)
library(MultiScholaR)

test_that("formatDIANN correctly transforms a mock EList into DIA-NN format", {
  # Create a mock EList
  mock_intensity <- matrix(
    c(10, 11, 12, 13), 
    nrow = 2, 
    dimnames = list(c("PEPTIDE1", "PEPTIDE2"), c("SampleA", "SampleB"))
  )
  
  mock_genes <- data.frame(
    Precursor.Id = c("PEPTIDE1", "PEPTIDE2"),
    Protein.Group = c("PROT1", "PROT2"),
    Protein.Names = c("Name1", "Name2"),
    Genes = c("Gene1", "Gene2"),
    Proteotypic = c(1, 1),
    stringsAsFactors = FALSE
  )
  
  mock_elist <- list(
    E = mock_intensity,
    genes = mock_genes
  )
  
  # Run formatDIANN
  result <- formatDIANN(mock_elist)
  
  # Check results
  expect_s3_class(result, "data.frame")
  expect_true(nrow(result) == 4)  # 2 peptides * 2 samples
  expect_true("Precursor.Id" %in% names(result))
  expect_true("Precursor.Quantity" %in% names(result))
  expect_true("Run" %in% names(result))
  
  # Check linear quantity conversion (2^10 = 1024)
  sample_a_peptide_1 <- result[result$Precursor.Id == "PEPTIDE1" & result$Run == "SampleA", ]
  expect_equal(sample_a_peptide_1$Precursor.Quantity, 2^10)
  
  # Verify both functions yield same result
  result_parquet <- formatDIANNParquet(mock_elist)
  expect_equal(result, result_parquet)
})

test_that("formatDIANN handles missing Precursor.Id in gene info", {
  mock_intensity <- matrix(
    c(10, 11, 12, 13), 
    nrow = 2, 
    dimnames = list(c("PEPTIDE1", "PEPTIDE2"), c("SampleA", "SampleB"))
  )
  
  # Precursor.Id is missing from gene_info, but rows match
  mock_genes <- data.frame(
    Protein.Group = c("PROT1", "PROT2"),
    Protein.Names = c("Name1", "Name2"),
    Genes = c("Gene1", "Gene2"),
    Proteotypic = c(1, 1),
    stringsAsFactors = FALSE
  )
  
  mock_elist <- list(
    E = mock_intensity,
    genes = mock_genes
  )
  
  # This should trigger the defensive logic in formatDIANN
  result <- formatDIANN(mock_elist)
  
  expect_true("Precursor.Id" %in% names(result))
  expect_true(all(result$Precursor.Id %in% c("PEPTIDE1", "PEPTIDE2")))
})
