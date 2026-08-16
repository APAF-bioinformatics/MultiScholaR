test_that("categorical colour rules cover categories, missing values, and machines", {
  metadata <- data.frame(
    Group = c("Control", "Treatment", NA),
    Machine = c("M01", "M02", "M01"),
    stringsAsFactors = FALSE
  )

  rules <- getCategoricalColourRules(
    metadata_tbl = metadata,
    metadata_column_labels = c("Study group", "Instrument"),
    metadata_column_selected = c("Group", "Machine"),
    categorical_columns = c("Group", "Machine"),
    ms_machine_column = "Machine",
    columns_to_exclude = character()
  )

  expect_named(rules, c("Study group", "Instrument"))
  expect_named(rules[["Study group"]], c("Control", "Treatment", "NA"))
  expect_identical(
    unname(rules[["Instrument"]][c("M01", "M02")]),
    unname(getCmriMachineColour()[c("M01", "M02")])
  )
})

test_that("excluded machine metadata does not create a colour rule", {
  rules <- getCategoricalColourRules(
    metadata_tbl = data.frame(Group = "Control", Machine = "M01"),
    metadata_column_labels = c("Study group", "Instrument"),
    metadata_column_selected = c("Group", "Machine"),
    categorical_columns = c("Group", "Machine"),
    ms_machine_column = "Machine",
    columns_to_exclude = "Machine"
  )

  expect_named(rules, "Study group")
})
