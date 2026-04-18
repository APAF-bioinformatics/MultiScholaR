# Shared assay input panel builder kept in the wrapper until the next staged UI wave.
buildLipidImportAssayInputPanel <- function(
  ns,
  assayNameInputId,
  assayNameLabel,
  assayNameValue = "",
  assayNamePlaceholder = NULL,
  shinyFilesButtonId,
  standardInputId,
  pathOutputId,
  useShinyFiles,
  fileAccept = c(".tsv", ".tab", ".txt", ".csv", ".xlsx", ".parquet"),
  buildTextInput = shiny::textInput,
  buildShinyFilesButton = shinyFiles::shinyFilesButton,
  buildPathOutput = shiny::verbatimTextOutput,
  buildFileInput = shiny::fileInput
) {
  file_input_ui <- if (useShinyFiles) {
    shiny::tagList(
      buildShinyFilesButton(
        id = ns(shinyFilesButtonId),
        label = "Select File",
        title = "Choose Data File",
        multiple = FALSE,
        icon = shiny::icon("file")
      ),
      shiny::br(),
      buildPathOutput(
        outputId = ns(pathOutputId),
        placeholder = TRUE
      )
    )
  } else {
    buildFileInput(
      inputId = ns(standardInputId),
      label = NULL,
      accept = fileAccept
    )
  }

  shiny::wellPanel(
    style = "background-color: #f8f9fa;",
    shiny::fluidRow(
      shiny::column(
        4,
        buildTextInput(
          ns(assayNameInputId),
          assayNameLabel,
          value = assayNameValue,
          placeholder = assayNamePlaceholder
        )
      ),
      shiny::column(
        8,
        file_input_ui
      )
    )
  )
}

buildLipidImportFileImportSection <- function(
  ns,
  useShinyFiles,
  buildAssayInputPanel = buildLipidImportAssayInputPanel,
  buildRadioButtons = shiny::radioButtons
) {
  shiny::tagList(
    shiny::h4("Step 1: Select Vendor Format"),
    buildRadioButtons(
      ns("vendor_format"),
      NULL,
      choices = c(
        "MS-DIAL" = "msdial",
        "Progenesis QI" = "progenesis",
        "XCMS" = "xcms",
        "Compound Discoverer" = "compound_discoverer",
        "LipidSearch" = "lipidsearch",
        "Other/Custom" = "custom"
      ),
      selected = "msdial",
      inline = TRUE
    ),
    shiny::hr(),
    shiny::h4("Step 2: Import Data Files"),
    shiny::p("Import one or more assay files (e.g., LCMS_Pos, LCMS_Neg, GCMS)"),
    buildAssayInputPanel(
      ns = ns,
      assayNameInputId = "assay1_name",
      assayNameLabel = "Assay Name",
      assayNameValue = "LCMS_Pos",
      shinyFilesButtonId = "assay1_file",
      standardInputId = "assay1_file_std",
      pathOutputId = "assay1_path",
      useShinyFiles = useShinyFiles
    ),
    buildAssayInputPanel(
      ns = ns,
      assayNameInputId = "assay2_name",
      assayNameLabel = "Assay Name (Optional)",
      assayNameValue = "",
      assayNamePlaceholder = "e.g., LCMS_Neg",
      shinyFilesButtonId = "assay2_file",
      standardInputId = "assay2_file_std",
      pathOutputId = "assay2_path",
      useShinyFiles = useShinyFiles
    )
  )
}

# Column-mapping shell kept in the wrapper until the next staged UI wave.
buildLipidImportColumnMappingSection <- function(ns) {
  shiny::tagList(
    shiny::h4("Step 3: Column Mapping"),
    shiny::p("Configure column mappings for your data. Dropdowns are auto-populated based on detected format."),
    shiny::conditionalPanel(
      condition = paste0("output['", ns("file_loaded"), "']"),
      shiny::wellPanel(
        # Format detection status (hidden for custom)
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] != 'custom'"),
          shiny::uiOutput(ns("format_detection_status")),
          shiny::hr()
        ),

        # Custom format instructions
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
          shiny::div(
            class = "alert alert-info",
            shiny::icon("edit"),
            shiny::strong(" Custom Format"),
            shiny::br(),
            "Type the exact column names from your file below."
          ),
          shiny::hr()
        ),

        # Lipid ID column - Dropdown mode (non-custom)
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] != 'custom'"),
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::selectInput(
                ns("lipid_id_col"),
                "Lipid ID Column",
                choices = NULL
              )
            ),
            shiny::column(
              6,
              shiny::uiOutput(ns("lipid_id_status"))
            )
          )
        ),

        # Lipid ID column - Text input mode (custom)
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::textInput(
                ns("lipid_id_col_custom"),
                "Lipid ID Column",
                value = "",
                placeholder = "e.g., Compound_ID"
              )
            ),
            shiny::column(
              6,
              shiny::uiOutput(ns("lipid_id_status_custom"))
            )
          )
        ),

        # Annotation column - Dropdown mode (non-custom)
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] != 'custom'"),
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::selectInput(
                ns("annotation_col"),
                "Annotation Column",
                choices = NULL
              )
            ),
            shiny::column(
              6,
              shiny::uiOutput(ns("annotation_status"))
            )
          )
        ),

        # Annotation column - Text input mode (custom)
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
          shiny::fluidRow(
            shiny::column(
              6,
              shiny::textInput(
                ns("annotation_col_custom"),
                "Annotation Column (optional)",
                value = "",
                placeholder = "e.g., Lipid_Name"
              )
            ),
            shiny::column(
              6,
              shiny::uiOutput(ns("annotation_status_custom"))
            )
          )
        ),

        # Sample columns pattern (custom mode only)
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
          shiny::textInput(
            ns("sample_cols_pattern"),
            "Sample Column Pattern (Regex)",
            value = "",
            placeholder = "e.g., ^Sample_|_Rep[0-9]+$"
          ),
          shiny::helpText("Leave blank to auto-detect numeric columns as samples")
        ),

        # Exclude normalized columns checkbox
        shiny::checkboxInput(
          ns("exclude_norm"),
          "Exclude Vendor-Normalized Columns",
          value = TRUE
        ),
        shiny::helpText("Excludes columns ending in _Norm, _Normalized, etc."),

        # Sanitization checkbox
        shiny::checkboxInput(
          ns("sanitize_names"),
          "Sanitize Sample Names",
          value = TRUE
        ),
        shiny::helpText("Clean sample IDs (e.g., '123-Sample!' -> 'x123_sample') for better compatibility with downstream analysis."),

        # Internal standard pattern
        shiny::textInput(
          ns("is_pattern"),
          "Internal Standard Pattern (Regex)",
          value = "",
          placeholder = "e.g., ^IS_|_d[0-9]+$|ISTD"
        ),
        shiny::helpText("Regular expression to identify internal standards"),
        shiny::hr(),

        # Sample column detection
        shiny::h5("Detected Sample Columns"),
        shiny::verbatimTextOutput(ns("sample_columns_display")),

        # Available columns helper (custom mode)
        shiny::conditionalPanel(
          condition = paste0("input['", ns("vendor_format"), "'] == 'custom'"),
          shiny::hr(),
          shiny::h5("Available Columns in File"),
          shiny::verbatimTextOutput(ns("available_columns_display"))
        )
      ),
      shiny::hr(),

      # Validation summary
      shiny::h4("Step 4: Validation"),
      shiny::uiOutput(ns("validation_summary"))
    ),
    shiny::conditionalPanel(
      condition = paste0("!output['", ns("file_loaded"), "']"),
      shiny::div(
        class = "alert alert-info",
        shiny::icon("info-circle"),
        " Import a data file to configure column mappings."
      )
    )
  )
}

# Process-button footer shell kept in the wrapper until the next staged UI wave.
buildLipidImportProcessFooterSection <- function(
  ns,
  buildActionButton = shiny::actionButton,
  buildUiOutput = shiny::uiOutput
) {
  shiny::fluidRow(
    shiny::column(
      12,
      buildActionButton(
        ns("process_import"),
        "Process Imported Data",
        class = "btn-success",
        width = "100%",
        icon = shiny::icon("check")
      ),
      shiny::br(),
      shiny::br(),
      buildUiOutput(ns("import_status"))
    )
  )
}

# Main module panel shell kept in the wrapper until the next staged UI wave.
buildLipidImportModulePanel <- function(
  ns,
  useShinyFiles,
  buildFileImportSection = buildLipidImportFileImportSection,
  buildColumnMappingSection = buildLipidImportColumnMappingSection,
  buildProcessFooterSection = buildLipidImportProcessFooterSection,
  buildWellPanel = shiny::wellPanel,
  buildHeader = shiny::h3,
  buildFluidRow = shiny::fluidRow,
  buildColumn = shiny::column,
  buildHr = shiny::hr
) {
  buildWellPanel(
    buildHeader("Lipidomics Data Import"),
    buildFluidRow(
      buildColumn(
        6,
        buildFileImportSection(
          ns = ns,
          useShinyFiles = useShinyFiles
        )
      ),
      buildColumn(
        6,
        buildColumnMappingSection(
          ns = ns
        )
      )
    ),
    buildHr(),
    buildProcessFooterSection(
      ns = ns
    )
  )
}

# Outer UI shell kept in the wrapper until the next staged UI review/apply.
buildLipidImportUiShell <- function(
  ns,
  useShinyFiles,
  buildModulePanel = buildLipidImportModulePanel,
  buildTagList = shiny::tagList,
  useShinyjs = shinyjs::useShinyjs,
  buildFluidRow = shiny::fluidRow,
  buildColumn = shiny::column
) {
  buildTagList(
    useShinyjs(),
    buildFluidRow(
      buildColumn(
        12,
        buildModulePanel(
          ns = ns,
          useShinyFiles = useShinyFiles
        )
      )
    )
  )
}

