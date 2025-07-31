# ============================================================================
# MultiOmics Workflow Project Setup Script
# ============================================================================
# Hi there! Welcome to MultiScholaR :)
# Instructions: 
# 1. Set your project name and optional directory below
# 2. Select ALL of this code (Ctrl+A or Cmd+A)
# 3. Run it (Ctrl+Enter or Cmd+Enter)
# ============================================================================

# SET YOUR PROJECT OPTIONS HERE:
my_project_name <- "my_analysis"  # Required: Name your project
my_project_dir <- NULL            # Optional: Set specific directory (leave as NULL for default Documents folder)


# Examples for my_project_dir:
# Windows: my_project_dir <- "C:/Users/username/Projects"
# Mac:     my_project_dir <- "~/Projects"

# SET YOUR WORKFLOW OPTIONS HERE:
omic_types <- c("proteomics", "metabolomics", "transcriptomics", "lipidomics", "integration") 
# Options: c("proteomics", "metabolomics", "transcriptomics", ...) 

# The following two options primarily apply when 'proteomics' is in omic_types:
workflow_type <- "DIANN-impute"  # Proteomics-specific: "DIA-NN", "DIANN-impute", "LFQ - FragPipe", "LFQ - MaxQuant", "TMT - MaxQuant", "TMT - FragPipe"
user_experience <- "beginner"  # Proteomics-specific: "experienced", "beginner"

# ============================================================================
# DO NOT MODIFY CODE BELOW THIS LINE
# ============================================================================

# Install required packages if not already installed
if (!requireNamespace("httr", quietly = TRUE)) {
  install.packages("httr")
}
if (!requireNamespace("rstudioapi", quietly = TRUE)) {
  install.packages("rstudioapi")
}
if (!requireNamespace("later", quietly = TRUE)) {
  install.packages("later")
}

# Validate user input
if (!is.character(my_project_name) || nchar(my_project_name) == 0) {
  stop("Please set a valid project name at the top of the script")
}

# --- Helper Functions for URL Generation ---

# Determine workflow file URL based on omic type, workflow type, and user experience
getWorkflowUrl <- function(wf_type, usr_exp, omic) {
  base_url <- paste0("https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/", omic)
  subfolder <- ifelse(usr_exp == "beginner", "starter", "standard")

  if (omic == "proteomics") {
    # Proteomics-specific logic
    if (wf_type == "DIA-NN") {
      # Use specific filenames for DIA-NN based on experience
      filename <- ifelse(usr_exp == "beginner", "DIA_workflow_starter.rmd", "DIA_workflow_experienced.rmd")
      return(paste0(base_url, "/", subfolder, "/", filename))
    } else if (wf_type == "DIANN-impute") {
      # Use specific filenames for DIA-NN based on experience
      filename <- ifelse(usr_exp == "beginner", "DIA_workflow_limpa_starter.rmd", "DIA_workflow_experienced.rmd")
      return(paste0(base_url, "/", subfolder, "/", filename))
    } else if (wf_type == "TMT - MaxQuant") {
      # Assuming TMT-MQ only has standard/experienced version for now
      if (usr_exp == "experienced") {
         # Note: TMT MQ workflow file name includes version 0.1
        return(paste0(base_url, "/standard/TMT_MQ_workflow0.1.rmd")) 
      } else {
         warning("Beginner TMT-MQ workflow not currently available.")
         return(NULL) # Or point to a default/starter if one exists
      }
      # Add other proteomics workflow_types here (LFQ etc.) when available
    }
    # Fallback warning if specific proteomics workflow not found
    warning("Proteomics workflow URL not found for: type: ", wf_type, ", experience: ", usr_exp)
    return(NULL)

  } else {
    # General Omics: Construct filename and path using omic type and user experience
    experience_suffix <- ifelse(usr_exp == "beginner", "starter", "experienced")
    workflow_filename <- paste0(omic, "_workflow_", experience_suffix, ".rmd")
    workflow_path <- paste0(base_url, "/", subfolder, "/", workflow_filename)
    
    # Optional: Add a check here to see if the constructed URL is valid before returning?
    # For now, we rely on the download step to report errors.
    return(workflow_path)
  }
}

# Get report template URL based on omic type and workflow type
getReportUrl <- function(wf_type, omic) {
  base_url <- paste0("https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/", omic, "/report")
  
  if (omic == "proteomics") {
    # Proteomics: Use wf_type to determine report
    if (wf_type == "DIA-NN") {
      return(paste0(base_url, "/DIANN_report.rmd"))
    } else if (wf_type == "TMT - MaxQuant" || wf_type == "TMT - FragPipe") {
      return(paste0(base_url, "/TMT_report.rmd"))
    } else if (wf_type == "LFQ - MaxQuant" || wf_type == "LFQ - FragPipe") {
      return(paste0(base_url, "/LFQ_report.rmd"))
    }
    warning("Proteomics report URL not found for: type: ", wf_type)
    return(NULL)
  } else {
    # Other Omics: Assume a standard report name, ignore wf_type
    standard_report_path <- paste0(base_url, "/", omic, "_report.rmd")
    return(standard_report_path)
  }
}

# --- Main Project Setup Function ---
setupOmicsProject <- function(root_dir = NULL, overwrite = FALSE, omic_types, workflow_type, user_experience) {
  # Set default root_dir based on OS
  if (is.null(root_dir)) {
    if (.Platform$OS.type == "windows") {
      docs_path <- file.path(Sys.getenv("USERPROFILE"), "Documents")
    } else if (Sys.info()["sysname"] == "Darwin") {  # Mac OS
      docs_path <- file.path(path.expand("~"), "Documents")
    } else {  # Linux or other Unix-like systems
      docs_path <- path.expand("~")
    }
    root_dir <- file.path(docs_path, "default_project")
    
    os_type <- if (.Platform$OS.type == "windows") {
      "Windows"
    } else if (Sys.info()["sysname"] == "Darwin") {
      # Get Mac OS version
      os_version <- system("sw_vers -productVersion", intern = TRUE)
      mac_name <- if (grepl("^14", os_version)) {
        "Sonoma (14)"
      } else if (grepl("^13", os_version)) {
        "Ventura (13)"
      } else if (grepl("^12", os_version)) {
        "Monterey (12)"
      } else if (grepl("^11", os_version)) {
        "Big Sur (11)"
      } else if (grepl("^10.15", os_version)) {
        "Catalina (10.15)"
      } else if (grepl("^10.14", os_version)) {
        "Mojave (10.14)"
      } else if (grepl("^10.13", os_version)) {
        "High Sierra (10.13)"
      } else if (grepl("^10.12", os_version)) {
        "Sierra (10.12)"
      } else if (grepl("^10.11", os_version)) {
        "El Capitan (10.11)"
      } else if (grepl("^10.10", os_version)) {
        "Yosemite (10.10)"
      } else {
        paste("MacOS", os_version)
      }
      message("MacOS version detected: ", mac_name)
      "MacOS"
    } else {
      "Unix-like system"
    }
    message("Operating system detected: ", os_type)
    message("Using default project location: ", root_dir)
  }
  
  # Create base directory structure (once)
  base_dirs <- list(
    root = root_dir,
    scripts = file.path(root_dir, "scripts"),
    data = file.path(root_dir, "data")
  )
  for (dir in base_dirs) {
    if (!dir.exists(dir)) {
      dir.create(dir, recursive = TRUE)
      message("Created directory: ", dir)
    }
  }
  
  # List to store paths of downloaded workflow/config files (primarily for return value)
  downloaded_files <- list()
  
  # --- Download central config.ini file --- 
  config_url <- "https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/config.ini"
  config_dest <- file.path(root_dir, "config.ini")
  message("\nAttempting to download central config.ini...")
  if (!overwrite && file.exists(config_dest)) {
    message("Skipping existing file: config.ini")
    downloaded_files$config <- normalizePath(config_dest)
  } else {
    response_config <- httr::GET(config_url)
    if (httr::status_code(response_config) == 200) {
      content_config <- httr::content(response_config, "raw")
      writeBin(content_config, config_dest)
      message("Successfully downloaded: config.ini to ", root_dir)
      downloaded_files$config <- normalizePath(config_dest)
    } else {
      warning("Failed to download central config.ini. Status code: ", 
              httr::status_code(response_config), ". URL: ", config_url)
      downloaded_files$config <- NULL # Indicate failure
    }
  }
  
  # --- Loop through each specified omic type ---
  for (current_omic_type in omic_types) {
    message(paste0("\n--- Setting up for: ", current_omic_type, " ---"))
    
    # Create omic-specific directories
    omic_dirs <- list(
      scripts_omic = file.path(base_dirs$scripts, current_omic_type),
      data_omic = file.path(base_dirs$data, current_omic_type)
    )
    for (dir in omic_dirs) {
      if (!dir.exists(dir)) {
        dir.create(dir, recursive = TRUE)
        message("Created directory: ", dir)
      }
    }
    
    # Special case for UniProt directory (if proteomics is selected)
    if (current_omic_type == "proteomics") {
      uniprot_dir <- file.path(base_dirs$data, "UniProt")
      if (!dir.exists(uniprot_dir)) {
        dir.create(uniprot_dir, recursive = TRUE)
        message("Created directory: ", uniprot_dir)
      }
    }
    
    # Get the appropriate workflow URL for the current omic (using global function)
    workflow_url <- getWorkflowUrl(workflow_type, user_experience, current_omic_type)
    
    # Get the appropriate report URL for the current omic (using global function)
    report_url <- getReportUrl(workflow_type, current_omic_type)
    
    # Check if workflow URL was found
    if (is.null(workflow_url)) {
      warning("Skipping workflow download for ", current_omic_type, " due to missing URL.")
    } else {
      # Determine workflow filename (potentially omic-dependent logic needed here too)
      workflow_filename <- basename(workflow_url) # Simple approach: get filename from URL
      
      # Determine report filename (potentially omic-dependent logic needed here too)
      report_filename <- if (!is.null(report_url)) basename(report_url) else NULL
      
      # Define GitHub raw content URLs for the current omic
      templates <- list() # Initialize empty list for this omic
      templates$workflow <- list(
        url = workflow_url,
        dest = file.path(omic_dirs$scripts_omic, workflow_filename)
      )
      # Always add config (assuming one config per omic type folder)
      # config_url <- paste0("https://raw.githubusercontent.com/APAF-bioinformatics/MultiScholaR/main/Workbooks/", current_omic_type, "/config.ini")
      # templates$config <- list(
      #  url = config_url,
      #  dest = file.path(omic_dirs$scripts_omic, "config.ini")
      # )
      
      # Add report template if URL exists
      if (!is.null(report_url) && !is.null(report_filename)) {
        templates$report <- list(
          url = report_url,
          dest = file.path(omic_dirs$scripts_omic, report_filename)
        )
      }
      
      # Check if files exist for the current omic
      dest_files_exist <- file.exists(sapply(templates, `[[`, "dest"))
      if (!overwrite && any(dest_files_exist)) {
        warning("Destination files already exist for ", current_omic_type, ". Skipping download for existing files (use overwrite=TRUE to replace).")
        # Optionally skip download loop entirely for this omic if any exist? Current logic skips individual files.
      }
      
      # Download files for the current omic
      message("\nDownloading template files for ", current_omic_type, "...")
      for (template_name in names(templates)) {
        template <- templates[[template_name]]
        # Skip download if overwrite is FALSE and file exists
        if (!overwrite && file.exists(template$dest)) {
          message("Skipping existing file: ", basename(template$dest))
          next
        }
        response <- httr::GET(template$url)
        if (httr::status_code(response) == 200) {
          content <- httr::content(response, "raw")
          # Ensure target directory exists before writing
          dir.create(dirname(template$dest), showWarnings = FALSE, recursive = TRUE)
          writeBin(content, template$dest)
          message("Successfully downloaded: ", basename(template$dest))
          # Store downloaded file path
          if (template_name == "workflow") downloaded_files$workflow[[current_omic_type]] <- normalizePath(template$dest)
          # if (template_name == "config") downloaded_files$config[[current_omic_type]] <- normalizePath(template$dest)
          if (template_name == "report") downloaded_files$report[[current_omic_type]] <- normalizePath(template$dest)
        } else {
          warning("Failed to download ", basename(template$dest), " for ", current_omic_type,
                  ". Status code: ", httr::status_code(response), ". URL: ", template$url)
        }
      }
      message("Finished downloads for ", current_omic_type, ".")
    }
  } # --- End of loop through omic types ---
  
  
  message("\nProject setup complete!")
  message("Project root: ", normalizePath(root_dir))
  # Print paths for the *first* omic type's files if they exist
  first_omic <- omic_types[1]
  if (!is.null(downloaded_files$workflow[[first_omic]])) {
    message("Workflow file (first omic): ", downloaded_files$workflow[[first_omic]])
  }
  if (!is.null(downloaded_files$config)) { # Check the central config download
    message("Config file: ", downloaded_files$config)
  }
  # You could add a loop here to print all downloaded files if needed
  
  # Return value (can be adjusted based on what's most useful)
  return(invisible(list(
    root_dir = normalizePath(root_dir),
    directories = lapply(base_dirs, normalizePath), # Return base dirs
    downloaded_files = downloaded_files # Return nested list of downloaded files by type/omic
  )))
}

# Determine project path based on user input
if (!is.null(my_project_dir)) {
  project_path <- file.path(my_project_dir, my_project_name)
} else {
  if (.Platform$OS.type == "windows") {
    project_path <- file.path(Sys.getenv("USERPROFILE"), "Documents", my_project_name)
  } else {
    project_path <- file.path(path.expand("~"), "Documents", my_project_name)
  }
}

# Create and setup the project
message("Creating project: ", project_path)
setup_result <- setupOmicsProject(project_path, overwrite = TRUE, omic_types, workflow_type, user_experience)

# Create and open R project
rproj_content <- c(
  "Version: 1.0",
  "",
  "RestoreWorkspace: Default",
  "SaveWorkspace: Default",
  "AlwaysSaveHistory: Default",
  "",
  "EnableCodeIndexing: Yes",
  "UseSpacesForTab: Yes",
  "NumSpacesForTab: 2",
  "Encoding: UTF-8",
  "",
  "RnwWeave: Sweave",
  "LaTeX: pdfLaTeX"
)
rproj_file <- file.path(project_path, paste0(my_project_name, ".Rproj"))
writeLines(rproj_content, rproj_file)
message("Created R project file: ", rproj_file)

# Create .Rprofile for automatic workflow opening
startup_script <- file.path(project_path, ".Rprofile")

# Generate .Rprofile content to open ALL downloaded workflows
startup_code_lines <- c(
  'if (interactive()) {',
  '  message("Initializing project...")',
  '  if (!requireNamespace("later", quietly = TRUE)) install.packages("later")',
  '  if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")',
  '  later::later(function() {',
  '    Sys.sleep(2) # Brief pause'
)

# Add commands to open each downloaded workflow
opened_workflows_message <- c()
if (!is.null(setup_result$downloaded_files$workflow) && length(setup_result$downloaded_files$workflow) > 0) {
  for (omic in names(setup_result$downloaded_files$workflow)) {
    workflow_path_full <- setup_result$downloaded_files$workflow[[omic]]
    if (!is.null(workflow_path_full) && file.exists(workflow_path_full)) {
      # Construct the relative path for use within the project
      workflow_filename <- basename(workflow_path_full)
      relative_workflow_path <- file.path("scripts", omic, workflow_filename)
      # Escape backslashes for Windows paths within the R string
      relative_workflow_path_escaped <- gsub("\\\\", "/", relative_workflow_path) 
      
      startup_code_lines <- c(
        startup_code_lines,
        paste0('    workflow_path_to_open <- "', relative_workflow_path_escaped, '"'),
        '    if (file.exists(workflow_path_to_open) && rstudioapi::isAvailable()) {',
        '      message("Attempting to open workflow: ", workflow_path_to_open)',
        '      try(rstudioapi::navigateToFile(workflow_path_to_open))',
        '    } else {',
        '      message("Could not automatically open workflow. Path checked: ", workflow_path_to_open)',
        '    }'
      )
      opened_workflows_message <- c(opened_workflows_message, relative_workflow_path_escaped)
    } else {
        startup_code_lines <- c(
            startup_code_lines,
            paste0('    message("Skipping auto-open for ', omic, ' workflow: File not found or path is NULL.")')
        )
    }
  }
} else {
    startup_code_lines <- c(
        startup_code_lines,
        '    message("No downloaded workflow files found to configure for auto-opening.")'
    )
}

startup_code_lines <- c(startup_code_lines, '  }, 3) # End later', '} # End interactive')
startup_content <- paste(startup_code_lines, collapse = '\n')

writeLines(startup_content, startup_script)
if (length(opened_workflows_message) > 0) {
    message("Created startup script (.Rprofile) to attempt opening the following workflows: ", paste(opened_workflows_message, collapse=", "))
} else {
    message("Created startup script (.Rprofile), but no workflows were found to configure for auto-opening.")
}


# Open the new project
if (rstudioapi::isAvailable()) {
  message("Opening new R project...")
  message("Note: If you have unsaved changes, you\'ll be prompted to save them")
  if (length(opened_workflows_message) > 0) {
      message("The following workflows will attempt to open automatically: ", paste(opened_workflows_message, collapse=", "))
  } else {
      message("No workflows configured to open automatically.")
  }
  Sys.sleep(2)  # Give user time to read messages
  rstudioapi::openProject(rproj_file)
}

message("\nSetup complete! Your new project is ready to use.")