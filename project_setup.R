# ============================================================================
# Proteomics Workflow Project Setup Script
# ============================================================================
# Hi there! Welcome to ProteomeScholaR :)
# Instructions: 
# 1. Set your project name and optional directory below
# 2. Select ALL of this code (Ctrl+A or Cmd+A)
# 3. Run it (Ctrl+Enter or Cmd+Enter)
# ============================================================================

# SET YOUR PROJECT OPTIONS HERE:
my_project_name <- "testing1234"  # Required: Name your project
my_project_dir <- NULL            # Optional: Set specific directory (leave as NULL for default Documents folder)


# Examples for my_project_dir:
# Windows: my_project_dir <- "C:/Users/username/Projects"
# Mac:     my_project_dir <- "~/Projects"

# SET YOUR WORKFLOW OPTIONS HERE:
omic_types <- c("proteomics", "metabolomics") # Options: c("proteomics", "metabolomics", "transcriptomics", ...) Add more as they become available

# The following two options primarily apply when 'proteomics' is in omic_types:
workflow_type <- "DIA-NN"  # Proteomics-specific: "DIA-NN", "LFQ - FragPipe", "LFQ - MaxQuant", "TMT - MaxQuant", "TMT - FragPipe"
user_experience <- "experienced"  # Proteomics-specific: "experienced", "beginner"

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
  base_url <- paste0("https://raw.githubusercontent.com/APAF-bioinformatics/ProteomeScholaR/main/Workbooks/", omic)
  
  if (omic == "proteomics") {
    # Proteomics: Use wf_type and usr_exp
    if (wf_type == "DIA-NN") {
      if (usr_exp == "experienced") {
        return(paste0(base_url, "/standard/DIA_workflow_experienced.rmd"))
      } else if (usr_exp == "beginner") {
        return(paste0(base_url, "/starter/DIA_workflow_starter.rmd"))
      }
    } else if (wf_type == "TMT - MaxQuant") {
      if (usr_exp == "experienced") {
        return(paste0(base_url, "/standard/TMT_MQ_workflow0.1.rmd"))
      }
      # Add beginner TMT-MQ if it exists
    } # Add other proteomics workflow_types here (LFQ etc.)
    
    warning("Proteomics workflow URL not found for: type: ", wf_type, ", experience: ", usr_exp)
    return(NULL)
    
  } else {
    # Other Omics: Assume a standard workflow name, ignore wf_type and usr_exp
    # Example: Assume workflow is always named '[omic]_workflow.rmd' in the 'standard' folder
    standard_workflow_path <- paste0(base_url, "/standard/", omic, "_workflow.rmd")
    return(standard_workflow_path)
  }
}

# Get report template URL based on omic type and workflow type
getReportUrl <- function(wf_type, omic) {
  base_url <- paste0("https://raw.githubusercontent.com/APAF-bioinformatics/ProteomeScholaR/main/Workbooks/", omic, "/report")
  
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
      config_url <- paste0("https://raw.githubusercontent.com/APAF-bioinformatics/ProteomeScholaR/main/Workbooks/", current_omic_type, "/config.ini")
      templates$config <- list(
        url = config_url,
        dest = file.path(omic_dirs$scripts_omic, "config.ini")
      )
      
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
          if (template_name == "config") downloaded_files$config[[current_omic_type]] <- normalizePath(template$dest)
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
  if (!is.null(downloaded_files$config[[first_omic]])) {
    message("Config file (first omic): ", downloaded_files$config[[first_omic]])
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

# Create .Rprofile for automatic workflow opening (opens the *first* specified omic workflow)
startup_script <- file.path(project_path, ".Rprofile")

# Determine workflow filename for the *first* omic type specified for .Rprofile
first_omic_type <- omic_types[1]
# NOW this call will work because getWorkflowUrl is global:
first_workflow_url_for_profile <- getWorkflowUrl(workflow_type, user_experience, first_omic_type)
first_workflow_filename_for_profile <- if (!is.null(first_workflow_url_for_profile)) basename(first_workflow_url_for_profile) else "UNKNOWN_WORKFLOW.rmd"

# Generate .Rprofile content
startup_content <- paste0(
  'if (interactive()) {\n',
  '  message("Initializing project...")\n',
  '  if (!requireNamespace("later", quietly = TRUE)) install.packages("later")\n',
  '  if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")\n',
  '  later::later(function() {\n',
  '    Sys.sleep(2)\n',
  # Use the first omic type and its derived workflow filename
  '    workflow_path <- file.path("scripts", "', first_omic_type, '", "', first_workflow_filename_for_profile, '")\n',
  '    if (file.exists(workflow_path) && rstudioapi::isAvailable()) {\n',
  '      message("Attempting to open workflow: ", workflow_path)\n',
  '      try(rstudioapi::navigateToFile(workflow_path))\n',
  '    } else {\n',
  '      message("Could not automatically open workflow. Please open manually. Path checked: ", workflow_path)\n',
  '    }\\n',
  '  }, 3)\\n',
  '}\\n'
)

writeLines(startup_content, startup_script)
message("Created startup script for automatic workflow opening (targets first omic type: ", first_omic_type, ")")

# Open the new project
if (rstudioapi::isAvailable()) {
  message("Opening new R project...")
  message("Note: If you have unsaved changes, you\'ll be prompted to save them")
  # Adjust message to reflect which workflow will attempt to open
  message(paste0("The workflow file for the first omic type (", first_omic_type, ") will attempt to open automatically at scripts/", first_omic_type, "/", first_workflow_filename_for_profile))
  Sys.sleep(2)  # Give user time to read messages
  rstudioapi::openProject(rproj_file)
}

message("\nSetup complete! Your new project is ready to use.")