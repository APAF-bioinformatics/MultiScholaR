# ============================================================================
# DIA Workflow Project Setup Script
# ============================================================================
# Hi there! Welcome to ProteomeScholaR :)
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

# Define the setup function
setup_dia_project <- function(root_dir = NULL, overwrite = FALSE) {
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
    
    # Create directory structure
    dirs <- list(
        root = root_dir,
        scripts = file.path(root_dir, "scripts"),
        scripts_proteomics = file.path(root_dir, "scripts", "proteomics"),
        data = file.path(root_dir, "data"),
        data_proteomics = file.path(root_dir, "data", "proteomics"),
        data_uniprot = file.path(root_dir, "data", "UniProt")
    )
    
    # Create directories
    for (dir in dirs) {
        if (!dir.exists(dir)) {
            dir.create(dir, recursive = TRUE)
            message("Created directory: ", dir)
        }
    }
    
    # Define GitHub raw content URLs
    templates <- list(
        workflow = list(
            url = "https://raw.githubusercontent.com/APAF-bioinformatics/ProteomeScholaR/dev-jr/Workbooks/DIA_workflow.rmd",
            dest = file.path(dirs$scripts_proteomics, "DIA_workflow.rmd")
        ),
        config = list(
            url = "https://raw.githubusercontent.com/APAF-bioinformatics/ProteomeScholaR/dev-jr/Workbooks/config.ini",
            dest = file.path(dirs$scripts_proteomics, "config.ini")
        )
    )
    
    # Check if files exist
    if (!overwrite && any(file.exists(sapply(templates, `[[`, "dest")))) {
        stop("Destination files already exist. Use overwrite=TRUE to replace them.")
    }
    
    # Download files
    message("\nDownloading template files from GitHub...")
    for (template in templates) {
        response <- httr::GET(template$url)
        if (httr::status_code(response) == 200) {
            content <- httr::content(response, "raw")  # Get raw content instead of text
            writeBin(content, template$dest)  # Write binary data directly
            message("Successfully downloaded: ", basename(template$dest))
        } else {
            stop("Failed to download ", basename(template$dest), 
                 ". Status code: ", httr::status_code(response))
        }
    }
    
    message("\nProject setup complete!")
    message("Project root: ", normalizePath(root_dir))
    message("Workflow file: ", normalizePath(file.path(dirs$scripts_proteomics, "DIA_workflow.rmd")))
    message("Config file: ", normalizePath(templates$config$dest))
    
    return(invisible(list(
        root_dir = normalizePath(root_dir),
        workflow_file = normalizePath(templates$workflow$dest),
        config_file = normalizePath(templates$config$dest),
        directories = lapply(dirs, normalizePath)
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
setup_result <- setup_dia_project(project_path, overwrite = TRUE)

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
startup_content <- paste0(
    'if (interactive()) {\n',
    '  message("Initializing project...")\n',
    '  if (!requireNamespace("later", quietly = TRUE)) install.packages("later")\n',
    '  if (!requireNamespace("rstudioapi", quietly = TRUE)) install.packages("rstudioapi")\n',
    '  later::later(function() {\n',
    '    Sys.sleep(2)\n',
    '    workflow_path <- file.path("scripts", "proteomics", "DIA_workflow.rmd")\n',
    '    if (file.exists(workflow_path) && rstudioapi::isAvailable()) {\n',
    '      try(rstudioapi::navigateToFile(workflow_path))\n',
    '    } else {\n',
    '      message("Please open the workflow file manually at: ", workflow_path)\n',
    '    }\n',
    '  }, 3)\n',
    '}\n'
)

writeLines(startup_content, startup_script)
message("Created startup script for automatic workflow opening")

# Open the new project
if (rstudioapi::isAvailable()) {
    message("Opening new R project...")
    message("Note: If you have unsaved changes, you'll be prompted to save them")
    message("The workflow file will open automatically in the new project at scripts/proteomics/DIA_workflow.rmd")
    Sys.sleep(2)  # Give user time to read messages
    rstudioapi::openProject(rproj_file)
}

message("\nSetup complete! Your new project is ready to use.")