#' Push Standard Project Files to GitHub
#'
#' @param base_dir Base directory of the project
#' @param source_dir Source directory containing workflow files (scripts/proteomics)
#' @param project_id Project identifier for the new repository and workflow file
#' @param github_org Your GitHub organization or username
#' @param github_user_email Email associated with GitHub account
#' @param github_user_name GitHub username for commits
#' @param commit_message Commit message (default: "Initial project commit")
#' @param private Boolean indicating if repository should be private (default: TRUE)
#'
#' @return Invisible TRUE if successful
#' @importFrom utils menu
#' @importFrom fs path_rel dir_ls dir_tree
#' @export
pushProjectToGithub <- function(base_dir,
                               source_dir,
                               project_id,
                               github_org = getOption("github_org", "your_organization_here"),
                               github_user_email = getOption("github_user_email", "your_email@example.com"),
                               github_user_name = getOption("github_user_name", "your_github_username"),
                               commit_message = "Initial project commit",
                               private = TRUE) {
  
  # === DEBUG: Print base_dir ===
  # message("DEBUG: Entered pushProjectToGithub. base_dir: ", base_dir)
  
  # Use gh package's token function
  github_token <- gh::gh_token()
  
  # === DEBUG: Check token ===
  if (is.null(github_token) || nchar(github_token) == 0) {
      stop("DEBUG: Failed to retrieve GitHub token from gh::gh_token()")
  } else {
      # message("DEBUG: Retrieved GitHub token (starts with: ", substr(github_token, 1, 6), "...)")
  }
  
  # Validate inputs
  if (missing(project_id)) {
    stop("project_id is required")
  }
  
  if (github_org == "your_organization_here") {
    stop("Please specify your GitHub organization or username using options() or as a parameter")
  }
  
  # Create new repository name
  repo_name <- project_id
  
  # Create temporary directory for files
  temp_dir <- file.path(tempdir(), repo_name)
  
  # Clean up existing temp directory if it exists
  if (dir.exists(temp_dir)) {
    unlink(temp_dir, recursive = TRUE, force = TRUE)
  }
  dir.create(temp_dir, recursive = TRUE)
  
  # === DEBUG: Print temp_dir ===
  # message("DEBUG: Created temp_dir: ", temp_dir)
  
  # Recursively find all R Markdown files in base_dir (case-insensitive) using fs
  # This tends to be more robust with paths and regex across platforms
  # STEP 1: List ALL files recursively
  all_files_found <- fs::dir_ls(
      path = base_dir
      # , regexp = "\\.[Rr][Mm][Dd]$" # Remove regex for now
      , recurse = TRUE
      , type = "file"
      , fail = FALSE # Don't error if path doesn't exist, just return empty
  )
  
  # === DEBUG: Print ALL files found by fs::dir_ls before ANY filtering ===
  # message("DEBUG: ALL files found by fs::dir_ls (no pattern):")
  # print(all_files_found)
  
  # STEP 2: Filter for Rmd files using grepl (case-insensitive)
  all_rmd_files <- all_files_found[grepl("\\.[Rr][Mm][Dd]$", all_files_found)]
  
  # === DEBUG: Print files after grepl filtering ===
  # message("DEBUG: Files after grepl filtering for Rmd:")
  print(all_rmd_files)
  
  # === DEBUG: Print files found by fs::dir_ls BEFORE filtering ===
  # message("DEBUG: Files found by fs::dir_ls before filtering:") # Original debug line
  # print(all_rmd_files) # Original debug line
  
  workflow_file <- NULL # Initialize variable
  
  if (length(all_rmd_files) == 0) {
    stop("No R Markdown files (.Rmd, .rmd, .RMD) found recursively in: ", base_dir)
  } else if (length(all_rmd_files) == 1) {
    workflow_file <- all_rmd_files[1]
    message("Using the only R Markdown file found: ", fs::path_rel(workflow_file, start = base_dir))
  } else {
    message("Multiple R Markdown files found. Please select one to use as the main project workflow:")
    # Display relative paths for clarity in the menu
    relative_paths <- fs::path_rel(all_rmd_files, start = base_dir)
    selection <- utils::menu(
        choices = relative_paths
        , title = "Select the main workflow file:"
    )
    if (selection == 0) {
        stop("User cancelled file selection.")
    }
    workflow_file <- all_rmd_files[selection]
    # Convert selected fs_path back to character for downstream use if needed, though file.copy handles it
    workflow_file_char <- as.character(workflow_file) 
    message("Selected workflow file: ", fs::path_rel(workflow_file_char, start = base_dir))
  }
  
  # Define the new workflow filename with project_id
  new_workflow_name <- paste0(project_id, ".Rmd") # Use consistent .Rmd extension
  
  # Create directory structure in temp folder
  scripts_proteomics_dir <- file.path(temp_dir, "scripts", "proteomics")
  dir.create(scripts_proteomics_dir, recursive = TRUE)
  
  # === Create .gitignore file ===
  gitignore_path <- file.path(temp_dir, ".gitignore")
  gitignore_content <- c(
    "# Ignore large data files",
    "scripts/proteomics/data_cln.tab",
    "",
    "# R temporary files and caches",
    ".Rproj.user/",
    ".Rhistory",
    ".Rdata",
    ".Ruserdata",
    "*.RData",
    ".httr-oauth",
    "*\\.Rproj.user",
    "*\\.Rhistory",
    "*\\.Rdata",
    "*\\.Ruserdata",
    "*\\*.RData",
    "*\\.httr-oauth"
  )
  writeLines(gitignore_content, gitignore_path)
  # message("DEBUG: Created .gitignore file at: ", gitignore_path)
  
  # Copy and rename workflow file to scripts/proteomics with new name
  file.copy(
    from = workflow_file,
    to = file.path(scripts_proteomics_dir, new_workflow_name),
    overwrite = TRUE
  )
  
  # Copy Rproj file
  rproj_file <- file.path(base_dir, paste0(project_id, ".Rproj"))
  if (!file.exists(rproj_file)) {
    stop("R project file does not exist: ", rproj_file)
  }
  file.copy(
    from = rproj_file,
    to = file.path(temp_dir, basename(rproj_file)),
    overwrite = TRUE
  )
  
  # Copy other files from source_dir (excluding any DIA workflow files)
  source_files <- list.files(source_dir, full.names = TRUE)
  # Ensure we don't accidentally copy the selected workflow file again if it happens to be in source_dir
  # Also keep the original exclusion logic for generic DIA_workflow files
  # *** AND explicitly exclude data_cln.tab ***
  normalized_workflow_file <- normalizePath(workflow_file_char, winslash = "/")
  excluded_files <- c("DIA_workflow\\.R[mM][dD]$", "data_cln\\.tab$") # Add data_cln.tab pattern (fixed escape)
  source_files_to_copy <- source_files[
    !grepl(paste(excluded_files, collapse="|"), basename(source_files)) & # Check basename against patterns
    normalizePath(source_files, winslash = "/") != normalized_workflow_file
  ]

  if (length(source_files_to_copy) > 0) {
    file.copy(
      from = source_files_to_copy,
      to = scripts_proteomics_dir,
      recursive = TRUE,
      overwrite = TRUE
    )
  }
  
  # Initialize git repository
  repo <- git2r::init(temp_dir)
  
  # === DEBUG: List files in temp_dir BEFORE git add ===
  # message("DEBUG: Contents of temp_dir before git add:")
  # print(fs::dir_tree(temp_dir))
  
  # Configure git user
  git2r::config(repo, user.name = github_user_name, user.email = github_user_email)
  
  # Stage all files
  git2r::add(repo, "*")
  
  # === DEBUG: Check git status AFTER git add ===
  # message("DEBUG: Git status after git add:")
  # print(git2r::status(repo))
  
  # Create initial commit
  git2r::commit(repo, message = commit_message)
  
  # === DEBUG: Show commit details ===
  # This requires the commit object, let's get the last commit
  last_commit <- git2r::last_commit(repo)
  # Get the current branch AFTER the commit
  current_branch_obj <- git2r::repository_head(repo)
  # Check if it's a valid branch head (not detached)
  if (!git2r::is_branch(current_branch_obj)) {
      stop("Repository is in a detached HEAD state after commit.")
  }
  current_branch_name <- current_branch_obj$name
  # message("DEBUG: Branch name after commit: ", current_branch_name)
  
  # Rename branch to 'main' if it's not already 'main'
  target_branch_name <- "main"
  if (current_branch_name != target_branch_name) {
    # message("DEBUG: Renaming local branch from '", current_branch_name, "' to '", target_branch_name, "'")
    # Find the branch object and rename it
    branch_to_rename <- git2r::branches(repo, "local")[[current_branch_name]]
    if (is.null(branch_to_rename)) {
        stop("Could not find local branch object: ", current_branch_name)
    }
    git2r::branch_rename(branch_to_rename, target_branch_name)
    # message("DEBUG: Branch renamed successfully.")
    # Update current branch name variable for push
    current_branch_name <- target_branch_name 
  }
  
  # message("DEBUG: Last commit details:")
  # print(last_commit)
  
  # Create repository in the organization
  tryCatch({
    gh::gh(
      "POST /orgs/{org}/repos",
      org = github_org,
      .token = github_token,
      name = repo_name,
      private = TRUE,
      auto_init = FALSE,
      description = sprintf("Project repository for %s", project_id)
    )
  }, error = function(e) {
    stop(sprintf("Failed to create repository: %s", e$message))
  })
  
  # Add remote
  remote_url <- sprintf("https://x-access-token:%s@github.com/%s/%s.git",
                       github_token, github_org, repo_name)
  git2r::remote_add(repo, "origin", remote_url)
  # message("DEBUG: Added remote 'origin' with URL: ", sprintf("https://github.com/%s/%s.git", github_org, repo_name)) # Hide token in debug
  
  # Push to GitHub using the potentially renamed 'main' branch via system() call
  tryCatch({
    # Construct the git command to run in the temp directory
    git_cmd <- sprintf("git -C %s push origin %s", 
                       shQuote(normalizePath(temp_dir, winslash = "/", mustWork = FALSE)), # Normalize and quote path
                       paste0("refs/heads/", current_branch_name))
                       
    # message("DEBUG: Attempting push via system() command: git push origin ", paste0("refs/heads/", current_branch_name), " (in dir: ", normalizePath(temp_dir, winslash = "/", mustWork = FALSE), ")")
    
    # Execute the command
    push_output <- system(git_cmd, intern = TRUE, ignore.stderr = FALSE)
    
    # Check the system command exit status (0 means success)
    exit_status <- attr(push_output, "status")
    # message("DEBUG: system('git push') exit status: ", exit_status)
    # message("DEBUG: system('git push') output:")
    # print(push_output)
    
    if (!is.null(exit_status) && exit_status != 0) {
        stop("system('git push') failed with exit code: ", exit_status)
    }
    
    # === DEBUG: Push successful ===
    # message("DEBUG: system('git push') command appeared successful (exit code 0).")
  }, error = function(e) {
    # === DEBUG: Push FAILED ===
    message("DEBUG: Push operation failed. Error:")
    print(e)
    # Include original error message if available
    stop("Failed to push to repository. Error during push step: ", e$message)
  })
  
  # Clean up
  # === DEBUG: Temporarily disabled cleanup to allow inspection ===
  # message("DEBUG: Skipping cleanup of temp_dir: ", temp_dir)
  unlink(temp_dir, recursive = TRUE)
  
  # Create success message with clickable URL
  repo_url <- sprintf("https://github.com/%s/%s", github_org, repo_name)
  message("Successfully pushed files to ", repo_url)
  message("Repository created at: ", repo_url)
  message("Workflow file renamed to: ", new_workflow_name)
  
  invisible(TRUE)
}