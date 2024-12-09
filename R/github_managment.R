#' Push Standard Project Files to GitHub
#'
#' @param base_dir Base directory of the project
#' @param source_dir Source directory containing workflow files
#' @param project_id Project identifier for the new repository and workflow file
#' @param github_org Your GitHub organization or username
#' @param github_user_email Email associated with GitHub account
#' @param github_user_name GitHub username for commits
#' @param commit_message Commit message (default: "Initial project commit")
#' @param private Boolean indicating if repository should be private (default: TRUE)
#'
#' @return Invisible TRUE if successful
#' @export
pushProjectToGithub <- function(base_dir,
                               source_dir,
                               project_id,
                               github_org = getOption("github_org", "your_organization_here"),
                               github_user_email = getOption("github_user_email", "your_email@example.com"),
                               github_user_name = getOption("github_user_name", "your_github_username"),
                               commit_message = "Initial project commit",
                               private = TRUE) {
  
  # Load required packages
  requireNamespace("gh", quietly = TRUE)
  requireNamespace("git2r", quietly = TRUE)
  
  # Use gh package's token function
  github_token <- gh::gh_token()
  
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
  
  # Define standard files to copy
  workflow_file <- file.path(base_dir, paste0(project_id, ".rmd"))
  rproj_file <- file.path(base_dir, paste0(project_id, ".Rproj"))
  
  # Verify files exist
  if (!dir.exists(source_dir)) {
    stop("source_dir does not exist: ", source_dir)
  }
  if (!file.exists(workflow_file)) {
    stop("Workflow file does not exist: ", workflow_file)
  }
  if (!file.exists(rproj_file)) {
    stop("R project file does not exist: ", rproj_file)
  }
  
  # Copy source directory with proper recursive copying
  source_dir_name <- "study_parameters" 
  dir.create(file.path(temp_dir, source_dir_name), recursive = TRUE)
  file.copy(
    from = list.files(source_dir, full.names = TRUE),
    to = file.path(temp_dir, source_dir_name),
    recursive = TRUE,
    overwrite = TRUE
  )
  
  # Copy workflow file maintaining directory structure
  workflow_dir <- file.path(temp_dir, "Workbooks")
  dir.create(workflow_dir, recursive = TRUE)
  file.copy(
    from = workflow_file,
    to = file.path(workflow_dir, basename(workflow_file)),
    overwrite = TRUE
  )
  
  # Copy Rproj file
  file.copy(
    from = rproj_file,
    to = file.path(temp_dir, basename(rproj_file)),
    overwrite = TRUE
  )
  
    # Initialize git repository
  repo <- git2r::init(temp_dir)
  
  # Configure git user
  git2r::config(repo, user.name = github_user_name, user.email = github_user_email)
  
  # Stage all files
  git2r::add(repo, "*")
  
  # Create initial commit
  git2r::commit(repo, message = commit_message)
  
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
  
  # Get current branch name
  current_branch <- git2r::repository_head(repo)$name
  
  # Push to GitHub using the current branch name
  tryCatch({
    git2r::push(repo, "origin", paste0("refs/heads/", current_branch))
  }, error = function(e) {
    stop("Failed to push to repository: ", e$message)
  })
  
  # Clean up
  unlink(temp_dir, recursive = TRUE)
  
  # Create success message with clickable URL
  repo_url <- sprintf("https://github.com/%s/%s", github_org, repo_name)
  message("Successfully pushed files to ", repo_url)
  message("Repository created at: ", repo_url)
  
  invisible(TRUE)
}

### Example Usage
###  options(
###  github_org = "your org",
###  github_user_email = "your email",
###  github_user_name = "your username"
### )

#Example usage:
### pushProjectToGithub(
### base_dir = base_dir,
### source_dir = source_dir,
### project_id = "your project"
### )
