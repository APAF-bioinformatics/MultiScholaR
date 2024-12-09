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
#' @export
pushProjectToGithub <- function(base_dir,
                               source_dir,
                               project_id,
                               github_org = getOption("github_org", "your_organization_here"),
                               github_user_email = getOption("github_user_email", "your_email@example.com"),
                               github_user_name = getOption("github_user_name", "your_github_username"),
                               commit_message = "Initial project commit",
                               private = TRUE) {
  
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
  
  # Find workflow file with any case
  possible_workflow_files <- c(
    file.path(source_dir, "DIA_workflow.rmd"),
    file.path(source_dir, "DIA_workflow.Rmd"),
    file.path(source_dir, "DIA_workflow.RMD")
  )
  workflow_file <- possible_workflow_files[file.exists(possible_workflow_files)][1]
  
  if (is.na(workflow_file)) {
    stop("Could not find DIA workflow file in: ", source_dir)
  }
  
  # Define the new workflow filename with project_id
  new_workflow_name <- paste0(project_id, ".rmd")
  
  # Create directory structure in temp folder
  scripts_proteomics_dir <- file.path(temp_dir, "scripts", "proteomics")
  dir.create(scripts_proteomics_dir, recursive = TRUE)
  
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
  source_files <- source_files[!grepl("DIA_workflow\\.R[mM][dD]$", source_files)]
  if (length(source_files) > 0) {
    file.copy(
      from = source_files,
      to = scripts_proteomics_dir,
      recursive = TRUE,
      overwrite = TRUE
    )
  }
  
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
  message("Workflow file renamed to: ", new_workflow_name)
  
  invisible(TRUE)
}