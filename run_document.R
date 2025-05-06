# Simple script to run devtools::document()
library(devtools)

# Set working directory to the package root
# This line will be run if the script is executed from a different directory
# setwd("put/your/package/path/here/if/needed")

# Run document() and capture any errors
tryCatch({
  message("Starting documentation...")
  document()
  message("Documentation completed successfully.")
}, error = function(e) {
  message("Error during documentation process:")
  message(e$message)
  if(!is.null(e$call)) message("Call: ", deparse(e$call))
}) 