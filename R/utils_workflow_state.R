#' @title WorkflowState
#' @description An R6 class to manage the state of a multi-step workflow.
#'
#' @details This class is designed to track the state of data objects (specifically S4 objects)
#' as they are processed through a series of transformations in a Shiny app. It allows for
#' saving named states, retrieving them, and reverting to previous states, which is essential
#' for creating a non-linear, interactive analysis workflow.
#'
#' It works by storing copies of data objects in a list, acting as a "save game"
#' feature for the analysis.
#'
#' @importFrom R6 R6Class
#' @export
WorkflowState <- R6::R6Class("WorkflowState",
  public = list(
    #' @field states A list to store the different data states. Each element is a
    #' list containing the data object, timestamp, configuration, and other metadata.
    states = list(),
    
    #' @field current_state A character string holding the name of the current active state.
    current_state = NULL,
    
    #' @field state_history A character vector that logs the sequence of saved states.
    state_history = list(),

    #' @field workflow_type A character string indicating the type of workflow (e.g., "LFQ", "TMT").
    workflow_type = "LFQ",

    #' @description
    #' Initialize the WorkflowState object.
    #' Sets up the initial state.
    initialize = function() {
      self$current_state <- "initial"
      self$states[["initial"]] <- list(
        timestamp = Sys.time(),
        data = NULL,
        config = NULL,
        description = "Initial empty state"
      )
      self$state_history <- append(self$state_history, "initial")
    },

    #' @description
    #' Save a new state (a snapshot of a data object).
    #' @param state_name A character string for the name of the state (e.g., "raw_data_s4").
    #' @param s4_data_object The S4 data object to save. A copy will be made.
    #' @param config_object A list containing the configuration parameters used to generate this state.
    #' @param description A character string describing the state.
    saveState = function(state_name, s4_data_object, config_object, description) {
      # Due to R's copy-on-modify semantics for S4, assigning the object
      # to the list creates an independent copy.
      self$states[[state_name]] <- list(
        timestamp = Sys.time(),
        data = s4_data_object,
        config = config_object,
        description = description,
        previous_state = self$current_state,
        s4_class = class(s4_data_object)[1]
      )
      # Add to history only if it's a new state
      if (!state_name %in% self$state_history) {
          self$state_history <- append(self$state_history, state_name)
      }
      self$current_state <- state_name
    },

    #' @description
    #' Retrieve a data object from a specific state.
    #' @param state_name The name of the state to retrieve. If NULL, retrieves the current state.
    #' @return The data object (e.g., an S4 object) stored in that state.
    getState = function(state_name = NULL) {
      if (is.null(state_name)) {
        state_name <- self$current_state
      }
      return(self$states[[state_name]]$data)
    },

    #' @description
    #' Revert the current active state to a previously saved state.
    #' @param state_name The name of the state to revert to.
    #' @return The data object from the state that is now active.
    revertToState = function(state_name) {
      if (state_name %in% names(self$states)) {
        # Remove any states that came after the one we are reverting to
        revert_index <- which(self$state_history == state_name)
        self$state_history <- self$state_history[1:revert_index]
        self$current_state <- state_name
        return(self$states[[state_name]]$data)
      }
      stop("State not found: ", state_name)
    },

    #' @description
    #' Get the list of all saved state names.
    #' @return A character vector of state names.
    getHistory = function() {
      return(unlist(self$state_history))
    },
    
    #' @description
    #' Set the workflow type for the session.
    #' Workflow type determines which set of modules to load.
    #' @param type A character string for the workflow type.
    #'   Proteomics: "LFQ", "TMT", "DIA"
    #'   Metabolomics: "metabolomics_standard"
    #'   Lipidomics: "lipidomics_standard"
    setWorkflowType = function(type) {
      allowed_types <- c(
        # Proteomics workflows
        "LFQ", "TMT", "DIA"
        # Metabolomics workflows
        , "metabolomics_standard"
        # Lipidomics workflows
        , "lipidomics_standard"
      )
      if (type %in% allowed_types) {
        self$workflow_type <- type
      } else {
        stop(paste("Invalid workflow type. Must be one of:", paste(allowed_types, collapse = ", ")))
      }
    }
  )
)

