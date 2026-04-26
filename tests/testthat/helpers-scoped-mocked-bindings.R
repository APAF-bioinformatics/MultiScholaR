scoped_rebind_bindings <- function(env, bindings, .local_envir = parent.frame()) {
  rebound <- character()

  for (name in names(bindings)) {
    if (!exists(name, envir = env, inherits = FALSE)) {
      next
    }

    old_value <- get(name, envir = env, inherits = FALSE)
    was_locked <- bindingIsLocked(name, env)

    if (was_locked) {
      unlockBinding(name, env)
    }
    assign(name, bindings[[name]], envir = env)
    if (was_locked) {
      lockBinding(name, env)
    }

    restore_name <- name
    restore_env <- env
    restore_value <- old_value
    restore_was_locked <- was_locked
    withr::defer({
      if (
        exists(restore_name, envir = restore_env, inherits = FALSE) &&
          bindingIsLocked(restore_name, restore_env)
      ) {
        unlockBinding(restore_name, restore_env)
      }
      assign(restore_name, restore_value, envir = restore_env)
      if (
        restore_was_locked &&
          exists(restore_name, envir = restore_env, inherits = FALSE)
      ) {
        lockBinding(restore_name, restore_env)
      }
    }, envir = .local_envir)

    rebound <- c(rebound, name)
  }

  rebound
}

scoped_mocked_bindings <- function(..., .package = NULL, .env = parent.frame(), .local_envir = parent.frame()) {
  bindings <- list(...)
  binding_names <- names(bindings)
  found <- stats::setNames(rep(FALSE, length(binding_names)), binding_names)

  if (is.null(.package)) {
    rebound <- scoped_rebind_bindings(.env, bindings, .local_envir = .local_envir)
    found[rebound] <- TRUE
  } else {
    namespace_env <- asNamespace(.package)
    attached_name <- paste0("package:", .package)
    candidate_envs <- list(
      namespace_env,
      parent.env(namespace_env),
      globalenv(),
      namespace_env[[".__S3MethodsTable__."]],
      if (attached_name %in% search()) as.environment(attached_name) else NULL,
      testthat:::the$testing_env
    )

    for (candidate_env in candidate_envs) {
      if (!is.environment(candidate_env)) {
        next
      }
      rebound <- scoped_rebind_bindings(candidate_env, bindings, .local_envir = .local_envir)
      found[rebound] <- TRUE
    }
  }

  if (!all(found)) {
    missing <- names(found)[!found]
    stop(sprintf("Can't find binding for %s", paste(missing, collapse = ", ")), call. = FALSE)
  }

  invisible()
}

capture_binding_state <- function(env, names) {
  stats::setNames(lapply(names, function(name) {
    had_binding <- exists(name, envir = env, inherits = FALSE)
    list(
      had_binding = had_binding,
      value = if (had_binding) get(name, envir = env, inherits = FALSE) else NULL,
      was_locked = had_binding && bindingIsLocked(name, env)
    )
  }), names)
}

restore_binding_state <- function(env, state) {
  for (name in names(state)) {
    entry <- state[[name]]
    if (exists(name, envir = env, inherits = FALSE) && bindingIsLocked(name, env)) {
      unlockBinding(name, env)
    }
    if (isTRUE(entry$had_binding)) {
      assign(name, entry$value, envir = env)
    } else if (exists(name, envir = env, inherits = FALSE)) {
      rm(list = name, envir = env)
    }
    if (isTRUE(entry$was_locked) && exists(name, envir = env, inherits = FALSE)) {
      lockBinding(name, env)
    }
  }
}

register_binding_teardown <- function(env, names, .local_envir = parent.frame()) {
  state <- capture_binding_state(env, names)
  teardown_env <- environment()
  testthat:::teardown(
    restore_binding_state(env, state),
    env = teardown_env
  )
  invisible(NULL)
}
