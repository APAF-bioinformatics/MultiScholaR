resolveWorkflowGitInfo <- function() {
    tryCatch(
        {
            pkg_path <- tryCatch(
                {
                    pkg_ns_path <- system.file(package = "MultiScholaR")
                    if (nzchar(pkg_ns_path) && dir.exists(file.path(dirname(pkg_ns_path), ".git"))) {
                        dirname(pkg_ns_path)
                    } else {
                        dev_paths <- c(
                            getwd(),
                            Sys.getenv("MULTISCHOLAR_PATH", unset = NA),
                            file.path(Sys.getenv("HOME"), "Documents", "MultiScholaR")
                        )
                        found_path <- NA_character_
                        for (p in dev_paths) {
                            if (!is.na(p) && nzchar(p) && dir.exists(file.path(p, ".git"))) {
                                found_path <- p
                                break
                            }
                        }
                        found_path
                    }
                },
                error = function(e) NA_character_
            )

            if (!is.na(pkg_path) && dir.exists(file.path(pkg_path, ".git"))) {
                cat(sprintf("WORKFLOW ARGS: Found local git repo at %s\n", pkg_path))

                branch <- tryCatch(
                    {
                        result <- system2("git",
                            args = c("-C", pkg_path, "branch", "--show-current"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(result) > 0 && nzchar(result[1])) result[1] else NA_character_
                    },
                    error = function(e) NA_character_
                )

                commit_sha <- tryCatch(
                    {
                        result <- system2("git",
                            args = c("-C", pkg_path, "rev-parse", "HEAD"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(result) > 0 && nzchar(result[1])) result[1] else NA_character_
                    },
                    error = function(e) NA_character_
                )

                timestamp <- tryCatch(
                    {
                        result <- system2("git",
                            args = c("-C", pkg_path, "log", "-1", "--format=%cI"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(result) > 0 && nzchar(result[1])) result[1] else NA_character_
                    },
                    error = function(e) NA_character_
                )

                repo_name <- tryCatch(
                    {
                        remote <- system2("git",
                            args = c("-C", pkg_path, "remote", "get-url", "origin"),
                            stdout = TRUE, stderr = FALSE
                        )
                        if (length(remote) > 0 && nzchar(remote[1])) {
                            basename(gsub("\\.git$", "", remote[1]))
                        } else {
                            basename(pkg_path)
                        }
                    },
                    error = function(e) basename(pkg_path)
                )

                list(
                    commit_sha = commit_sha,
                    branch = branch,
                    repo = repo_name,
                    timestamp = timestamp,
                    source = "local_git"
                )
            } else {
                cat("WORKFLOW ARGS: No local git repo found, using package info\n")
                pkg_version <- tryCatch(
                    {
                        as.character(utils::packageVersion("MultiScholaR"))
                    },
                    error = function(e) "unknown"
                )

                list(
                    commit_sha = NA_character_,
                    branch = NA_character_,
                    repo = "MultiScholaR",
                    timestamp = NA_character_,
                    package_version = pkg_version,
                    source = "package_info"
                )
            }
        },
        error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error getting git info: %s\n", e$message))
            list(commit_sha = NA_character_, branch = NA_character_, repo = NA_character_, timestamp = NA_character_)
        }
    )
}

formatWorkflowParameterValue <- function(param_value) {
    tryCatch(
        {
            if (is.null(param_value)) {
                "NULL"
            } else if (is.data.frame(param_value)) {
                sprintf(
                    "[Data frame: %d rows x %d cols - omitted for brevity]",
                    nrow(param_value), ncol(param_value)
                )
            } else if (is.logical(param_value)) {
                if (length(param_value) == 1) {
                    ifelse(param_value, "TRUE", "FALSE")
                } else if (length(param_value) > 50) {
                    true_count <- sum(param_value, na.rm = TRUE)
                    total_count <- length(param_value)
                    sprintf(
                        "logical vector [%d TRUE, %d FALSE out of %d total]",
                        true_count, total_count - true_count, total_count
                    )
                } else {
                    preview <- ifelse(utils::head(param_value, 5), "TRUE", "FALSE")
                    if (length(param_value) > 5) {
                        paste0("c(", paste(preview, collapse = ", "), ", ...)")
                    } else {
                        paste0("c(", paste(preview, collapse = ", "), ")")
                    }
                }
            } else if (is.numeric(param_value)) {
                if (length(param_value) == 1) {
                    as.character(param_value)
                } else if (length(param_value) > 5) {
                    sprintf(
                        "numeric vector [%d values: %s, ...]",
                        length(param_value),
                        paste(as.character(utils::head(param_value, 3)), collapse = ", ")
                    )
                } else {
                    paste0("c(", paste(as.character(param_value), collapse = ", "), ")")
                }
            } else if (is.character(param_value)) {
                if (length(param_value) == 1) {
                    param_value
                } else if (length(param_value) > 5) {
                    sprintf(
                        "character vector [%d values: %s, ...]",
                        length(param_value),
                        paste(shQuote(utils::head(param_value, 3)), collapse = ", ")
                    )
                } else {
                    paste0("c(", paste(shQuote(param_value), collapse = ", "), ")")
                }
            } else {
                paste0("[", class(param_value)[1], " object]")
            }
        },
        error = function(e) {
            "[SERIALIZATION ERROR]"
        }
    )
}

buildWorkflowRuvOptimizationLines <- function(ruv_optimization_result, s4_params) {
    if (isTRUE(ruv_optimization_result$ruv_skipped)) {
        cat("WORKFLOW ARGS: RUV was skipped - writing skip section\n")
        return(c(
            "RUV-III Batch Correction:",
            "-------------------------",
            "* Status: Not Applied",
            "* Reason: User determined RUV was not appropriate due to dataset constraints",
            ""
        ))
    }

    cat("WORKFLOW ARGS: Formatting RUV optimization results\n")
    ruv_lines <- c(
        "Automatic RUV Optimization Results:",
        "-----------------------------------"
    )

    tryCatch(
        {
            best_percentage <- if (!is.null(ruv_optimization_result$best_percentage)) {
                sprintf("%.1f%%", ruv_optimization_result$best_percentage)
            } else {
                "N/A"
            }
            best_k <- if (!is.null(ruv_optimization_result$best_k)) {
                as.character(ruv_optimization_result$best_k)
            } else {
                "N/A"
            }
            separation_score <- if (!is.null(ruv_optimization_result$best_separation_score)) {
                sprintf("%.4f", ruv_optimization_result$best_separation_score)
            } else {
                "N/A"
            }
            composite_score <- if (!is.null(ruv_optimization_result$best_composite_score)) {
                sprintf("%.4f", ruv_optimization_result$best_composite_score)
            } else {
                "N/A"
            }
            control_genes_count <- if (!is.null(ruv_optimization_result$best_control_genes_index)) {
                as.character(sum(ruv_optimization_result$best_control_genes_index, na.rm = TRUE))
            } else {
                "N/A"
            }
            separation_metric <- if (!is.null(ruv_optimization_result$separation_metric_used)) {
                as.character(ruv_optimization_result$separation_metric_used)
            } else {
                "N/A"
            }
            k_penalty_weight <- if (!is.null(ruv_optimization_result$k_penalty_weight)) {
                sprintf("%.1f", ruv_optimization_result$k_penalty_weight)
            } else {
                "N/A"
            }
            adaptive_penalty <- if (!is.null(ruv_optimization_result$adaptive_k_penalty_used)) {
                ifelse(ruv_optimization_result$adaptive_k_penalty_used, "TRUE", "FALSE")
            } else {
                "N/A"
            }
            sample_size <- if (!is.null(ruv_optimization_result$sample_size)) {
                as.character(ruv_optimization_result$sample_size)
            } else {
                "N/A"
            }

            ruv_grouping_variable <- "N/A"
            if (!is.null(s4_params$ruvIII_C_Varying$ruv_grouping_variable)) {
                ruv_grouping_variable <- s4_params$ruvIII_C_Varying$ruv_grouping_variable
            } else if (!is.null(s4_params$getNegCtrlProtAnova$ruv_grouping_variable)) {
                ruv_grouping_variable <- s4_params$getNegCtrlProtAnova$ruv_grouping_variable
            }

            ruv_lines <- c(
                ruv_lines,
                paste("* Best percentage:", best_percentage),
                paste("* Best k value:", best_k),
                paste("* Separation score:", separation_score),
                paste("* Composite score:", composite_score),
                paste("* Control genes:", control_genes_count),
                paste("* RUV grouping variable:", ruv_grouping_variable),
                paste("* Separation metric:", separation_metric),
                paste("* K penalty weight:", k_penalty_weight),
                paste("* Adaptive penalty:", adaptive_penalty),
                paste("* Sample size:", sample_size),
                ""
            )
            cat("WORKFLOW ARGS: Successfully formatted RUV optimization results\n")
        },
        error = function(e) {
            cat(sprintf("WORKFLOW ARGS: Error formatting RUV results: %s\n", e$message))
            ruv_lines <- c(
                ruv_lines,
                paste("* [Error formatting RUV optimization results:", e$message, "]"),
                ""
            )
        }
    )

    ruv_lines
}
