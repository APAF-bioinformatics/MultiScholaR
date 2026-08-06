# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU Lesser General Public License for more details.
#
# You should have received a copy of the GNU Lesser General Public License
# along with this program.  If not, see <https://www.gnu.org/licenses/>.

# ============================================================================
# func_prot_qc_peptide.R
# ============================================================================
# Purpose: Peptide-level quality control and filtering functions
# 
# This file contains functions for peptide-level QC filtering, including
# intensity filtering, missing value filtering, replicate filtering, and
# q-value filtering. Functions in this file are used by mod_prot_qc_peptide.R
# and related QC modules.
#
# NOTE: Peptides are part of the proteomics workflow, hence "prot" prefix.
# This file contains peptide-specific QC functions within the proteomics context.
#
# Ownership: shared peptide-QC compatibility helpers. Filtering responsibilities
# live in func_prot_qc_peptide_*; peptide S4 behavior lives in func_pept_s4_*.
# ============================================================================

#' @keywords internal
resolvePeptideQcColumnName <- function(column_expr, env = parent.frame()) {
  if (is.null(column_expr)) {
    return(NULL)
  }

  if (is.character(column_expr) && length(column_expr) == 1) {
    return(column_expr)
  }

  column_value <- tryCatch(
    eval(column_expr, envir = env),
    error = function(e) NULL
  )
  if (is.character(column_value) && length(column_value) == 1) {
    return(column_value)
  }
  if (is.symbol(column_value)) {
    return(as.character(column_value))
  }

  if (is.symbol(column_expr)) {
    return(as.character(column_expr))
  }

  rlang::expr_text(column_expr)
}

#' @keywords internal
resolvePeptideQcColumnArgument <- function(column_expr, column_value, candidates = character(), env = parent.frame()) {
  if (is.character(column_expr) && length(column_expr) == 1) {
    return(column_expr)
  }

  if (is.symbol(column_expr)) {
    symbol_name <- as.character(column_expr)
    if (symbol_name %in% candidates) {
      return(symbol_name)
    }
  }

  forced_value <- tryCatch(
    force(column_value),
    error = function(e) NULL
  )
  if (is.character(forced_value) && length(forced_value) == 1) {
    return(forced_value)
  }
  if (is.symbol(forced_value)) {
    return(as.character(forced_value))
  }

  resolvePeptideQcColumnName(column_expr, env = env)
}































