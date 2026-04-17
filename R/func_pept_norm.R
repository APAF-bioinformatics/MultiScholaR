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
# func_pept_norm.R
# ============================================================================
# Breadcrumb for the peptide-normalization surface that now lives in
# R/func_pept_s4_norm_methods.R. Keep the standalone exported helper here
# because it still owns the log2Transformation() public symbol.

#' @export
#' @title Log2 Transformation with Pseudo-count
#' @description Log 2 transformation with pseudo count
log2Transformation <- function(input_matrix) {
  pseudo_count <- min(input_matrix[input_matrix > 0], na.rm = TRUE) / 100
  positive_values <- input_matrix > 0 & !is.na(input_matrix)
  input_matrix[positive_values] <- input_matrix[positive_values] + pseudo_count
  log2(input_matrix)
}
