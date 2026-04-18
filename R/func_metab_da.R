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
# func_metab_da.R
# ============================================================================
# Purpose: Metabolomics differential abundance analysis wrapper breadcrumb
#
# This file intentionally remains in `Collate:` as the stable metabolomics DA
# source identity after the bounded helper extractions. Live implementations
# now collate from the dedicated `R/func_metab_da_*` files, while metabolite S4
# DA methods and the shared `getCountsTable()` helper live in
# `R/func_metab_s4_da_results.R` and `R/func_metab_s4_objects.R`.
# ============================================================================
