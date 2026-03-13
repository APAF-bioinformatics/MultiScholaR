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

# R/zzz.R

.onAttach <- function(libname, pkgname) {
    # Construct the message dynamically using the pkgname variable
    welcome_message <- paste0("Welcome to ", pkgname, "!\n\n")
    load_deps_message <- paste0(
        "IMPORTANT: Please run ", pkgname, "::loadDependencies() \n"
        , "           to ensure all required packages are installed and loaded.\n"
    )

    packageStartupMessage(
        "-------------------------------------------------------------------\n"
        , welcome_message
        , load_deps_message
        , "-------------------------------------------------------------------"
    )
}
