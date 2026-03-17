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
# mod_lipid_norm_scale.R
# ============================================================================
# Purpose: Lipidomics scaling and normalization sub-module (Proteomics UX)
# ============================================================================

#' @title Lipidomics Scaling/Normalization Module
#' @description A Shiny module for between-sample normalization and visualization.
#' @name mod_lipid_norm_scale
NULL

#' @rdname mod_lipid_norm_scale
#' @export
#' @importFrom shiny NS tagList tabPanel br fluidRow column h5 textOutput imageOutput hr tags
#' @importFrom shinyjqui jqui_resizable
mod_lipid_norm_scale_ui <- function(id, plot_type = c("pca", "density", "rle", "correlation")) {
    ns <- shiny::NS(id)
    plot_type <- match.arg(plot_type)
    
    icon_map <- c(pca = "project-diagram", density = "chart-area", rle = "chart-bar", correlation = "th")
    label_map <- c(pca = "PCA", density = "Density", rle = "RLE", correlation = "Correlation")
    
    shiny::tabPanel(
        label_map[plot_type]
        , icon = shiny::icon(icon_map[plot_type])
        , shiny::br()
        # --- Row 1: Assay 1 ---
        , shiny::fluidRow(
            shiny::column(12
                , shiny::h5(shiny::textOutput(ns("assay1_label")), style = "margin-bottom: 5px;")
            )
        )
        , shiny::fluidRow(
            shiny::column(4
                , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                , shinyjqui::jqui_resizable(shiny::imageOutput(ns(paste0(plot_type, "_post_filter_assay1")), height = "300px"))
            )
            , shiny::column(4
                , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                , shinyjqui::jqui_resizable(shiny::imageOutput(ns(paste0(plot_type, "_post_norm_assay1")), height = "300px"))
            )
            , shiny::column(4
                , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                , shinyjqui::jqui_resizable(shiny::imageOutput(ns(paste0(plot_type, "_ruv_corrected_assay1")), height = "300px"))
            )
        )
        , shiny::hr()
        # --- Row 2: Assay 2 ---
        , shiny::fluidRow(
            shiny::column(12
                , shiny::h5(shiny::textOutput(ns("assay2_label")), style = "margin-bottom: 5px;")
            )
        )
        , shiny::fluidRow(
            shiny::column(4
                , shiny::tags$small("Post-Filtering", style = "display: block; text-align: center;")
                , shinyjqui::jqui_resizable(shiny::imageOutput(ns(paste0(plot_type, "_post_filter_assay2")), height = "300px"))
            )
            , shiny::column(4
                , shiny::tags$small("Post-Normalization", style = "display: block; text-align: center;")
                , shinyjqui::jqui_resizable(shiny::imageOutput(ns(paste0(plot_type, "_post_norm_assay2")), height = "300px"))
            )
            , shiny::column(4
                , shiny::tags$small("RUV-Corrected", style = "display: block; text-align: center;")
                , shinyjqui::jqui_resizable(shiny::imageOutput(ns(paste0(plot_type, "_ruv_corrected_assay2")), height = "300px"))
            )
        )
    )
}

#' @rdname mod_lipid_norm_scale
#' @export
mod_lipid_norm_scale_server <- function(id, workflow_data, experiment_paths, norm_data, plot_type) {
    shiny::moduleServer(id, function(input, output, session) {
        
        # Helper to render QC images
        render_qc_image <- function(assay_slot, stage_prefix) {
            shiny::renderImage({
                norm_data$plot_refresh_trigger
                assay_names <- norm_data$assay_names
                
                if (is.null(assay_names) || length(assay_names) < assay_slot) {
                    return(list(src = "", alt = "N/A"))
                }
                
                assay_name <- assay_names[[assay_slot]]
                safe_name <- gsub("[^A-Za-z0-9]", "_", tolower(assay_name))
                filename <- paste0(safe_name, "_", stage_prefix, "_", plot_type, ".png")
                
                img_path <- file.path(experiment_paths$lipid_qc_dir, filename)
                
                if (file.exists(img_path)) {
                    list(src = img_path, contentType = "image/png", width = "100%", height = "auto")
                } else {
                    list(src = "", alt = "Not generated")
                }
            }, deleteFile = FALSE)
        }
        
        # Bind outputs for 2 assays across 3 stages
        output[[paste0(plot_type, "_post_filter_assay1")]] <- render_qc_image(1, "pre_norm")
        output[[paste0(plot_type, "_post_norm_assay1")]]   <- render_qc_image(1, "post_norm")
        output[[paste0(plot_type, "_ruv_corrected_assay1")]] <- render_qc_image(1, "ruv_corrected")
        
        output[[paste0(plot_type, "_post_filter_assay2")]] <- render_qc_image(2, "pre_norm")
        output[[paste0(plot_type, "_post_norm_assay2")]]   <- render_qc_image(2, "post_norm")
        output[[paste0(plot_type, "_ruv_corrected_assay2")]] <- render_qc_image(2, "ruv_corrected")
        
        # Labels
        output$assay1_label <- shiny::renderText({
            if (length(norm_data$assay_names) >= 1) paste("Assay:", norm_data$assay_names[1]) else "Assay 1"
        })
        output$assay2_label <- shiny::renderText({
            if (length(norm_data$assay_names) >= 2) paste("Assay:", norm_data$assay_names[2]) else "Assay 2"
        })
    })
}

# <!-- APAF Bioinformatics | mod_lipid_norm_scale.R | Approved | 2026-03-17 -->
