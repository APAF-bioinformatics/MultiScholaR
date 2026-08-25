# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactLipidomicsCodecDeclarations <- function() {
    list(
        "multischolar.s4.lipidomics_assay_data" = list(
            codec_id = "multischolar.s4.lipidomics_assay_data",
            codec_version = 1L,
            class_name = "LipidomicsAssayData",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        ),
        "multischolar.s4.lipidomics_da_results" = list(
            codec_id = "multischolar.s4.lipidomics_da_results",
            codec_version = 1L,
            class_name = "LipidomicsDifferentialAbundanceResults",
            payload_schema_id = "multischolar.rectangular",
            payload_schema_version = 1L
        )
    )
}
