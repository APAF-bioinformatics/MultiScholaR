# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

artifactWorkflowStateDehydrate <- function(
    dehydrate_fn,
    state_object,
    identity,
    generation_id
) {
    arguments <- list(state_object)
    if ("logical_key" %in% names(formals(dehydrate_fn))) {
        arguments$logical_key <- list(
            project_id = identity$project_id,
            omic_type = identity$omic_type,
            workflow_slug = identity$workflow_slug,
            stage_id = paste0("state_", substr(generation_id, 5L, 20L)),
            state_role = "s4_bundle",
            generation_id = generation_id
        )
    }
    do.call(dehydrate_fn, arguments)
}
