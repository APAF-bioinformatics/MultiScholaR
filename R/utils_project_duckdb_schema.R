# MultiScholaR: Interactive Multi-Omics Analysis
# Copyright (C) 2024-2026 Ignatius Pang, William Klare, and APAF-bioinformatics
#
# This program is free software: you can redistribute it and/or modify
# it under the terms of the GNU Lesser General Public License as published by
# the Free Software Foundation, either version 3 of the License, or
# (at your option) any later version.

projectRegistrySchemaStatements <- function() {
    c(
        paste(
            "CREATE TABLE registry_metadata (",
            "singleton INTEGER PRIMARY KEY CHECK (singleton = 1),",
            "schema_id VARCHAR NOT NULL,",
            "schema_version INTEGER NOT NULL CHECK (schema_version >= 0),",
            "project_id VARCHAR NOT NULL,",
            "created_at VARCHAR NOT NULL,",
            "updated_at VARCHAR NOT NULL",
            ")"
        ),
        paste(
            "CREATE TABLE registry_migrations (",
            "migration_version INTEGER PRIMARY KEY CHECK (migration_version >= 1),",
            "migration_name VARCHAR NOT NULL UNIQUE,",
            "migration_checksum VARCHAR NOT NULL,",
            "applied_at VARCHAR NOT NULL,",
            "package_version VARCHAR NOT NULL",
            ")"
        ),
        paste(
            "CREATE TABLE registry_resource_settings (",
            "session_id VARCHAR NOT NULL,",
            "setting_name VARCHAR NOT NULL,",
            "requested_value VARCHAR NOT NULL,",
            "effective_value VARCHAR NOT NULL,",
            "unit VARCHAR NOT NULL,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (session_id, setting_name)",
            ")"
        ),
        paste(
            "CREATE TABLE projects (",
            "project_id VARCHAR PRIMARY KEY,",
            "root_locator VARCHAR NOT NULL CHECK (root_locator = '.'),",
            "status VARCHAR NOT NULL CHECK (status IN ('active', 'archived')),",
            "created_at VARCHAR NOT NULL,",
            "updated_at VARCHAR NOT NULL",
            ")"
        ),
        paste(
            "CREATE TABLE workflows (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "omic_type VARCHAR NOT NULL,",
            "omic_label VARCHAR NOT NULL,",
            "workflow_slug VARCHAR NOT NULL,",
            paste0(
                "status VARCHAR NOT NULL CHECK (status IN (",
                "'created', 'active', 'completed', 'failed', 'stale', 'archived')) ,"
            ),
            "created_at VARCHAR NOT NULL,",
            "updated_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id),",
            "UNIQUE (project_id, omic_type, workflow_slug),",
            "FOREIGN KEY (project_id) REFERENCES projects(project_id)",
            ")"
        ),
        paste(
            "CREATE TABLE workflow_runs (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "run_id VARCHAR NOT NULL,",
            paste0(
                "status VARCHAR NOT NULL CHECK (status IN (",
                "'created', 'running', 'completed', 'failed', 'cancelled')) ,"
            ),
            "action_id VARCHAR,",
            "started_at VARCHAR NOT NULL,",
            "completed_at VARCHAR,",
            "created_at VARCHAR NOT NULL,",
            "updated_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, run_id),",
            "UNIQUE (project_id, workflow_id, action_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id) REFERENCES ",
                "workflows(project_id, workflow_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE provenance_sources (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "run_id VARCHAR NOT NULL,",
            "source_id VARCHAR NOT NULL,",
            "source_role VARCHAR NOT NULL,",
            "source_uri VARCHAR,",
            "source_digest VARCHAR NOT NULL,",
            "source_size_bytes BIGINT,",
            "parser_id VARCHAR NOT NULL,",
            "parser_version VARCHAR NOT NULL,",
            "format_id VARCHAR NOT NULL,",
            "data_level VARCHAR NOT NULL,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, run_id, source_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id, run_id) REFERENCES ",
                "workflow_runs(project_id, workflow_id, run_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE provenance_parameters (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "run_id VARCHAR NOT NULL,",
            "parameter_id VARCHAR NOT NULL,",
            "parameter_name VARCHAR NOT NULL,",
            "value_json VARCHAR NOT NULL,",
            "value_digest VARCHAR NOT NULL,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, run_id, parameter_id),",
            "UNIQUE (project_id, workflow_id, run_id, parameter_name),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id, run_id) REFERENCES ",
                "workflow_runs(project_id, workflow_id, run_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE provenance_software (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "run_id VARCHAR NOT NULL,",
            "software_id VARCHAR NOT NULL,",
            "software_name VARCHAR NOT NULL,",
            "software_version VARCHAR NOT NULL,",
            "software_source VARCHAR NOT NULL,",
            "software_digest VARCHAR,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, run_id, software_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id, run_id) REFERENCES ",
                "workflow_runs(project_id, workflow_id, run_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE artifacts (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "artifact_id VARCHAR NOT NULL,",
            "run_id VARCHAR,",
            "generation_id VARCHAR NOT NULL,",
            "stage_id VARCHAR NOT NULL,",
            "state_role VARCHAR NOT NULL,",
            "hydration_ordinal INTEGER NOT NULL CHECK (hydration_ordinal >= 0),",
            "relative_path VARCHAR NOT NULL,",
            "codec_id VARCHAR NOT NULL,",
            "codec_version INTEGER NOT NULL CHECK (codec_version >= 1),",
            "payload_schema_id VARCHAR NOT NULL,",
            "payload_schema_version INTEGER NOT NULL CHECK (payload_schema_version >= 1),",
            "semantic_digest VARCHAR NOT NULL,",
            "byte_digest VARCHAR NOT NULL,",
            "row_count BIGINT,",
            "column_count BIGINT,",
            "payload_bytes BIGINT NOT NULL CHECK (payload_bytes >= 0),",
            paste0(
                "status VARCHAR NOT NULL CHECK (status IN (",
                "'staged', 'validated', 'committed', 'trashed')) ,"
            ),
            "created_at VARCHAR NOT NULL,",
            "updated_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, artifact_id),",
            paste0(
                "UNIQUE (project_id, workflow_id, generation_id, hydration_ordinal),"
            ),
            paste0(
                "FOREIGN KEY (project_id, workflow_id) REFERENCES ",
                "workflows(project_id, workflow_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE workflow_states (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "generation_id VARCHAR NOT NULL,",
            "parent_generation_id VARCHAR,",
            "logical_name VARCHAR NOT NULL,",
            "manifest_relative_path VARCHAR NOT NULL,",
            paste0(
                "status VARCHAR NOT NULL CHECK (status IN (",
                "'staged', 'current', 'historical', 'stale', 'failed')) ,"
            ),
            "created_at VARCHAR NOT NULL,",
            "updated_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, generation_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id) REFERENCES ",
                "workflows(project_id, workflow_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE artifact_dependencies (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "artifact_id VARCHAR NOT NULL,",
            "depends_on_artifact_id VARCHAR NOT NULL,",
            "relationship_type VARCHAR NOT NULL,",
            "ordinal INTEGER NOT NULL CHECK (ordinal >= 0),",
            "recorded_at VARCHAR NOT NULL,",
            paste0(
                "PRIMARY KEY (project_id, workflow_id, artifact_id, ",
                "depends_on_artifact_id, relationship_type),"
            ),
            paste0(
                "FOREIGN KEY (project_id, workflow_id, artifact_id) REFERENCES ",
                "artifacts(project_id, workflow_id, artifact_id),"
            ),
            paste0(
                "FOREIGN KEY (project_id, workflow_id, depends_on_artifact_id) ",
                "REFERENCES artifacts(project_id, workflow_id, artifact_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE run_artifacts (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "run_id VARCHAR NOT NULL,",
            "artifact_id VARCHAR NOT NULL,",
            "direction VARCHAR NOT NULL CHECK (direction IN ('input', 'output')),",
            "artifact_role VARCHAR NOT NULL,",
            "ordinal INTEGER NOT NULL CHECK (ordinal >= 0),",
            "recorded_at VARCHAR NOT NULL,",
            paste0(
                "PRIMARY KEY (project_id, workflow_id, run_id, direction, ",
                "artifact_role, ordinal),"
            ),
            paste0(
                "FOREIGN KEY (project_id, workflow_id, run_id) REFERENCES ",
                "workflow_runs(project_id, workflow_id, run_id),"
            ),
            paste0(
                "FOREIGN KEY (project_id, workflow_id, artifact_id) REFERENCES ",
                "artifacts(project_id, workflow_id, artifact_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE workflow_events (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "event_id VARCHAR NOT NULL,",
            "generation_id VARCHAR,",
            "run_id VARCHAR,",
            "event_type VARCHAR NOT NULL,",
            "event_status VARCHAR NOT NULL,",
            "details_json VARCHAR NOT NULL,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, event_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id) REFERENCES ",
                "workflows(project_id, workflow_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE workflow_figures (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "figure_id VARCHAR NOT NULL,",
            "generation_id VARCHAR NOT NULL,",
            "artifact_id VARCHAR,",
            "figure_role VARCHAR NOT NULL,",
            "relative_path VARCHAR NOT NULL,",
            "content_digest VARCHAR NOT NULL,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, figure_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id) REFERENCES ",
                "workflows(project_id, workflow_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE workflow_metrics (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "metric_id VARCHAR NOT NULL,",
            "generation_id VARCHAR,",
            "run_id VARCHAR,",
            "metric_name VARCHAR NOT NULL,",
            "numeric_value DOUBLE,",
            "value_json VARCHAR,",
            "unit VARCHAR,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, metric_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id) REFERENCES ",
                "workflows(project_id, workflow_id)"
            ),
            ")"
        ),
        paste(
            "CREATE TABLE registry_revisions (",
            "project_id VARCHAR NOT NULL,",
            "workflow_id VARCHAR NOT NULL,",
            "revision_id VARCHAR NOT NULL,",
            "generation_id VARCHAR,",
            "action_id VARCHAR NOT NULL,",
            "expected_parent_generation_id VARCHAR,",
            "revision_status VARCHAR NOT NULL,",
            "details_json VARCHAR NOT NULL,",
            "recorded_at VARCHAR NOT NULL,",
            "PRIMARY KEY (project_id, workflow_id, revision_id),",
            "UNIQUE (project_id, workflow_id, action_id),",
            paste0(
                "FOREIGN KEY (project_id, workflow_id) REFERENCES ",
                "workflows(project_id, workflow_id)"
            ),
            ")"
        ),
        paste0(
            "CREATE INDEX artifacts_generation_idx ON artifacts ",
            "(project_id, workflow_id, generation_id)"
        ),
        paste0(
            "CREATE INDEX events_workflow_idx ON workflow_events ",
            "(project_id, workflow_id, recorded_at)"
        )
    )
}
