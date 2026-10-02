#!/usr/bin/env python3
"""Validate the constrained Markdown ticket schema used by this repository."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path
from typing import Any


REQUIRED_FIELDS: dict[str, type] = {
    "id": str,
    "title": str,
    "description": str,
    "status": str,
    "time_estimate": str,
    "dependencies": list,
}
VALID_STATUSES = {"pending", "in_progress", "completed", "blocked"}
REQUIRED_CHECKLIST_SECTIONS = (
    "Simulation Setup",
    "Acceptance Criteria",
    "Regression Checks",
    "Backward Compatibility",
    "Rollback Note",
)


def parse_scalar(value: str) -> Any:
    value = value.strip()
    if value in {"true", "false"}:
        return value == "true"
    if value in {"null", "~"}:
        return None
    if value.startswith("["):
        return json.loads(value)
    if len(value) >= 2 and value[0] == value[-1] and value[0] in {'"', "'"}:
        return value[1:-1]
    return value


def parse_frontmatter(text: str) -> tuple[dict[str, Any], str]:
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        raise ValueError("ticket must start with YAML frontmatter")
    try:
        closing = next(index for index in range(1, len(lines)) if lines[index].strip() == "---")
    except StopIteration as exc:
        raise ValueError("frontmatter closing delimiter is missing") from exc

    metadata: dict[str, Any] = {}
    frontmatter_lines = lines[1:closing]
    index = 0
    while index < len(frontmatter_lines):
        raw = frontmatter_lines[index]
        if not raw.strip() or raw.lstrip().startswith("#"):
            index += 1
            continue
        match = re.match(r"^([A-Za-z0-9_-]+):\s*(.*)$", raw)
        if not match:
            raise ValueError(f"unsupported frontmatter line: {raw}")
        key, raw_value = match.groups()
        if raw_value:
            metadata[key] = parse_scalar(raw_value)
            index += 1
            continue

        values: list[Any] = []
        index += 1
        while index < len(frontmatter_lines):
            item = re.match(r"^\s+-\s+(.+)$", frontmatter_lines[index])
            if not item:
                break
            values.append(parse_scalar(item.group(1)))
            index += 1
        metadata[key] = values

    return metadata, "\n".join(lines[closing + 1 :])


def section_body(content: str, heading: str) -> str | None:
    pattern = re.compile(
        rf"^##\s+{re.escape(heading)}\s*$\n(.*?)(?=^##\s+|\Z)",
        re.MULTILINE | re.DOTALL,
    )
    match = pattern.search(content)
    return match.group(1) if match else None


def validate(ticket_file: Path, tickets_index_file: Path | None) -> dict[str, Any]:
    errors: list[str] = []
    warnings: list[str] = []

    if not ticket_file.is_file():
        return {"valid": False, "errors": [f"ticket file not found: {ticket_file}"], "warnings": []}

    try:
        metadata, content = parse_frontmatter(ticket_file.read_text(encoding="utf-8"))
    except (OSError, ValueError, json.JSONDecodeError) as exc:
        return {"valid": False, "errors": [f"failed to parse ticket: {exc}"], "warnings": []}

    for field, expected_type in REQUIRED_FIELDS.items():
        if field not in metadata:
            errors.append(f"missing required field: {field}")
            continue
        if not isinstance(metadata[field], expected_type):
            errors.append(f"field '{field}' must be {expected_type.__name__}")
        elif expected_type is str and not metadata[field].strip():
            errors.append(f"field '{field}' cannot be empty")

    status = metadata.get("status")
    if isinstance(status, str) and status not in VALID_STATUSES:
        errors.append(f"invalid status '{status}'")
    estimate = metadata.get("time_estimate")
    if isinstance(estimate, str) and not re.fullmatch(r"\d+(?:\.\d+)?[hm]", estimate.lower()):
        errors.append(f"invalid time_estimate '{estimate}'")
    dependencies = metadata.get("dependencies", [])
    if isinstance(dependencies, list) and any(not isinstance(item, str) for item in dependencies):
        errors.append("all dependencies must be ticket-id strings")

    for heading in REQUIRED_CHECKLIST_SECTIONS:
        body = section_body(content, heading)
        if body is None:
            errors.append(f"missing required section: {heading}")
        elif not re.search(r"^\s*-\s+\[[ xX]\]", body, re.MULTILINE):
            errors.append(f"section '{heading}' must contain a checklist item")

    plan = metadata.get("plan")
    if isinstance(plan, str) and plan:
        plan_path = ticket_file.parent / plan.split("#", 1)[0]
        if not plan_path.is_file():
            errors.append(f"referenced plan does not exist: {plan_path}")

    if tickets_index_file is not None and tickets_index_file.is_file():
        try:
            index_data = json.loads(tickets_index_file.read_text(encoding="utf-8"))
            indexed = {item.get("id"): item for item in index_data.get("tickets", [])}
            ticket_id = metadata.get("id")
            if ticket_id not in indexed:
                errors.append(f"ticket '{ticket_id}' is missing from the index")
            else:
                entry = indexed[ticket_id]
                if entry.get("status") != status:
                    errors.append("frontmatter status does not match index status")
                if entry.get("file") != ticket_file.name:
                    errors.append("frontmatter ticket file does not match index file")
                if entry.get("dependencies", []) != dependencies:
                    errors.append("frontmatter dependencies do not match index dependencies")
            known_ids = set(indexed)
            for dependency in dependencies if isinstance(dependencies, list) else []:
                if dependency not in known_ids:
                    errors.append(f"dependency '{dependency}' is missing from the index")
        except (OSError, json.JSONDecodeError) as exc:
            errors.append(f"could not parse tickets index: {exc}")
    elif tickets_index_file is not None:
        warnings.append(f"tickets index not found: {tickets_index_file}")

    return {"valid": not errors, "errors": errors, "warnings": warnings}


def main() -> int:
    if len(sys.argv) not in {2, 3}:
        print("Usage: validate-ticket-schema.py <ticket-file> [tickets-index.json]", file=sys.stderr)
        return 2
    ticket_file = Path(sys.argv[1])
    tickets_index = Path(sys.argv[2]) if len(sys.argv) == 3 else None
    result = validate(ticket_file, tickets_index)
    print(json.dumps(result))
    return 0 if result["valid"] else 1


if __name__ == "__main__":
    raise SystemExit(main())
