#!/usr/bin/env python3
"""Report whether a ticket still needs a separate implementation plan."""

from __future__ import annotations

import json
import re
import sys
from pathlib import Path


def frontmatter(text: str) -> str:
    lines = text.splitlines()
    if not lines or lines[0].strip() != "---":
        raise ValueError("ticket must start with frontmatter")
    for index in range(1, len(lines)):
        if lines[index].strip() == "---":
            return "\n".join(lines[1:index])
    raise ValueError("frontmatter closing delimiter is missing")


def scalar(metadata: str, key: str) -> str | None:
    match = re.search(rf"^{re.escape(key)}:\s*(.*?)\s*$", metadata, re.MULTILINE)
    if not match:
        return None
    return match.group(1).strip().strip('"\'')


def main() -> int:
    if len(sys.argv) != 2:
        print("Usage: check-planning-needed.py <ticket-file>", file=sys.stderr)
        return 2
    path = Path(sys.argv[1])
    try:
        metadata = frontmatter(path.read_text(encoding="utf-8"))
    except (OSError, ValueError) as exc:
        print(json.dumps({"error": str(exc)}), file=sys.stderr)
        return 1

    explicit = scalar(metadata, "needs_planning")
    if explicit in {"true", "false"}:
        needed = explicit == "true"
        result = {
            "needs_planning": needed,
            "reason": f"Explicit frontmatter field: needs_planning={str(needed).lower()}",
            "confidence": "explicit",
        }
    else:
        plan = scalar(metadata, "plan")
        needed = not bool(plan)
        result = {
            "needs_planning": needed,
            "reason": "No plan reference" if needed else f"Plan referenced: {plan}",
            "confidence": "schema",
        }
    print(json.dumps(result))
    return 0


if __name__ == "__main__":
    raise SystemExit(main())
