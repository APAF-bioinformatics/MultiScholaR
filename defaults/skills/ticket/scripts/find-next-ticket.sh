#!/usr/bin/env bash
# Select the first pending ticket whose dependencies are completed.

set -euo pipefail

if [[ $# -ne 1 ]]; then
  echo "Usage: find-next-ticket.sh <tickets-index.json>" >&2
  exit 2
fi

tickets_index="$1"
if [[ ! -f "$tickets_index" ]]; then
  echo "ERROR: tickets index not found: $tickets_index" >&2
  exit 1
fi
if ! command -v jq >/dev/null 2>&1; then
  echo "ERROR: jq is required" >&2
  exit 1
fi

completed_ids="$(jq -c '[.tickets[] | select(.status == "completed" or .status == "complete") | .id]' "$tickets_index")"

jq -r --argjson completed_ids "$completed_ids" '
  .tickets[]
  | select(.status == "pending")
  | select(
      (.dependencies | length == 0) or
      (.dependencies | all(. as $dependency | $completed_ids | index($dependency) != null))
    )
  | .id
' "$tickets_index" | head -n 1
