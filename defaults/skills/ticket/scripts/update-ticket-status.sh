#!/usr/bin/env bash
# Atomically update one ticket status in tickets-index.json.

set -euo pipefail

if [[ $# -ne 3 ]]; then
  echo "Usage: update-ticket-status.sh <tickets-index.json> <ticket-id> <status>" >&2
  exit 2
fi

tickets_index="$1"
ticket_id="$2"
new_status="$3"

if [[ "$new_status" == "complete" ]]; then
  new_status="completed"
fi
if [[ ! "$new_status" =~ ^(pending|in_progress|completed|blocked)$ ]]; then
  echo "ERROR: invalid ticket status: $new_status" >&2
  exit 2
fi
if [[ ! -f "$tickets_index" ]]; then
  echo "ERROR: tickets index not found: $tickets_index" >&2
  exit 1
fi
if ! jq -e --arg id "$ticket_id" '.tickets[] | select(.id == $id)' "$tickets_index" >/dev/null; then
  echo "ERROR: ticket not found in index: $ticket_id" >&2
  exit 1
fi

temporary_file="$(mktemp)"
trap 'rm -f "$temporary_file"' EXIT
jq --arg id "$ticket_id" --arg status "$new_status" \
  '(.tickets[] | select(.id == $id) | .status) = $status' \
  "$tickets_index" > "$temporary_file"
mv "$temporary_file" "$tickets_index"
trap - EXIT

printf 'Updated %s -> %s\n' "$ticket_id" "$new_status"
