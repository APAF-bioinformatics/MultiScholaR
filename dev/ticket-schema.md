# Ticket Schema

All tickets in this repository use YAML frontmatter for machine-readable metadata.

## Required Fields

```yaml
---
ticket_id: string       # Unique ID. Format: PREFIX-NNN (e.g., RUV-001, LIP-003)
title: string            # Short title
status: string           # One of: pending, in_progress, completed, blocked, cancelled
priority: string         # One of: P0, P1, P2
depends_on: list         # List of ticket_id strings this ticket depends on (empty list if none)
---
```

## Optional Fields

```yaml
verified_at: date        # Date the ticket's current-state claims were last verified against live code
source_plan: string      # Path to the plan document this ticket derives from
series: string           # Series identifier for grouping related tickets (e.g., RUV, LIP)
created_at: date         # Date ticket was created
completed_at: date       # Date ticket was marked done
```

## Status Lifecycle

```
pending → in_progress → completed
  ↓           ↓
blocked     cancelled
```

- **pending**: All dependencies satisfied, can be picked up
- **in_progress**: Actively being worked on
- **completed**: Acceptance criteria met, code merged
- **blocked**: Waiting on external dependency or decision
- **cancelled**: No longer needed

## File Locations

Tickets live alongside their context — there is no single `tickets/` directory. The central index at `dev/tickets.md` aggregates across all locations.

| Series | Location | Source |
|--------|----------|--------|
| LIP    | `dev/TICKET-NNN_*.md` | Lipidomics scaffolding plan |
| RUV    | `dev/audit/RUV/tickets/RUV-NNN-*.md` | RUV math audit fix plan |
