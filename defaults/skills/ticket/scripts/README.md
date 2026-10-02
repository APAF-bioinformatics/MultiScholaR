# Local ticket helpers

These four provider-neutral metadata scripts were restored on 2026-09-11 from
`MultiScholaR-performance-archive-20260911-5ywjc71f/source/defaults/skills/ticket/scripts/`
so the explicitly requested `$openai-ticket` skill works on restored janitor.
They were reviewed and restored without modification. They are not runtime code,
benchmark drivers, or a restoration of the performance architecture.

- `find-next-ticket.sh`: select the first pending ticket with completed dependencies.
- `validate-ticket-schema.py`: validate the repository's constrained frontmatter,
  required checklist sections and index consistency (not arbitrary YAML).
- `check-planning-needed.py`: report the explicit planning field/plan reference.
- `update-ticket-status.sh`: update the index; use the skill's wrapper to update
  both index and frontmatter together.

Requirements: Bash, jq and Python 3. Invoke through the installed skill's
`next-openai-ticket.sh` and `update-openai-ticket-status.sh` wrappers.
Gate outcomes, scientific validity and performance acceptance require the ticket's
own evidence checks; metadata validation does not establish them.
