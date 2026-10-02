# DIA-SEQ-01 — standalone development unnecessarily gated

Severity: medium. State blocked: implementation, not promotion.

DIA-005 revision 1.1 depends on DIA-004 and requires G0:iq; DIA-006 also requires
G0:iq. This prevents even isolated save/hydrate development before complete GUI
baseline qualification, contrary to the user's newly approved execution sequence.
The user explicitly approved early standalone DIA-005/006 work. Preserve all
scientific/performance gates before integration or acceptance.

Smallest repair: make DIA-005 depend on completed DIA-001/002, retain DIA-006's
DIA-005 dependency, and remove only their development-time G0 prerequisites.
All later ticket gates and thresholds remain unchanged. Closure requires schema
and graph checks plus explicit standalone-only restrictions. No empirical pass
or promotion follows from this amendment. Dirty worktree changes and all earlier
failed evidence remain untouched; no dependency updates, deletion or GUI wiring.
