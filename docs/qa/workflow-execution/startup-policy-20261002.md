# Startup policy amendment — 2026-10-02

The user explicitly set the production startup host-memory minimum to **4 GiB**.
This supersedes the 16 GiB startup field in the frozen WFX programme/qualification
manifest for subsequent WFX-004 members. Those historical authorities and failed
launch receipts remain unchanged; comparison requests bind this successor policy.

Only production `startup_headroom` changes: 17,179,869,184 -> 4,294,967,296 bytes.
Retain 16 GiB application hard cap, 15.5 GiB early stop, zero swap, eight threads,
12 GiB upper engine budget and 2 GiB continuous host reserve. Synthetic diagnostic
policies are unchanged. Apply the same amended policy to control and candidate.
No timed measurement preceded this amendment. Preserve both refused launch attempts;
the user's instruction authorizes a replacement first member under this policy.
