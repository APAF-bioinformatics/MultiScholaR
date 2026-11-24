# KIT-193 statements gate: approved commits and `SPEC.sha256`

Ticket section 12 (the approved commit and `SPEC.sha256`); section 5.3. Registers in `tickets/kit/audits/`.

| Step | Commit | Spec files | Approved by |
|---|---|---|---|
| First freeze | c783e1ae (hashes of the Spec files at 4183eb10) | `FilterCount`, `ReplayRule`, `GoClosure`, `QvalueBh`, `Permutation`, `RuvAlgebra`, `RuvDownstream`, `RuvIIIC` | `AUDIT-KIT-193-statements-20261124.md` (S01 to S07 applied; closure check equal to its values) |
| Amendment (5.3 item 4) | 5474e512 | `QvalueBh` (one line replaced: `bh_sorted_eq_stmt` gains `0 ≤ p j`; two statements added) | statements gate re-review 1 (same register); the auditor confirmed the commit |
| New file (P03) | b8a62ed3 | `SharperWitnesses` (one line added) | `AUDIT-KIT-193-proofs-20261124.md` P03, adopted by the coordinator; confirmed in `AUDIT-KIT-193-proofs2-20261124.md` A.1 row 1a |
| Drafts for review | b32c94ca | `RemovedFactors`, `RuvIIICPattern`, `RuvControlRules` (three lines added, not approved) | none: reviewed in `AUDIT-KIT-193-proofs2-20261124.md` section B |
| Freeze of R-I1, R-J3, R-N1 (Q09) | 2ab4c189 | the three drafts after the edits Q06 to Q08 (the optional statements included); their three lines replaced | **awaiting** the auditor's confirmation by reading (Q09) |

`SPEC.sha256` at 2ab4c189 (`sha256sum -c`: twelve OK):

```
6dc7b530b0b032c849d1b8aa25f8e4107bb80335e51a44659fc88f151b7578ce  MultiScholaRProofs/Spec/FilterCount.lean
615d39d0afeeccbbb0d49e76560e926cded79f9c2f7b7d5a64c8414ad17502a2  MultiScholaRProofs/Spec/GoClosure.lean
bad08e62064ce0b417c666c06541f453fcd85c588ca4b893f4d7b1f2cb657929  MultiScholaRProofs/Spec/Permutation.lean
114d9534374c150f7df92fafc6b6de95124ea53b4a8d65040e7947ebe6599ace  MultiScholaRProofs/Spec/QvalueBh.lean
f9ebc52e44b5d6ddf31e738d45dad6aad32ee9fc0b6cedf3cceaa1d0c40f1962  MultiScholaRProofs/Spec/ReplayRule.lean
9918d44dce7274cd81abea4bb2b63d44ba5bd00fb9564a39fc8d22e825c9bfc5  MultiScholaRProofs/Spec/RuvAlgebra.lean
3890bfdaf81ba73c31f342771ca897559a8a7c6d8588b3fc97feadce36b04c6f  MultiScholaRProofs/Spec/RuvDownstream.lean
69dff74fd1a70f1b1e55d857d730cb6dbbf6a94f349b8df5bf4246b5e9661744  MultiScholaRProofs/Spec/RuvIIIC.lean
1c6a728b81d468d588e6d46a974ea86d763b47e2e580c4b61c5a2d5bb62598da  MultiScholaRProofs/Spec/SharperWitnesses.lean
d1736b7e540626a1c8f3a0b63663a8533fb71201f922dbe2c7d8bfbfb6543cfd  MultiScholaRProofs/Spec/RemovedFactors.lean
51a14205c7a75ec12c8e758aab46300a9e3633b1bffd37097bb0201444077150  MultiScholaRProofs/Spec/RuvIIICPattern.lean
a67a6d26ccb8acab328b4a342ab537efec11da00337b18d1b9dfe28bb3700538  MultiScholaRProofs/Spec/RuvControlRules.lean
```

The last three lines replace `f8b60b8c…`, `3b6b8e8c…` and `8383aac4…` (the drafts at b32c94ca). The first nine are
unchanged since their approval. No proof of R-I1, R-J3 or R-N1 starts before the confirmation of 2ab4c189.
