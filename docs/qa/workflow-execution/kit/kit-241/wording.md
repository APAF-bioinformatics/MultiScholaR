# KIT-241 wording (drafts for the owner's wording pass, VAL-41)

No text in the application changes. One NEWS entry, under `### Fixes` of the unreleased heading:

> Proteomics: the enrichment tab's 'Download All Results' now writes its archive itself. Before, it needed a
> separate `zip` program, which the installed Windows application does not have (and which some Linux computers
> lack), so on those computers the download produced no archive. This change does not alter the files in the
> archive (KIT-241).

- Last sentence: the implementation audit's I07 text, replacing "The files in the archive are unchanged", which read
  beside KIT-215's entry ('Not tested' now shown in the download) as if nothing in the download changed.
- The Linux clause is kept because G0 showed the download failing in a bare ubuntu:22.04 container without
  `/usr/bin/zip` (run 45577719573; stop rule 8 not met).
