# KIT-125 rev 1.1 (M-03): availability survey of public data sets per workflow

Date: 2026-11-17. Method: web research only. Every accession below was opened at its source (PRIDE archive API and FTP listing, Zenodo API, GEO, ENA portal, MetaboLights FTP and investigation files, Metabolomics Workbench REST, GitHub API). The files listed were seen in the listing; where a column header or a sample table is quoted, it was read from the file itself (a range request into the archive, or a small metadata file). Nothing large was downloaded. MassIVE could not be reached (every request timed out), so anything held only at MassIVE is marked "not confirmed at source". No KC data set was read.

## What the importers read (from the code)

- DIA and DIA_LIMPA: a DIA-NN main report, `.tsv` or `.parquet` (`importDIANNData`, R/func_prot_import_readers.R:145). Required columns `Protein.Group`, `Stripped.Sequence`, `Run`, and `Precursor.Normalised` or `Precursor.Quantity`; `Q.Value` and `PG.Q.Value` are mapped. The `pg_matrix.tsv` alone is not enough.
- LFQ MaxQuant: `proteinGroups.txt` with per-experiment `LFQ intensity <sample>` or `Intensity <sample>` columns (`importMaxQuantData`, :448). A search that was not split into experiments (only a summed `Intensity` column) cannot be read.
- LFQ FragPipe: `combined_protein.tsv` with `Protein ID` and `<sample> MaxLFQ Intensity` or `<sample> Intensity` columns (`importFragPipeData`, :294).
- TMT: a Proteome Discoverer protein export (`.xlsx`, `.csv`, `.tsv`, or a zip of several plexes) with an `Accession` column and `Abundance: F<n>: <channel>, <sample>` columns, or the `Abundance.` fallback (`importProteomeDiscovererTMTData`, R/func_prot_import_tmt.R:27). A PSM table does not qualify (it has `Master Protein Accessions`, not `Accession`).
- Lipidomics: vendor tables from MS-DIAL, LipidSearch, Progenesis QI, XCMS or Compound Discoverer, or "Other/Custom": a wide table with an id column and numeric sample columns (R/utils_omic_gui_config.R:1756, R/func_lipid_import_detection.R).
- Metabolomics: the same, without LipidSearch (R/utils_omic_gui_config.R:1981).
- Transcriptomics: one integer count table, an id column, a gene name column (required) and one column per sample; whole, non-negative numbers only (`previewTrxCountsImport`, `checkTrxImportCounts`, R/func_trx_counts.R).

Criteria applied: (a) quantified output in one of these formats is deposited; (b) at least 12 replicates of one condition, ideally in more than one recorded batch; (c) a licence or terms that allow reuse in a paper. "Technical" below means repeated measurements of one material (a digest, a pool, a reference RNA); the ticket's split unit is the biological unit, so a null split over technical replicates measures technical noise only, which is stated for each such set.

## DIA and DIA_LIMPA (one set serves both)

| Rank | Accession | Source | Processed files (format) | Replicates of one condition | Batches | Licence | Still to do | Verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | PXD024043 (Brunner et al. 2022, single HeLa cells, timsTOF) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2022/02/PXD024043/ | `DIANN1.8_SingleCells_CellCycle.zip` (628 MB; holds `20210919_DIANN_SingleCellOutput.tsv`, 1.53 GB, the DIA-NN 1.8 main report, plus `pg_matrix`, `pr_matrix`, `stats`). Also `DIANN1.8_diaPASEF_1ng_repetitions.zip` (110 MB; main report 296 MB, 15 runs of 1 ng HeLa). Report header read: has every required column. | From the stats file (435 runs), by the label in the run name: UB 203, G1 93, G2 52, NB 45, TB 41. The labels are cell-cycle states in the paper; not resolved further here. | Yes, as acquisition date in the run name: UB 65 runs on 2020-10-10 and 138 on 2020-11-26; NB on two dates. | CC0 | Unzip; build the design from run names. Single-cell data: very sparse, so expect heavy filtering. | usable |
| 2 | Zenodo 10.5281/zenodo.23615038 (mirror of MassIVE MSV000098622 / PXD066486, the MSstats+ K562 and CSF dilution series, Orbitrap Astral) | https://zenodo.org/records/23615038 | `K562_DIANN__report.parquet` (580 MB, DIA-NN 2, 40 runs), `K562_DIANN__report.pg_matrix.tsv.gz`, `K562_DIANN__K562_annotation.csv.gz`; the same for CSF (`CSF_DIANN__report.parquet`, 127 MB). | K562: 15 neat runs ("Condition1 Good") and 15 runs at 1:2 ("Condition2 Good"), plus 10 diluted runs. CSF: 15 neat, 16 at 1:2. All technical (one digest). | None recorded (one acquisition series with an order column). | CC BY 4.0 (the record); a third-party mirror, not the authors' release. The MassIVE original was not reachable. | None for format (parquet is read). Design from the annotation file. | usable (technical replicates only) |
| 3 | PXD054790 (The Pan-Cancer Proteome Atlas, Q Exactive HF DIA, three labs) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/05/PXD054790/ ; gzip mirror in Zenodo 23615038 (`TPCPA__diann-report.tsv.gz`, 11.3 GB) | `diann-report.zip` (11 GB; DIA-NN 1.8 report, 42 GB uncompressed, 1,236 runs; Deflate64, which common unzip tools do not read), `sample_info.txt`, `sample_labels.txt`, `protein_report.txt`. | Biological: e.g. colon primary cancer, lab Jimenez, 192 samples; breast 100; HeLa QC 38. | Yes: MS run date per sample (colon: 38, 16, 14, 123 and 1 samples on five dates); lab. | CC0 | Size: 42 GB report; subset by runs before import. Use the gzip mirror or a Deflate64-capable unzip. | usable with conversion (subsetting) |

Spike-in with known truth (DIA):

| Accession | Source | Processed files | Design | Licence | Verdict |
|---|---|---|---|---|---|
| Zenodo 10.5281/zenodo.8063173 (Gotti et al. 2021 benchmark; raw data PXD026600 at MassIVE MSV000087597) | https://zenodo.org/records/8063173 | `benchmarkingDIA.tsv` (531 MB): a DIA-NN `report.tsv` (header read: required columns present), FASTA and pipeline files. | UPS1 in an E. coli background; per the record and the paper, 8 UPS1 amounts with 3 replicates each; the record searched the "Overlapped" subset only (run names such as `RD139_Overlap_UPS1_0_1fmol_inj2`). | CC BY 4.0 | usable |
| Zenodo 10.5281/zenodo.14823673 (Bruderer et al. 2015 data, DIA-NN 1.8.1) | https://zenodo.org/records/14823673 | `report.zip` (435 MB, holds `report.tsv`, 1.85 GB) | Spiked-protein design of Bruderer 2015 (not re-read here). | CC BY 4.0 | usable (design not verified) |

Not usable: PXD070049 (LFQ benchmark "Generation Beta", CC0): raw, mgf, mzid and SDRF only; the DIA-NN analysis named in its protocol is not deposited. PXD028735 (Van Puyvelde 2022, CC0): raw files and one `identifications.tsv` only.

## LFQ MaxQuant

| Rank | Accession | Source | Processed files (format) | Replicates of one condition | Batches | Licence | Still to do | Verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | PXD002854 (Geyer et al. 2016, plasma proteome profiling, Q Exactive) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2016/02/PXD002854/ | Twelve MaxQuant search zips, each with `proteinGroups.txt`, `peptides.txt`, `summary.txt`, `parameters.txt`: e.g. `Search_15WorkflowReplicates.zip` (4.9 MB; 48 `LFQ intensity` columns), `Search_15RepeatedInjections.zip` (4.8 MB; 48), `Search_96LH.zip` (10 MB; 129), `Search_10Fingers.zip`. | 15 workflow replicates of one plasma sample (columns `LFQ intensity 1` to `15`); 15 repeated injections in the other zip. Technical. The same files carry 33 further columns from other experiments (F1_1 to M5_3, Top6_Top20), which are dropped. | None. | EBI terms of use (no explicit open licence). | Select the 15 columns. | usable (technical replicates only) |
| 2 | PXD042233 (Webel et al. 2024, 7,444 HeLa QC runs, MaxQuant 1.6.12) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2023/12/PXD042233/ | `proteinGroups_single_dumps.zip` (2.5 GB; one CSV per run, each run searched alone; columns read: `Protein IDs`, `Majority protein IDs`, `Gene names`, `Q-value`, `Score`, `Intensity`, ...), `geneGroups_aggregated.zip` (176 MB), `pride_metadata.csv` (instrument, serial number, date per run), SDRF. | Thousands of HeLa runs. Technical (one cell line, standard digest). | Yes: instrument and date per run (metadata read in part: e.g. 353 of the first 514 rows on one Q Exactive, 39, 34, 23 on others). | CC0 | Merge per-run files into one wide `Intensity <run>` table in proteinGroups shape (no LFQ across runs exists: every run was searched separately); select runs by instrument and date. | usable with conversion |

Spike-in with known truth (MaxQuant): PXD000279 (Cox et al. 2014, MaxLFQ benchmarks; https://ftp.pride.ebi.ac.uk/pride/data/archive/2014/09/PXD000279/). `dynamicrangebenchmark.zip` (11 MB) holds `proteinGroups.txt` with `LFQ intensity UPS1_01..04` and `UPS2_01..04`; `proteomebenchmark.zip` (37 MB) holds `proteinGroups.txt` with `LFQ intensity H1..H3`, `L1..L3` (human with E. coli at two levels). EBI terms of use. Usable as is; small (4 vs 4, 3 vs 3).

Not usable for MaxQuant: PXD057758's `*_Results-MaxQuant.zip` (Ecoli, Lim, SingleCell) contain `proteinGroups.txt` with only a summed `Intensity` column (no per-sample columns; the search was not split into experiments). PXD007145 holds MaxQuant output in `.rar` archives for phospho benchmarks (not opened).

## LFQ FragPipe

| Rank | Accession | Source | Processed files (format) | Replicates of one condition | Batches | Licence | Still to do | Verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | Zenodo 10.5281/zenodo.17607826 (FragPipe reprocessing of the COMBAT blood atlas, PXD023175, timsTOF) | https://zenodo.org/records/17607826 | `combined_protein.tsv` (4.2 MB; 92 `MaxLFQ Intensity` columns), `complete92.fp-manifest`, `sdrf.tsv`, `timstof_combat.workflow`. | 30 healthy volunteers (HC), 62 sepsis. Biological. | One plate (all runs named `plate5`); none usable. | CC BY 4.0 (the record); source PXD023175 is CC0. A third-party reprocessing, not the authors' own. | None for format. | usable |
| 2 | PXD057758 (PIP-ECHO paper; reprocessing of PXD014415 "Lim" and an E. coli two-proteome set) | https://ftp.pride.ebi.ac.uk/pride/data/archive/2025/10/PXD057758/ | `Lim_Results-FragPipe.zip` (494 MB) and `Ecoli_Results-FragPipe.zip` (247 MB): six FragPipe output folders each, each with `combined_protein.tsv`, `experiment_annotation.tsv`, `fragpipe.workflow`. Lim: 40 `MaxLFQ Intensity` columns. | Lim: 20 human-only runs (`a11849` to `a11868`) and 20 human+yeast runs. Technical. The annotation labels all 40 runs as condition A, so the design comes from the run names. E. coli set: 10 human-only and 10 human+E. coli runs (fewer than 12). | None. | CC0 | Rebuild the design from file names; pick one of the six output folders. | usable (technical replicates only) |

Spike-in with known truth (FragPipe): the same PXD057758. Lim: yeast present in 20 runs and absent in 20 (presence truth rather than a ratio). E. coli set: 10 human vs 10 human+E. coli at 10:1. `IonStar_Results-FragPipe.zip` (3.7 GB; IonStar is E. coli with human spike at five levels, PXD003881) was listed but not opened. CC0. Usable with design rebuilding.

## TMT (Proteome Discoverer)

| Rank | Accession | Source | Processed files (format) | Replicates of one condition | Batches | Licence | Still to do | Verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | Zenodo 10.5281/zenodo.14767905 (msqrob2TMT data; `spikein1` is the MSstatsTMT controlled mixture of Huang et al. 2020, said to come from MassIVE RMSV000000265) | https://zenodo.org/records/14767905 | `spikein1_psms.txt` (239 MB): Proteome Discoverer PSM export, columns `Master Protein Accessions`, `Abundance: 126` to `Abundance: 131`, `Spectrum File`, ...; `spikein1_annotations.csv` (150 rows); `processed.zip` (271 MB, not opened). | Annotation read: 5 mixtures (TMT10 plexes) x 10 channels x 3 technical-replicate runs; conditions Norm, 0.125, 0.5, 0.667, 1 (30 channel-runs each). Run names `161117_SILAC_HeLa_UPS1_TMT10_Mixture<k>_0<r>.raw`: a HeLa background with UPS1 spiked at known levels, so the HeLa proteins are a null across all 50 channels. Technical (one HeLa digest). | Yes: the mixture (plex). | CC BY 4.0 (the record). The MassIVE original was not reachable. | Conversion: roll PSMs up to proteins and write a PD-style protein table (`Accession`, `Abundance: F<n>: <channel>, <sample>`), one file per plex. | usable with conversion |

Spike-in with known truth (TMT): the same set (UPS1 at four levels against a constant HeLa background, five plexes). It serves both the null and the known-truth check, both after conversion.

No second TMT candidate passed. PXD070151 (UPS1/UPS2 TMT10 benchmark, Proteome Discoverer 2.4, CC0) has `..._ProteinGroups.txt` (150 KB), but its header read here has no abundance columns, and the design is one plex with two channels per condition. PXD073725 (PD 2.4 `Proteins.txt`, CC0) has three dose groups in one plex (fewer than 12 per group). PXD036799 (inter-batch isobaric benchmark, CC0) has only `InterBatch_Peptides.csv.csv` (69 MB, Comet); not opened, and not a PD protein export. PXD027248 has PD `.msf` files only.

## Lipidomics

No public set was found in a vendor format the importer reads (MS-DIAL, LipidSearch) that also has 12 or more replicates of one condition and an open licence. The sets that do have the replicate structure are at the Metabolomics Workbench, whose every study (4,714 checked) carries CC BY-NC-ND 4.0; their quantified tables come from the Workbench REST service as generic matrices.

| Rank | Accession | Source | Processed files (format) | Replicates of one condition | Batches | Licence | Still to do | Verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | ST002402 (4D TIMS lipidomics of human blood, Mainz) | https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002402 ; table from `rest/study/analysis_id/AN003915/untarg_data` | Untargeted feature table, positive mode: 421 rows (samples), 2,260 features named `m/z_RT`; negative and a third analysis (AN003914, AN003916) also exist. The named-metabolite datatable is empty. | NIST human plasma: 32 extractions per plate; plates 1 (intraday), 2 and 3 (interday). NIST serum: 32. Technical (one reference material). | Yes: extraction plate (3 plates of plasma). | CC BY-NC-ND 4.0 | Transpose to features x samples and import as Custom. Licence: see Limits. | usable with conversion; licence needs a decision |
| 2 | ST003456 (sepsis serum lipidomics, Leibniz HKI) | https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST003456 ; `untarg_data` of AN005677 | Untargeted table, positive mode: 452 rows, 4,810 features; named datatable 443 lipids (pos) and 54 (neg). | Controls at "Pre": 38; at "-1": 38; at "-3": 35. Biological. QC injections: 108. | Not recorded. | CC BY-NC-ND 4.0 | Transpose; import as Custom. | usable with conversion; licence needs a decision |

Checked and rejected: MetaboLights MTBLS15755 (CC0, targeted lipids, 16 NIST SRM 1950 and 32 pooled QC samples) and MTBLS13538 (CC0, 57 controls): their MAF files have no values. MTBLS10722 (EBI terms): 10 volunteers only. Zenodo 7097871 is a genuine MS-DIAL multi-batch export but of wastewater metabolomics, with no 12-replicate condition.

## Metabolomics

| Rank | Accession | Source | Processed files (format) | Replicates of one condition | Batches | Licence | Still to do | Verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | MTBLS404 (Thevenot et al., "sacurine", human urine, LC-HRMS negative mode) | https://ftp.ebi.ac.uk/pub/databases/metabolights/studies/public/MTBLS404/ | `m_MTBLS404_sacurine_v2_maf.tsv` (227 KB; 120 identified metabolites x 234 sample columns, 210 with values), `s_MTBLS404.txt`, `a_MTBLS404_LC-MS_sacurine.txt`. | 184 adult volunteers (biological; one population, which varies in age, BMI and sex, all recorded); 26 pooled QC injections. | Yes: `Factor Value[Batch]` ne1 / ne2 (90 and 94 subjects; 14 and 12 QCs), and injection order. | EMBL-EBI terms of use (no explicit open licence). | Drop the 27 MAF annotation columns and the 24 empty blank columns; import as Custom with `database_identifier` or `metabolite_identification` as the id. | usable with conversion |
| 2 | ST002092 (Addressing batch effects in large-scale metabolomics, C. elegans) | https://www.metabolomicsworkbench.org/data/DRCCMetadata.php?Mode=Study&StudyID=ST002092 ; `untarg_data` of AN003416 (also AN003417, AN003418) | Untargeted table, RP positive: 147 rows, 3,953 features (`m/z_RT`). | Strain PD1074: 42 samples (biological). Pools: 12 PD1074, 12 mutant. | Yes: the batch prefix of the sample name (PD1074: 6, 7, 8, 6, 8, 7 samples in b1 to b6). | CC BY-NC-ND 4.0 | Transpose; import as Custom. | usable with conversion; licence needs a decision |

Rejected: MTBLS79 (multi-batch DIMS benchmark, 8 batches): one species has 10 biological samples (fewer than 12; the rest are repeats). ST002233 (batch variation of plasma QC, 268 samples) and ST004266: their Workbench tables are empty.

## Transcriptomics

| Rank | Accession | Source | Processed files (format) | Replicates of one condition | Batches | Licence | Still to do | Verdict |
|---|---|---|---|---|---|---|---|---|
| 1 | ENA PRJEB5348 (Schurch et al. 2016, 48-replicate yeast RNA-seq) with counts in github.com/bartongroup/profDGE48 | https://www.ebi.ac.uk/ena/browser/view/PRJEB5348 ; https://github.com/bartongroup/profDGE48/tree/master/Preprocessed_data | `WT_countdata.tar.gz` (1.5 MB; 48 files `WT_repNN_..._tophat2.0.5.bam.gbgout`, two columns gene and count, 7,131 rows each), `Snf2_countdata.tar.gz`; `Bad_replicate_identification/exclude.lst` names 6 WT replicates as bad. ENA: 672 runs (48 + 48 samples x 7 lanes). | WT: 48 biological replicates (42 after the authors' exclusion list). | None recorded (every sample is pooled over the 7 lanes). | Repository MIT (code); the counts carry no separate licence; ENA data are open. | Join the 48 files into one table; add a gene name column (the importer requires one; the id can be repeated). | usable with conversion |
| 2 | GSE47774 (SEQC/MAQC-III) | https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE47774 | Per-site count tables, e.g. `GSE47774_SEQC_ILM_BGI.txt.gz` (16 MB), `_AGR` (13 MB), `_CNL` (15 MB), `_COH` (7.1 MB), `_MAY` (16 MB), `_NVS` (14 MB), and LIF/ROC files; integer counts by RefSeq transcript, one column per library x lane x flow cell. | Sample A (Universal Human Reference RNA): libraries per Illumina site read from the headers: AGR 4, BGI 5, COH 4, MAY 5, NVS 4 (CNL not read) = 22 libraries. Technical (one reference RNA). | Yes: sequencing site. | GEO (NCBI places no restrictions on use). | Sum lanes per library; aggregate transcripts to genes or keep transcript ids; add a name column. | usable with conversion (technical replicates only) |

## Summary

| Workflow | Usable public set | Best candidate | Biological or technical replicates | Batch recorded | Licence | Spike-in (proteomics) |
|---|---|---|---|---|---|---|
| DIA | yes | PXD024043 (DIA-NN 1.8 report) | single cells (biological) | yes (date) | CC0 | Zenodo 8063173 (Gotti, DIA-NN report): usable |
| DIA_LIMPA | yes (same report) | PXD024043 | as DIA | yes | CC0 | as DIA |
| LFQ MaxQuant | yes, technical only as deposited; batched only after conversion | PXD002854 (15 workflow replicates); PXD042233 (per-run dumps, conversion) | technical | PXD042233 only | EBI terms; CC0 | PXD000279: usable |
| LFQ FragPipe | yes | Zenodo 17607826 (30 healthy volunteers) | biological | no | CC BY 4.0 | PXD057758: usable after rebuilding the design |
| TMT | only with conversion | Zenodo 14767905 spikein1 (PD PSMs to PD protein table) | technical | yes (plex) | CC BY 4.0 | the same set, after conversion |
| Lipidomics | only with conversion and under CC BY-NC-ND | ST002402 (NIST plasma, 3 plates); ST003456 | technical; biological | ST002402 only | CC BY-NC-ND | n/a |
| Metabolomics | yes, with conversion | MTBLS404 (2 batches); ST002092 (6 batches) | biological | yes | EBI terms; CC BY-NC-ND | n/a |
| Transcriptomics | yes, with conversion | Schurch WT (48); SEQC sample A (22 libraries, sites) | biological; technical | SEQC only | open | n/a |

No workflow is without any candidate. Two have no set that is usable as deposited: TMT (PSM-level only, so a PSM-to-protein conversion is needed, which is part of what is tested) and lipidomics (no native-vendor set; the candidates are generic matrices under CC BY-NC-ND). If the licence is not accepted, lipidomics has no usable public set and is reported as not checked on real data. LFQ MaxQuant has no batched set as deposited.

## Limits of the survey

- Web research only; nothing was imported or run. A file counts as present when it was seen in the listing and, where stated, its header or annotation was read; no file was loaded end to end, so an import can still fail on content.
- MassIVE timed out on every request. PXD026600, PXD066486 and RMSV000000265 are known only through their ProteomeCentral records and the Zenodo mirrors; the mirrors are third-party copies (Zenodo 23615038 and 23168773 by one depositor; 17607826 a reprocessing; 14767905 the msqrob2TMT authors).
- Repository searches were keyword-based (PRIDE v3 search, Zenodo, EBI search for MetaboLights, the Workbench study list). Sets that do not name their output files or their design in searchable text can have been missed, in particular TMT sets with PD protein exports and lipidomics sets with MS-DIAL tables.
- Licences: the Metabolomics Workbench licenses every study CC BY-NC-ND 4.0. Using the data to compute and report a benchmark is not redistribution of an adapted table, but the "no derivatives" term and redistribution of converted tables need a decision by the user. EMBL-EBI terms of use (MTBLS404, PXD002854, PXD000279) are not an explicit open licence.
- Several best candidates offer only technical replicates of one material (K562 digest, plasma workflow replicates, HeLa QC, the TMT HeLa background, NIST plasma, SEQC sample A). Under the ticket's split unit these measure technical noise only. The sets with biological replicates of one condition are PXD024043 (single cells), Zenodo 17607826, ST003456, MTBLS404, ST002092 and the Schurch set.
- Batch-aligned splits (rev 1.1, second arm) are possible where batches are recorded: PXD024043 (two dates), PXD054790 (run dates), PXD042233 (instrument, date), Zenodo 14767905 (five plexes), ST002402 (three plates), MTBLS404 (two batches), ST002092 (six batches) and SEQC (sites).
- Scratch listings and the small metadata files read are under the session scratchpad (`k125-survey/data`); they are not part of the repository.
