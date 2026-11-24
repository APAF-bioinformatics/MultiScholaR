# KIT-167, pass 2 (coordinator's rulings of 2026-11-24): the public inputs the implementation still needs, fetched on a
# runner only. Sections, each recorded whatever the others do (a failed section writes its error and the job goes on):
#   A  which PRIDE project holds each strain's search FASTA (DIA and DDA projects listed; size and SHA-256 of each
#      `BPA_SEP_<strain>_Contaminants.fasta` found) - ruling (a);
#   F4 the GEO annotation files of the four strains (SHA-256, size; GTF feature types and attribute keys; GFF3 counts);
#   CUT a cut of the B36 GEO GFF3 and GTF (the loci of the parser's cases) with the expected protein per locus derived
#      by awk, independently of the package's parser - ruling (e);
#   F6 the reference order: each shipped reference's INSDC assembly as UniProt names it, its GFF from NCBI, the share of
#      entries whose `reference_locus` joins a protein-coding locus (`locus_tag`, else `old_locus_tag`); pass >= 90 %;
#   F3 the authors' accessions in `Results<strain>.xlsx` resolved to UniParc ids through their UniProtKB entries;
#   DEMO the demo mapping tables of the four strains with the pinned DIAMOND and the shipped references (the
#      real-engine lane: KPC2 also mapped at threads 1 and 4, decision tables compared).
# Writes only into <out>/ (tables of identifiers, counts and hashes; the cut annotation lines; the demo tables).
# Usage: Rscript pass2-inputs.R <out dir>

arguments <- commandArgs(trailingOnly = TRUE)
out <- normalizePath(if (length(arguments)) arguments[[1L]] else "pass2-out", mustWork = FALSE)
dir.create(out, recursive = TRUE, showWarnings = FALSE)
work <- file.path(tempdir(), "k167-pass2")
dir.create(work, recursive = TRUE, showWarnings = FALSE)
options(timeout = 600)
sha256 <- function(path) digest::digest(file = path, algo = "sha256")
log_line <- function(...) cat(sprintf("KIT167 %s\n", paste(..., collapse = " ")))
write_tsv <- function(table, name) utils::write.table(table, file.path(out, name), sep = "\t", quote = FALSE,
    row.names = FALSE, na = "NA")
fetch <- function(url, dest) {
    for (attempt in 1:3) {
        ok <- tryCatch({
            utils::download.file(url, dest, mode = "wb", quiet = TRUE, method = "libcurl")
            file.exists(dest)
        }, error = function(e) FALSE)
        if (ok) return(TRUE)
        Sys.sleep(30)
    }
    FALSE
}
read_text <- function(url) tryCatch(paste(readLines(url, warn = FALSE), collapse = "\n"), error = function(e) "")
listing <- function(folder) {
    text <- read_text(folder)
    files <- unique(sub("^href=\"(.*)\"$", "\\1", regmatches(text, gregexpr("href=\"[^\"?]+\"", text))[[1L]]))
    files[!grepl("^/|^\\.\\.", files)]
}
section <- function(name, code) {
    status <- tryCatch({
        code
        "ok"
    }, error = function(e) {
        log_line(name, "error", conditionMessage(e))
        paste("error:", conditionMessage(e))
    })
    log_line(name, "status", status)
    write_tsv(data.frame(section = name, status = status), sprintf("status-%s.tsv", name))
}
suppressMessages(pkgload::load_all(".", quiet = TRUE, export_all = TRUE, helpers = FALSE))
strains <- data.frame(strain = c("B36", "KPC2", "BPH2900", "HKU419"),
    dia = c("PXD016840", "PXD016672", "PXD016504", "PXD016913"),
    dda = c("PXD020545", "PXD020839", "PXD020791", "PXD020863"),
    geo = c("GSE152966", "GSE152843", "GSE152835", "GSE152824"),
    reference = c("UP000000625", "UP000000265", "UP000008816", "UP000000750"), stringsAsFactors = FALSE)
located <- list()

# --- A: the projects' folders and the search FASTAs ---------------------------------------------------------------
pride_folder <- function(accession) {
    text <- read_text(sprintf("https://www.ebi.ac.uk/pride/ws/archive/v3/projects/%s", accession))
    date <- regmatches(text, regexpr("\"publicationDate\"[[:space:]]*:[[:space:]]*\"[0-9]{4}-[0-9]{2}", text))
    candidates <- character()
    if (length(date)) {
        ym <- regmatches(date, regexpr("[0-9]{4}-[0-9]{2}$", date))
        candidates <- sprintf("https://ftp.pride.ebi.ac.uk/pride/data/archive/%s/%s/%s/", substr(ym, 1, 4),
            substr(ym, 6, 7), accession)
    }
    candidates <- c(candidates, as.vector(outer(sprintf("%02d/%s/", 1:12, accession), 2019:2024,
        function(tail, year) sprintf("https://ftp.pride.ebi.ac.uk/pride/data/archive/%d/%s", year, tail))))
    for (folder in candidates) {
        files <- listing(folder)
        if (length(files)) return(list(folder = folder, files = files))
    }
    list(folder = NA_character_, files = character())
}
section("A", {
    rows <- list()
    for (project in c(strains$dia, strains$dda)) {
        found <- pride_folder(project)
        keep <- found$files[grepl("\\.fasta$|^Results.*\\.xlsx$|MQtxt\\.zip$", found$files, ignore.case = TRUE)]
        log_line("A", project, found$folder, length(found$files), "files;", paste(keep, collapse = " "))
        for (file in keep) {
            head <- tryCatch(curlGetHeaders(paste0(found$folder, file)), error = function(e) character())
            size <- suppressWarnings(as.numeric(sub("(?i)^content-length:[[:space:]]*", "",
                grep("(?i)^content-length:", head, value = TRUE, perl = TRUE)[1L], perl = TRUE)))
            rows[[length(rows) + 1L]] <- data.frame(project = project, folder = found$folder, file = file,
                bytes = size, stringsAsFactors = FALSE)
        }
        if (!length(keep)) rows[[length(rows) + 1L]] <- data.frame(project = project, folder = found$folder,
            file = NA_character_, bytes = NA_real_)
    }
    files <- do.call(rbind, rows)
    fasta <- list()
    for (i in seq_len(nrow(strains))) {
        name <- sprintf("BPA_SEP_%s_Contaminants.fasta", strains$strain[i])
        hits <- files[!is.na(files$file) & tolower(files$file) == tolower(name), , drop = FALSE]
        for (j in seq_len(nrow(hits))) {
            dest <- file.path(work, sprintf("%s-%s.fasta", strains$strain[i], hits$project[j]))
            ok <- fetch(paste0(hits$folder[j], hits$file[j]), dest)
            fasta[[length(fasta) + 1L]] <- data.frame(strain = strains$strain[i], project = hits$project[j],
                url = paste0(hits$folder[j], hits$file[j]), bytes = if (ok) file.size(dest) else NA_real_,
                sha256 = if (ok) sha256(dest) else NA_character_, stringsAsFactors = FALSE)
            if (ok && is.null(located[[strains$strain[i]]])) located[[strains$strain[i]]] <<- dest
        }
        if (!nrow(hits)) fasta[[length(fasta) + 1L]] <- data.frame(strain = strains$strain[i], project = "none found",
            url = NA_character_, bytes = NA_real_, sha256 = NA_character_)
    }
    write_tsv(files, "A-project-files.tsv")
    write_tsv(do.call(rbind, fasta), "A-strain-fasta.tsv")
})

# --- F4: the GEO annotations; CUT: the B36 cut ---------------------------------------------------------------------
geo <- list()
section("F4", {
    rows <- list()
    attributes <- list()
    for (i in seq_len(nrow(strains))) {
        series <- strains$geo[i]
        folder <- sprintf("https://ftp.ncbi.nlm.nih.gov/geo/series/%snnn/%s/suppl/", substr(series, 1, 6), series)
        files <- listing(folder)
        log_line("F4", series, paste(files, collapse = " "))
        for (file in files[grepl("\\.(gff|gtf)\\.gz$", files)]) {
            dest <- file.path(work, file)
            if (!fetch(paste0(folder, file), dest)) next
            format <- if (grepl("\\.gtf\\.gz$", file)) "gtf" else "gff3"
            geo[[paste(strains$strain[i], format)]] <<- dest
            lines <- readLines(gzfile(dest), warn = FALSE)
            body <- lines[nzchar(lines) & !startsWith(lines, "#")]
            fields <- strsplit(body, "\t", fixed = TRUE)
            types <- table(vapply(fields, `[`, "", 3L))
            attr_text <- vapply(fields, `[`, "", 9L)
            keys <- if (format == "gtf") unlist(regmatches(attr_text, gregexpr("(^|; ?)[A-Za-z_]+ ", attr_text))) else
                unlist(regmatches(attr_text, gregexpr("(^|;)[A-Za-z_]+=", attr_text)))
            keys <- table(gsub("[; =]", "", keys))
            cds <- vapply(fields, `[`, "", 3L) == "CDS"
            rows[[length(rows) + 1L]] <- data.frame(strain = strains$strain[i], series = series, file = file,
                format = format, bytes = file.size(dest), sha256 = sha256(dest), lines = length(body),
                cds_rows = sum(cds), cds_with_protein_id = sum(cds & grepl("protein_id", attr_text)),
                feature_types = paste(sprintf("%s:%d", names(types), as.integer(types)), collapse = ";"),
                stringsAsFactors = FALSE)
            attributes[[length(attributes) + 1L]] <- data.frame(strain = strains$strain[i], file = file,
                attribute = names(keys), rows = as.integer(keys), stringsAsFactors = FALSE)
        }
    }
    write_tsv(do.call(rbind, rows), "F4-geo-annotation.tsv")
    write_tsv(do.call(rbind, attributes), "F4-geo-attributes.tsv")
})
section("CUT", {
    gff <- geo[["B36 gff3"]]
    gtf <- geo[["B36 gtf"]]
    if (is.null(gff)) stop("the B36 GFF3 was not fetched")
    lines <- readLines(gzfile(gff), warn = FALSE)
    body <- which(nzchar(lines) & !startsWith(lines, "#"))
    attr_of <- function(text, key) ifelse(grepl(paste0("(^|;)", key, "="), text),
        sub(paste0("^.*(^|;)", key, "=([^;]*).*$"), "\\2", text), NA_character_)
    fields <- strsplit(lines[body], "\t", fixed = TRUE)
    type <- vapply(fields, `[`, "", 3L)
    locus <- attr_of(vapply(fields, `[`, "", 9L), "locus_tag")
    protein <- attr_of(vapply(fields, `[`, "", 9L), "protein_id")
    cds_loci <- locus[type == "CDS" & !is.na(locus)]
    multi_row <- unique(cds_loci[duplicated(cds_loci)])
    pseudo <- unique(locus[!is.na(locus) & grepl("pseudo=true", lines[body])])
    rna <- unique(locus[type %in% c("rRNA", "tRNA", "ncRNA") & !is.na(locus)])
    shared <- unique(protein[type == "CDS" & !is.na(protein)][duplicated(paste(protein, locus)[type == "CDS" &
        !is.na(protein)]) == FALSE])
    by_protein <- tapply(locus[type == "CDS" & !is.na(protein)], protein[type == "CDS" & !is.na(protein)],
        function(l) length(unique(l)))
    several <- names(by_protein)[by_protein > 1L][1:2]
    several_loci <- unique(locus[type == "CDS" & protein %in% several])
    plain <- unique(cds_loci)[1:40]
    chosen <- unique(stats::na.omit(c(plain, utils::head(multi_row, 2L), utils::head(pseudo, 3L), utils::head(rna, 3L),
        several_loci)))
    keep <- body[!is.na(locus) & locus %in% chosen]
    cut <- c(lines[startsWith(lines, "##gff-version")][1L], lines[keep])
    writeLines(cut, file.path(out, "CUT-b36.gff"))
    cut_gtf <- character()
    if (!is.null(gtf)) {
        gtf_lines <- readLines(gzfile(gtf), warn = FALSE)
        gtf_locus <- ifelse(grepl("locus_tag \"", gtf_lines), sub("^.*locus_tag \"([^\"]*)\".*$", "\\1", gtf_lines), NA)
        gtf_gene <- ifelse(grepl("gene_id \"", gtf_lines), sub("^.*gene_id \"([^\"]*)\".*$", "\\1", gtf_lines), NA)
        cut_gtf <- c(gtf_lines[startsWith(gtf_lines, "#")], gtf_lines[!startsWith(gtf_lines, "#") &
            ((!is.na(gtf_locus) & gtf_locus %in% chosen) | (!is.na(gtf_gene) & gtf_gene %in% chosen))])
        writeLines(cut_gtf, file.path(out, "CUT-b36.gtf"))
    }
    # The expected protein accessions per locus, by awk on the cut file (not the package's parser): CDS rows, the
    # locus_tag and protein_id attributes, distinct protein ids joined with ";" in C-locale order.
    awk <- paste("awk -F'\\t' '$3==\"CDS\" { lt=\"\"; pid=\"\"; n=split($9,a,\";\");",
        "for(i=1;i<=n;i++){ split(a[i],kv,\"=\"); if(kv[1]==\"locus_tag\") lt=kv[2]; if(kv[1]==\"protein_id\") pid=kv[2] }",
        "if(lt!=\"\" && pid!=\"\") print lt \"\\t\" pid }'", shQuote(file.path(out, "CUT-b36.gff")),
        "| LC_ALL=C sort -u")
    pairs <- utils::read.delim(text = system(awk, intern = TRUE), header = FALSE, col.names = c("locus_tag", "protein"),
        colClasses = "character")
    expected <- data.frame(locus_tag = sort(chosen, method = "radix"), stringsAsFactors = FALSE)
    expected$protein_accession <- vapply(expected$locus_tag, function(l) {
        p <- sort(unique(pairs$protein[pairs$locus_tag == l]), method = "radix")
        if (length(p)) paste(p, collapse = ";") else NA_character_
    }, "")
    write_tsv(expected, "CUT-b36-expected.tsv")
    log_line("CUT", length(chosen), "loci;", length(cut), "GFF lines;", length(cut_gtf), "GTF lines")
})

# --- F6: the reference order ------------------------------------------------------------------------------------------
section("F6", {
    rows <- list()
    organisms <- readReferenceOrganisms()
    for (i in seq_len(nrow(organisms))) {
        upid <- organisms$uniprot_upid[i]
        json <- read_text(sprintf("https://rest.uniprot.org/proteomes/%s.json", upid))
        assembly <- regmatches(json, regexpr("\"assemblyId\"[[:space:]]*:[[:space:]]*\"GC[AF]_[0-9]+\\.[0-9]+\"", json))
        assembly <- sub("^.*\"(GC[AF]_[0-9]+\\.[0-9]+)\"$", "\\1", assembly)
        log_line("F6", upid, "assembly", if (length(assembly)) assembly else "none")
        if (!length(assembly)) next
        digits <- sub("^GC[AF]_([0-9]+)\\..*$", "\\1", assembly)
        base <- sprintf("https://ftp.ncbi.nlm.nih.gov/genomes/all/%s/%s/%s/%s/", substr(assembly, 1, 3),
            substr(digits, 1, 3), substr(digits, 4, 6), substr(digits, 7, 9))
        folder <- grep(paste0("^", assembly, "_"), listing(base), value = TRUE)[1L]
        folder <- sub("/$", "", folder)
        gff <- file.path(work, paste0(upid, "_genomic.gff.gz"))
        if (is.na(folder) || !fetch(sprintf("%s%s/%s_genomic.gff.gz", base, folder, folder), gff)) {
            log_line("F6", upid, "genome GFF not fetched")
            next
        }
        lines <- readLines(gzfile(gff), warn = FALSE)
        body <- lines[nzchar(lines) & !startsWith(lines, "#")]
        fields <- strsplit(body, "\t", fixed = TRUE)
        type <- vapply(fields, `[`, "", 3L)
        text <- vapply(fields, `[`, "", 9L)
        attr_of <- function(key) ifelse(grepl(paste0("(^|;)", key, "="), text),
            sub(paste0("^.*(^|;)", key, "=([^;]*).*$"), "\\2", text), NA_character_)
        coding <- type == "gene" & grepl("gene_biotype=protein_coding", text)
        loci <- unique(c(stats::na.omit(attr_of("locus_tag")[coding]), unlist(strsplit(stats::na.omit(
            attr_of("old_locus_tag")[coding]), ","))))
        xref <- parseReferenceCrossReferences(readReferenceCrossReferenceTable(referenceOrganismFiles(upid)$xref),
            organisms$kegg_code[i], organisms$string_species[i], organisms$locus_pattern[i])
        present <- !is.na(xref$reference_locus)
        joined <- vapply(strsplit(ifelse(present, xref$reference_locus, ""), ";", fixed = TRUE),
            function(l) length(l) > 0L && any(l %in% loci), logical(1))
        rows[[length(rows) + 1L]] <- data.frame(upid = upid, assembly = assembly, gff = basename(gff),
            gff_sha256 = sha256(gff), protein_coding_loci = length(loci), entries = nrow(xref),
            with_reference_locus = sum(present), joined = sum(joined), share_joined = sum(joined) / nrow(xref),
            pass = sum(joined) / nrow(xref) >= 0.9, stringsAsFactors = FALSE)
    }
    write_tsv(do.call(rbind, rows), "F6-reference-order.tsv")
})

# --- F3: the authors' accessions to UniParc -------------------------------------------------------------------------
section("F3", {
    pattern <- getUniprotRegexPatterns()$entry
    rows <- list()
    summary <- list()
    for (i in seq_len(nrow(strains))) {
        name <- sprintf("Results%s.xlsx", strains$strain[i])
        found <- pride_folder(strains$dia[i])
        file <- found$files[tolower(found$files) == tolower(name)][1L]
        if (is.na(file)) {
            log_line("F3", strains$strain[i], "no", name, "in", strains$dia[i])
            next
        }
        dest <- file.path(work, file)
        if (!fetch(paste0(found$folder, file), dest)) next
        sheets <- readxl::excel_sheets(dest)
        accessions <- character()
        for (sheet in sheets) {
            table <- readxl::read_excel(dest, sheet = sheet, col_types = "text")
            log_line("F3", strains$strain[i], "sheet", sheet, nrow(table), "rows; columns", paste(names(table), collapse = "|"))
            for (column in names(table)) {
                values <- trimws(sub("^(sp|tr)\\|([^|]+)\\|.*$", "\\2", stats::na.omit(table[[column]])))
                if (length(values) && mean(grepl(pattern, values)) > 0.9 && !any(grepl("^WP_", values))) {
                    accessions <- c(accessions, values[grepl(pattern, values)])
                    log_line("F3", strains$strain[i], "accession column", sheet, column)
                }
            }
        }
        accessions <- sort(unique(sub("-[0-9]+$", "", accessions)), method = "radix")
        upi <- character(length(accessions))
        status <- character(length(accessions))
        release <- NA_character_
        for (k in seq_along(accessions)) {
            response <- tryCatch(httr2::request(sprintf("https://rest.uniprot.org/uniprotkb/%s.json", accessions[k])) |>
                httr2::req_retry(max_tries = 3L) |> httr2::req_error(is_error = function(r) FALSE) |>
                httr2::req_perform(), error = function(e) NULL)
            if (is.null(response) || httr2::resp_status(response) != 200L) {
                status[k] <- "not_resolvable"
                upi[k] <- NA_character_
                next
            }
            release <- httr2::resp_header(response, "x-uniprot-release") %||% release
            body <- httr2::resp_body_string(response)
            status[k] <- if (grepl("\"entryType\"[[:space:]]*:[[:space:]]*\"Inactive\"", body)) "inactive" else "active"
            hit <- regmatches(body, regexpr("UPI[0-9A-F]{10}", body))
            upi[k] <- if (length(hit)) hit else NA_character_
        }
        rows[[length(rows) + 1L]] <- data.frame(strain = strains$strain[i], accession = accessions,
            uniprotkb_status = status, upi = upi, release = release, stringsAsFactors = FALSE)
        summary[[length(summary) + 1L]] <- data.frame(strain = strains$strain[i], accessions = length(accessions),
            resolved = sum(!is.na(upi)), share = if (length(upi)) mean(!is.na(upi)) else NA_real_,
            evaluable = length(upi) > 0 && mean(!is.na(upi)) >= 0.9, release = release)
    }
    if (length(rows)) write_tsv(do.call(rbind, rows), "F3-authors-accession-upi.tsv")
    if (length(summary)) write_tsv(do.call(rbind, summary), "F3-summary.tsv")
})

# --- DEMO: the demo tables and the real-engine lane ----------------------------------------------------------------
section("DEMO", {
    archive <- file.path(work, "diamond-linux64.tar.gz")
    pin <- readDiamondPin()
    linux <- pin[pin$platform == "linux", ]
    if (!fetch(linux$archive_url, archive) || sha256(archive) != linux$archive_sha256) stop("engine archive not verified")
    utils::untar(archive, exdir = file.path(work, "engine"))
    Sys.setenv(MSR_DIAMOND = file.path(work, "engine", "diamond"))
    Sys.chmod(Sys.getenv("MSR_DIAMOND"), "0755")
    rows <- list()
    for (i in seq_len(nrow(strains))) {
        fasta <- located[[strains$strain[i]]]
        if (is.null(fasta)) {
            log_line("DEMO", strains$strain[i], "no strain FASTA located")
            next
        }
        started <- Sys.time()
        result <- mapStrainToReference(fasta, strains$reference[i], threads = 4L)
        seconds <- as.numeric(difftime(Sys.time(), started, units = "secs"))
        write_tsv(result$mapping, sprintf("DEMO-%s.tsv", strains$strain[i]))
        decision <- result$mapping[, c("strain_accession", "reference_accession", "map_status", "reference_accession_used")]
        rows[[length(rows) + 1L]] <- data.frame(strain = strains$strain[i], reference = strains$reference[i],
            threads = 4L, seconds = round(seconds, 1), strain_proteins = result$counts$strain_proteins,
            t(unlist(result$counts$by_status)), identity_median = result$identity_summary$median,
            identity_p05 = result$identity_summary$p05, mapping_sha256 = result$provenance$mapping_sha256,
            decision_sha256 = canonicalReferenceMappingSha256(decision), check.names = FALSE)
        if (identical(strains$strain[i], "KPC2")) {
            again <- mapStrainToReference(fasta, strains$reference[i], threads = 1L)
            decision_again <- again$mapping[, names(decision)]
            rows[[length(rows) + 1L]] <- data.frame(strain = "KPC2", reference = strains$reference[i], threads = 1L,
                seconds = NA_real_, strain_proteins = again$counts$strain_proteins, t(unlist(again$counts$by_status)),
                identity_median = again$identity_summary$median, identity_p05 = again$identity_summary$p05,
                mapping_sha256 = again$provenance$mapping_sha256,
                decision_sha256 = canonicalReferenceMappingSha256(decision_again), check.names = FALSE)
        }
    }
    write_tsv(do.call(rbind, rows), "DEMO-summary.tsv")
})
