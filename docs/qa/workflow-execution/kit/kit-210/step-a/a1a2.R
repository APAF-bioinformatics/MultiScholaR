# KIT-210 step A, tasks A1 and A2, in one runner job (evidence code under docs/qa, not package code).
#
# Usage (from the package root, after the package commit is installed):
#   Rscript docs/qa/workflow-execution/kit/kit-210/step-a/a1a2.R <out dir>
#
# A1 (KIT-210 4.1): fetch the pinned RUVIIIC 1.0.19 source twice, the CRAN archive RUVIIIC_1.0.19.tar.gz and the
#    GitHub mirror cran/RUVIIIC at c534151f (the commit packaging/multischolar.lock.json locks); SHA-256 of every
#    file in each; the two compared file by file and with the hashes recorded in KIT-210 2.2; the installed
#    package's DESCRIPTION (Version, RemoteSha). Writes source-hashes.tsv, tarballs.tsv, md5-check.tsv,
#    installed.tsv, functions.tsv.
# A2 (KIT-210 4.1): the installed compiled library zero-fills, and only in the direction step. A 12 x 60
#    log-scale matrix (one target, 10 potential controls observed everywhere, 49 features that are neither),
#    k = 2, three replicate sets of four runs. (a) one NA of a feature that is neither target nor control, on a row
#    where the target is observed, replaced by 0: the output is identical(); (b) replaced by 20: the output
#    differs by more than 1e-6 (the negative control of (a)); (c) with withExtra, withW, withAlpha: (1) W and
#    alpha(target) reproduce the corrected values to 1e-10; (2) W and the controls' alpha are reproduced from U,
#    the top-k eigenvectors of P Sel (Y0 Y0') Sel' P computed here from the zero-filled matrix, to 1e-8 after
#    aligning the basis; (3) the replacement by 20 changes U and the output, by 0 changes nothing. Writes a2.tsv.
#
# Network: only the two archives named above. Nothing else is downloaded or called.
# Every check is wrapped so that an error is recorded as "error: <message>" (an instrument fault, no result)
# rather than as a fail.

args <- commandArgs(trailingOnly = TRUE)
out <- if (length(args)) args[[1L]] else "kit210-step-a"
dir.create(out, recursive = TRUE, showWarnings = FALSE)
options(warn = 1L, stringsAsFactors = FALSE)

writeTsv <- function(x, name) {
  utils::write.table(x, file.path(out, name), sep = "\t", quote = FALSE, row.names = FALSE, na = "NA")
}
sha256 <- function(paths) {
  if (!length(paths)) return(character())
  res <- system2("sha256sum", shQuote(paths), stdout = TRUE)
  sub(" .*$", "", res)
}
fmt <- function(x) formatC(x, format = "e", digits = 3)

SHA <- "c534151f6b2265b36294ff630f4dd13c2e86a9fa"
work <- file.path(tempdir(), "k210-src")
dir.create(work, showWarnings = FALSE)

## ---- A1: fetch -------------------------------------------------------------------------------------------------
candidates <- list(
  cran = c("https://cran.r-project.org/src/contrib/RUVIIIC_1.0.19.tar.gz",
    "https://cran.r-project.org/src/contrib/Archive/RUVIIIC/RUVIIIC_1.0.19.tar.gz"),
  mirror = sprintf("https://github.com/cran/RUVIIIC/archive/%s.tar.gz", SHA))
fetchOne <- function(urls, dest) {
  for (u in urls) {
    unlink(dest)
    ok <- tryCatch({
      utils::download.file(u, dest, mode = "wb", quiet = TRUE)
      TRUE
    }, error = function(e) FALSE)
    if (isTRUE(ok) && file.exists(dest) && file.size(dest) > 0) return(u)
  }
  NA_character_
}
tarballs <- data.frame(source = names(candidates), url = NA_character_, sha256 = NA_character_, bytes = NA_real_,
  recorded_sha256 = c(NA_character_, "846cc1442ad90711efade03c9e2f6e90c963753023860a91efe9ef63958d4009"))
roots <- list()
for (i in seq_along(candidates)) {
  src <- names(candidates)[i]
  dest <- file.path(work, paste0(src, ".tar.gz"))
  used <- fetchOne(candidates[[i]], dest)
  tarballs$url[i] <- used
  if (!is.na(used)) {
    tarballs$sha256[i] <- sha256(dest)
    tarballs$bytes[i] <- file.size(dest)
    exdir <- file.path(work, src)
    utils::untar(dest, exdir = exdir)
    top <- list.dirs(exdir, recursive = FALSE)
    roots[[src]] <- if (length(top) == 1L) top else exdir
  }
}
writeTsv(tarballs, "tarballs.tsv")

listHashes <- function(root) {
  files <- sort(list.files(root, recursive = TRUE, all.files = TRUE, no.. = TRUE))
  data.frame(file = files, sha256 = sha256(file.path(root, files)))
}
hashes <- lapply(roots, listHashes)
cran <- if (!is.null(hashes$cran)) hashes$cran else data.frame(file = character(), sha256 = character())
mirror <- if (!is.null(hashes$mirror)) hashes$mirror else data.frame(file = character(), sha256 = character())
names(cran)[2] <- "cran_sha256"
names(mirror)[2] <- "mirror_sha256"
src_tab <- merge(cran, mirror, by = "file", all = TRUE)
src_tab$archives_equal <- !is.na(src_tab$cran_sha256) & !is.na(src_tab$mirror_sha256) &
  src_tab$cran_sha256 == src_tab$mirror_sha256

# KIT-210 2.2 (r1.1/r1.2): prefix ... suffix as written there (src/RUVIIIC.cpp in full).
ticket <- data.frame(
  file = c("src/RUVIIIC_Varying.cpp", "R/RUVIII_C_Varying.R", "tests/testthat/test-compare_RUVIII_Varying.R",
    "tests/testthat/test-RUVIII_C_Varying_CPP_implementation.R", "src/RUVIIIC.cpp", "R/RUVIII_C.R"),
  prefix = c("a6d8c74d", "6809fe96", "36aa0989", "413ab716",
    "57615f70c68ba1357a6543e9533e5714fad855e1ecb78cba0acd5bc48c0e1797", "2d4859d9"),
  suffix = c("29fc743", "3e9314ab", "34722e", "f668ee", "", "3183e253"))
src_tab$ticket_2_2 <- NA_character_
src_tab$ticket_match <- NA
for (j in seq_len(nrow(ticket))) {
  r <- which(src_tab$file == ticket$file[j])
  rec <- if (nzchar(ticket$suffix[j])) paste0(ticket$prefix[j], "...", ticket$suffix[j]) else ticket$prefix[j]
  if (!length(r)) {
    src_tab <- rbind(src_tab, data.frame(file = ticket$file[j], cran_sha256 = NA, mirror_sha256 = NA,
      archives_equal = FALSE, ticket_2_2 = rec, ticket_match = FALSE))
    next
  }
  h <- c(src_tab$cran_sha256[r], src_tab$mirror_sha256[r])
  src_tab$ticket_2_2[r] <- rec
  src_tab$ticket_match[r] <- all(!is.na(h) & startsWith(h, ticket$prefix[j]) & endsWith(h, ticket$suffix[j]))
}
writeTsv(src_tab, "source-hashes.tsv")

# The archives' own MD5 lists against their files (reported).
md5rows <- list()
for (src in names(roots)) {
  f <- file.path(roots[[src]], "MD5")
  if (!file.exists(f)) {
    md5rows[[src]] <- data.frame(source = src, listed = NA, matching = NA, note = "no MD5 file")
    next
  }
  l <- readLines(f)
  l <- l[nzchar(l)]
  listed <- sub("^[0-9a-f]+ \\*", "", l)
  want <- sub(" .*$", "", l)
  got <- unname(tools::md5sum(file.path(roots[[src]], listed)))
  md5rows[[src]] <- data.frame(source = src, listed = length(l), matching = sum(!is.na(got) & got == want),
    note = paste(listed[is.na(got) | got != want], collapse = ","))
}
writeTsv(do.call(rbind, md5rows), "md5-check.tsv")

# The installed package.
desc <- utils::packageDescription("RUVIIIC")
fields <- c("Package", "Version", "RemoteType", "RemoteHost", "RemoteUsername", "RemoteRepo", "RemoteSha",
  "RemoteRef", "Packaged", "Date/Publication", "Repository", "Built")
installed <- data.frame(field = fields, value = vapply(fields, function(f) {
  v <- desc[[f]]
  if (is.null(v)) NA_character_ else paste(as.character(v), collapse = " ")
}, ""))
so <- file.path(find.package("RUVIIIC"), "libs", "RUVIIIC.so")
installed <- rbind(installed,
  data.frame(field = "libs/RUVIIIC.so sha256", value = if (file.exists(so)) sha256(so) else NA_character_),
  data.frame(field = "R", value = R.version.string))
writeTsv(installed, "installed.tsv")

# The installed R functions against the pinned source's (reported): formals and deparsed bodies identical.
fn_rows <- list()
if (!is.null(roots$mirror)) {
  srcEnv <- new.env()
  for (f in c("R/RUVIII_C.R", "R/RUVIII_C_Varying.R")) {
    sys.source(file.path(roots$mirror, f), envir = srcEnv, keep.source = FALSE)
  }
  ns <- asNamespace("RUVIIIC")
  for (fn in c("RUVIII_C", "RUVIIIC_R", "RUVIII_C_Varying", "RUVIIIC_Varying_R", "RUVIII_C_Varying_residual_dimension")) {
    same <- tryCatch({
      a <- get(fn, envir = ns)
      b <- get(fn, envir = srcEnv)
      identical(formals(a), formals(b)) && identical(deparse(body(a)), deparse(body(b)))
    }, error = function(e) paste("error:", conditionMessage(e)))
    fn_rows[[fn]] <- data.frame(function_name = fn, installed_equals_source = as.character(same))
  }
}
writeTsv(if (length(fn_rows)) do.call(rbind, fn_rows) else data.frame(function_name = NA, installed_equals_source = NA),
  "functions.tsv")

## ---- A2 --------------------------------------------------------------------------------------------------------
invisible(RUVIIIC::omp_set_num_threads(1L))
set.seed(210000L, kind = "Mersenne-Twister", normal.kind = "Inversion", sample.kind = "Rejection")
nr <- 12L
nf <- 60L
k <- 2L
runs <- sprintf("r%02d", seq_len(nr))
feats <- sprintf("f%02d", seq_len(nf))
target <- "f01"
ctl <- sprintf("f%02d", 2:11)
sets <- factor(rep(c("A", "B", "C"), each = 4L))
M <- stats::model.matrix(~ sets - 1)
dimnames(M) <- list(runs, levels(sets))
# Unwanted factors orthogonal to the replicate sets; a wanted set effect on the target and f40-f60.
w_true <- matrix(stats::rnorm(nr * k, sd = 1.5), nr, k)
w_true <- w_true - M %*% solve(crossprod(M), crossprod(M, w_true))
a_true <- matrix(stats::rnorm(k * nf, sd = 1), k, nf)
mu <- stats::rnorm(nf, 20, 2)
effect <- c(A = 0, B = 1, C = -1)[as.character(sets)]
bio <- outer(effect, as.numeric(feats %in% c(target, sprintf("f%02d", 40:60))))
Y <- matrix(mu, nr, nf, byrow = TRUE) + bio + w_true %*% a_true + matrix(stats::rnorm(nr * nf, sd = 0.2), nr, nf)
dimnames(Y) <- list(runs, feats)
# Missing values: the target at r03; six in features that are neither (one of them, r05 / f20, is the one replaced).
Y["r03", target] <- NA
for (cell in list(c("r05", "f20"), c("r07", "f25"), c("r10", "f33"), c("r02", "f41"), c("r11", "f52"), c("r03", "f30"))) {
  Y[cell[1], cell[2]] <- NA
}
Y_zero <- Y
Y_zero["r05", "f20"] <- 0
Y_twenty <- Y
Y_twenty["r05", "f20"] <- 20

callVarying <- function(y, extra = FALSE) {
  RUVIIIC::RUVIII_C_Varying(k = k, Y = y, M = M, toCorrect = target, potentialControls = ctl,
    withExtra = extra, withW = extra, withAlpha = extra, progress = FALSE)
}
a2 <- list()
add <- function(part, check, observed, rule, pass, deciding = "yes") {
  a2[[length(a2) + 1L]] <<- data.frame(part = part, check = check, observed = as.character(observed), rule = rule,
    pass = as.character(pass), deciding = deciding)
}
guard <- function(part, check, rule, deciding, expr) {
  tryCatch(expr, error = function(e) add(part, check, paste("error:", conditionMessage(e)), rule, NA, deciding))
}

S <- which(!is.na(Y[, target]))
fill0 <- function(y) {
  y[is.na(y)] <- 0
  y
}
# U from the zero-filled matrix restricted to the target's rows, as .cpp:241 forms it (eigen in place of Spectra).
directionsFrom <- function(y0) {
  ys <- y0[S, , drop = FALSE]
  ms <- M[S, , drop = FALSE]
  ms <- ms[, colSums(ms) > 0, drop = FALSE]
  p <- diag(length(S)) - ms %*% solve(crossprod(ms)) %*% t(ms)
  g <- p %*% tcrossprod(ys) %*% t(p)
  g <- (g + t(g)) / 2
  e <- eigen(g, symmetric = TRUE)
  list(u = e$vectors[, seq_len(k), drop = FALSE], values = e$values, ys = ys)
}
# The library's U recovered from its alpha = U' Y0[S, ] (.cpp:270): U' = alpha Y0S' (Y0S Y0S')^-1.
libraryU <- function(alpha, ys) t(alpha %*% t(ys) %*% solve(tcrossprod(ys)))
projDist <- function(u1, u2) max(abs(u1 %*% t(u1) - u2 %*% t(u2)))

out_na <- callVarying(Y)

guard("a", "NA -> 0 in f20 at r05 (neither target nor control; target observed there): output identical()",
  "identical", "yes", {
    same <- identical(out_na, callVarying(Y_zero))
    add("a", "NA -> 0 in f20 at r05 (neither target nor control; target observed there): output identical()",
      same, "identical", same)
  })
guard("b", "NA -> 20 in f20 at r05: largest output difference (negative control of (a))", "> 1e-6", "yes", {
  d <- max(abs(out_na - callVarying(Y_twenty)), na.rm = TRUE)
  add("b", "NA -> 20 in f20 at r05: largest output difference (negative control of (a))", fmt(d), "> 1e-6", d > 1e-6)
})

ext <- callVarying(Y, extra = TRUE)
guard("c1", "y - W alpha(target) against the returned corrected values, target's rows; NA elsewhere", "<= 1e-10", "yes", {
  w_lib <- ext$W[[target]][S, , drop = FALSE]
  alpha_lib <- ext$alpha[[target]]
  rec <- Y[S, target] - drop(w_lib %*% alpha_lib[, target])
  d <- max(abs(rec - ext$newY[S, target]))
  na_ok <- all(is.na(ext$newY[-S, target]))
  add("c1", "y - W alpha(target) against the returned corrected values, target's rows; NA elsewhere",
    paste0(fmt(d), "; NA rows ", na_ok, "; residual dimension ", ext$residualDimensions[[target]]),
    "<= 1e-10", d <= 1e-10 && na_ok)
})

guard("c2", "W and the controls' alpha from U of P Sel (Y0 Y0') Sel' P (zero-filled), basis aligned", "<= 1e-8", "yes", {
  dirs <- directionsFrom(fill0(Y))
  u <- dirs$u
  ys <- dirs$ys
  ctl_s <- ctl[colSums(is.na(Y[S, ctl, drop = FALSE])) == 0L]
  yc <- ys[, ctl_s, drop = FALSE]
  ac <- t(u) %*% yc
  w_r <- yc %*% t(ac) %*% solve(ac %*% t(ac))
  alpha_lib <- ext$alpha[[target]]
  w_lib <- ext$W[[target]][S, , drop = FALSE]
  u_lib <- libraryU(alpha_lib, ys)
  rot <- t(u) %*% u_lib
  d_rot <- max(abs(t(rot) %*% rot - diag(k)))
  d_w <- max(abs(w_lib - w_r %*% rot))
  d_a <- max(abs(alpha_lib[, ctl_s, drop = FALSE] - t(rot) %*% ac))
  d_free <- max(abs(drop(w_lib %*% alpha_lib[, target]) - drop(w_r %*% (t(u) %*% ys[, target]))))
  gap <- (dirs$values[k] - dirs$values[k + 1L]) / dirs$values[1L]
  worst <- max(d_rot, d_w, d_a, d_free)
  add("c2", "W and the controls' alpha from U of P Sel (Y0 Y0') Sel' P (zero-filled), basis aligned",
    paste0("W ", fmt(d_w), "; alpha(controls) ", fmt(d_a), "; W alpha(target) ", fmt(d_free),
      "; R'R - I ", fmt(d_rot), "; controls used ", length(ctl_s), "; (lambda_k - lambda_k+1) / lambda_1 ", fmt(gap)),
    "<= 1e-8", worst <= 1e-8)
  # Negative control (reported): the same reproduction from the matrix with every missing value filled by 20.
  y20 <- Y
  y20[is.na(y20)] <- 20
  d20 <- directionsFrom(y20)
  ac20 <- t(d20$u) %*% (d20$ys[, ctl_s, drop = FALSE])
  w20 <- d20$ys[, ctl_s, drop = FALSE] %*% t(ac20) %*% solve(ac20 %*% t(ac20))
  free20 <- max(abs(drop(w_lib %*% alpha_lib[, target]) - drop(w20 %*% (t(d20$u) %*% d20$ys[, target]))))
  add("c2-neg", "the same reproduction with every missing value filled by 20 instead of 0: W alpha(target) difference",
    fmt(free20), "> 1e-6", free20 > 1e-6, "reported")
})

guard("c3", "NA -> 0 in f20 at r05: W, alpha, residual dimension and output identical(); NA -> 20: U and output change",
  "0: identical; 20: > 1e-6", "yes", {
    same0 <- identical(ext, callVarying(Y_zero, extra = TRUE))
    ext20 <- callVarying(Y_twenty, extra = TRUE)
    u_lib <- libraryU(ext$alpha[[target]], fill0(Y)[S, , drop = FALSE])
    u_lib20 <- libraryU(ext20$alpha[[target]], fill0(Y_twenty)[S, , drop = FALSE])
    d_u <- projDist(u_lib, u_lib20)
    d_y <- max(abs(ext$newY - ext20$newY), na.rm = TRUE)
    add("c3", "NA -> 0 in f20 at r05: W, alpha, residual dimension and output identical(); NA -> 20: U and output change",
      paste0("0: identical ", same0, "; 20: U projector ", fmt(d_u), ", output ", fmt(d_y)),
      "0: identical; 20: > 1e-6", same0 && d_u > 1e-6 && d_y > 1e-6)
  })

# Reported: the fixed-control function RUVIII_C zero-fills the same way (S2; src/RUVIIIC.cpp:88-93).
guard("d", "RUVIII_C (fixed controls): NA -> 0 identical(); NA -> 20 differs", "0: identical; 20: > 1e-6", "reported", {
  callFixed <- function(y) RUVIIIC::RUVIII_C(k = k, Y = y, M = M, toCorrect = target, controls = ctl, progress = FALSE)
  base <- callFixed(Y)
  same0 <- identical(base, callFixed(Y_zero))
  d20 <- max(abs(base - callFixed(Y_twenty)), na.rm = TRUE)
  add("d", "RUVIII_C (fixed controls): NA -> 0 identical(); NA -> 20 differs",
    paste0("0: identical ", same0, "; 20: ", fmt(d20)), "0: identical; 20: > 1e-6", same0 && d20 > 1e-6, "reported")
})

a2_tab <- do.call(rbind, a2)
writeTsv(a2_tab, "a2.tsv")

## ---- verdict lines ---------------------------------------------------------------------------------------------
ga1 <- nrow(src_tab) > 0 && all(src_tab$archives_equal) && all(src_tab$ticket_match[!is.na(src_tab$ticket_2_2)]) &&
  isTRUE(unname(installed$value[installed$field == "Version"]) == "1.0.19") &&
  isTRUE(unname(installed$value[installed$field == "RemoteSha"]) == SHA)
dec <- a2_tab[a2_tab$deciding == "yes", ]
ga2 <- if (any(is.na(dec$pass) | dec$pass == "NA")) "instrument error" else if (all(dec$pass == "TRUE")) "pass" else "FAIL"
summary_lines <- c(
  paste("files compared:", nrow(src_tab), "; equal in both archives:", sum(src_tab$archives_equal)),
  paste("GA1:", if (ga1) "pass" else "FAIL"),
  paste("GA2:", ga2),
  paste("deciding A2 parts:", paste(paste(dec$part, dec$pass), collapse = ", ")))
writeLines(summary_lines, file.path(out, "summary.txt"))
writeLines(utils::capture.output(utils::sessionInfo()), file.path(out, "session.txt"))
cat(summary_lines, sep = "\n")
