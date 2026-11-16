# KIT-130 piece 1: R 4.6.1 runs from wherever its folder is, with no other R on the PATH,
# --vanilla, and library paths inside the folder only.
#
# Usage: <bundle>/.../Rscript --vanilla check-r.R <bundle directory>

source(file.path(dirname(sub("^--file=", "", grep("^--file=", commandArgs(), value = TRUE))), "common.R"))
bundle <- expectedBundle()

cat(sprintf("R: %s\nR home: %s\n", R.version.string, R.home()))
cat(sprintf("library paths: %s\n", paste(.libPaths(), collapse = " | ")))
cat(sprintf("command line: %s\n", paste(commandArgs(), collapse = " ")))
cat(sprintf("R or Rscript on the PATH: '%s'\n", paste(Sys.which(c("R", "Rscript")), collapse = "' '")))
cat(sprintf("PATH: %s\n", Sys.getenv("PATH")))
print(Sys.getenv(c("R_HOME", "R_LIBS", "R_LIBS_USER", "R_LIBS_SITE", "R_ENVIRON_USER", "R_PROFILE_USER")))
print(extSoftVersion()[c("BLAS", "PCRE", "ICU", "zlib")])
cat(sprintf("LAPACK: %s\n", La_library()))

stopifnot(
    "R is 4.6.1" = getRversion() == "4.6.1",
    "R home is inside the bundle" = insideBundle(R.home(), bundle),
    "every library path is inside the bundle" = all(insideBundle(.libPaths(), bundle)),
    "no other R on the PATH" = length(foreignOnPath(c("R", "Rscript"), bundle)) == 0L,
    "started with --vanilla" = "--vanilla" %in% commandArgs()
)

# The base system works, not only starts: linear algebra, a graphics device, compression.
x <- matrix(c(4, 1, 1, 3), 2L)
stopifnot(isTRUE(all.equal(solve(x) %*% x, diag(2L))))
plot_file <- tempfile(fileext = ".pdf")
grDevices::pdf(plot_file)
graphics::plot(1:10)
grDevices::dev.off()
stopifnot(file.size(plot_file) > 0, identical(memDecompress(memCompress(charToRaw("kit-130"), "gzip"), "gzip"), charToRaw("kit-130")))

checkPassed("R 4.6.1 in a folder")
