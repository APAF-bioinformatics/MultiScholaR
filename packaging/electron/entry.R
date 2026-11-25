# MultiScholaR bundle (KIT-134 P3): the R entry the shell starts, with the bundle's R:
#   Rscript --vanilla entry.R <port> <run folder> <shell pid>
#
# When the installed MultiScholaR has runPackagedApplication (the package side of P3, made after
# KIT-131 is merged), this file only calls it: the start-up probe's choice, the plain launcher and
# the systemd worker are the package's. Until then it starts the application in this R process
# (revision 2.0, 2.0.5, the third step of the fallback order) and says so on the MSR-LAUNCHER line.
# Either way R must not outlive the shell: a watchdog ends R when the shell's process is gone.
# Never calls loadDependencies(): nothing is installed or downloaded at start (F12).

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 3L) stop("usage: Rscript --vanilla entry.R <port> <run folder> <shell pid>")

# KIT-156: the test setting for data testing. MULTISCHOLAR_TEST_FEATURES=1 in the environment the application is started
# from switches on, for this session only, the pre-release feature that is off by default: the RUV design panel (option
# MultiScholaR.ruv_design_panel). KIT-240: the integration module is on by the package's own default
# (MultiScholaR.integration_gui), so every launcher shows it and this setting no longer touches it. The package's
# defaults are unchanged; without the variable this file does exactly what it did before. The option lives in this R
# process, so the application must run here: where the shell's probe chose the systemd worker (whose controller is a
# separate unit that would not see it), the plain launcher is used instead, and the reason says so.
if (identical(Sys.getenv("MULTISCHOLAR_TEST_FEATURES"), "1")) {
    options(MultiScholaR.ruv_design_panel = TRUE)
    cat(sprintf("MSR-TEST-FEATURES %s\n", jsonlite::toJSON(list(variable = "MULTISCHOLAR_TEST_FEATURES",
        ruv_design_panel = getOption("MultiScholaR.ruv_design_panel")), auto_unbox = TRUE)))
    probe <- tryCatch(jsonlite::fromJSON(Sys.getenv("MSR_LAUNCHER_PROBE"), simplifyVector = FALSE), error = function(error) NULL)
    if (is.list(probe) && identical(probe$launcher, "systemd")) {
        Sys.setenv(MSR_LAUNCHER_PROBE = jsonlite::toJSON(list(launcher = "process", reason = paste(
            "the test setting MULTISCHOLAR_TEST_FEATURES keeps the application in this R process; the systemd worker",
            "was available:", paste(as.character(probe$reason %||% ""), collapse = " "))), auto_unbox = TRUE))
    }
    flush(stdout())
}
namespace <- asNamespace("MultiScholaR")
if (exists("runPackagedApplication", envir = namespace, mode = "function", inherits = FALSE)) {
    get("runPackagedApplication", envir = namespace)(arguments)
    quit(save = "no", status = 0L)
}

port <- as.integer(arguments[[1L]])
run_folder <- arguments[[2L]]
shell_pid <- as.integer(arguments[[3L]])
reason <- paste("this MultiScholaR has no runPackagedApplication yet (the package side of KIT-134 P3",
    "waits for KIT-131), so steps run inside the application")
cat(sprintf("MSR-LAUNCHER %s\n", jsonlite::toJSON(list(launcher = "in-process", reason = reason,
    probe = jsonlite::fromJSON(Sys.getenv("MSR_LAUNCHER_PROBE", "null"))), auto_unbox = TRUE, null = "null")))
flush(stdout())

shellAlive <- function() {
    tryCatch(ps::ps_is_running(ps::ps_handle(shell_pid)), error = function(e) FALSE)
}
watch <- function() {
    if (!shellAlive()) {
        cat("the shell is gone: stopping R\n")
        quit(save = "no", status = 0L)
    }
    later::later(watch, 2)
}
later::later(watch, 2)
setwd(run_folder)
shiny::runApp(MultiScholaR::run_app(), host = "127.0.0.1", port = port, launch.browser = FALSE)
