# KIT-245 G0 probe (ticket r1.2 section 8): what "Push to GitHub" (pushProjectToGithub) does today under the
# installed bundle's R, with no real token, no GitHub call and no network path to github.com.
# Run with the installed Rscript --vanilla from a folder outside the checkout. Base code is never changed.
#
# Usage: Rscript --vanilla g0-probe.R <output folder> <host label> <main|p2n> <remote root>
#
# Guards (the calling step sets them; step 0 asserts them and stops without the END line if one fails):
#   the four token variables empty; GIT_ALLOW_PROTOCOL=file (unset for p2n); GIT_TERMINAL_PROMPT=0;
#   GIT_CONFIG_NOSYSTEM=1; GIT_CONFIG_GLOBAL = K245_GITCONFIG_REDIRECT (main) or K245_GITCONFIG_PLAIN (p2n);
#   GIT_CONFIG_COUNT and GIT_CONFIG_PARAMETERS unset; HOME = K245_EXPECT_HOME (no git configuration in it);
#   no environment value shaped like a GitHub token; gh::gh replaced by a recorder that never calls the original
#   and stops on any token other than the dummy.
# The dummy is the plainly fake string below; every line printed by this script has it replaced by <dummy>.
# Lines: "KIT245<TAB>G0<TAB><host><TAB><key><TAB><value>"; last line "G0 END <host>".

args <- commandArgs(trailingOnly = TRUE)
stopifnot(length(args) == 4L)
outDir <- args[[1L]]
host <- args[[2L]]
caseSet <- args[[3L]]
remoteRoot <- args[[4L]]
stopifnot(caseSet %in% c("main", "p2n"))
dir.create(outDir, recursive = TRUE, showWarnings = FALSE)
dir.create(remoteRoot, recursive = TRUE, showWarnings = FALSE)

DUMMY <- "kit245-dummy-not-a-token"
ORG <- "kit245-probe-no-such-org"
isWindows <- identical(.Platform$OS.type, "windows")

mask <- function(x) gsub(DUMMY, "<dummy>", x, fixed = TRUE)
flat <- function(x) gsub("[\t\r\n]+", " | ", paste(mask(as.character(x)), collapse = " ; "))
say <- function(key, value) cat("KIT245", "G0", host, key, paste0(flat(value), "\n"), sep = "\t")
yn <- function(x) if (isTRUE(x)) "yes" else "no"
fault <- function(what) {
    say("INSTRUMENT_FAULT", what)
    cat("KIT245 probe stopped before its END line (instrument fault, VAL-54)\n")
    quit(save = "no", status = 3L)
}

gitPath <- Sys.which("git")
hasGit <- nzchar(gitPath)
gitText <- function(args) {
    # system2 runs through a shell: every argument is quoted (G0 run 45587834788: the unquoted regex was a shell syntax error).
    out <- suppressWarnings(system2("git", shQuote(args), stdout = TRUE, stderr = TRUE))
    list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), text = out)
}

# ---- (0) Step 0: the guards, asserted before anything runs git or gh.
say("step0.case_set", caseSet)
for (v in c("GITHUB_PAT", "GITHUB_TOKEN", "GH_TOKEN", "GITHUB_PAT_GITHUB_COM")) {
    set <- nzchar(Sys.getenv(v))
    say(paste0("step0.", v), if (set) "set" else "unset or empty")
    if (set) fault(paste(v, "is set"))
}
envValues <- Sys.getenv()
shaped <- names(envValues)[grepl("^(gh[pousr]_|github_pat_)", envValues)]
say("step0.token_shaped_variables", if (length(shaped)) shaped else "<none>")
if (length(shaped)) fault(paste("token-shaped value in", paste(shaped, collapse = ", ")))
for (v in c("GIT_CONFIG_COUNT", "GIT_CONFIG_PARAMETERS")) {
    set <- !is.na(Sys.getenv(v, unset = NA))
    say(paste0("step0.", v), if (set) "SET" else "unset")
    if (set) fault(paste(v, "is set"))
}
expectAllow <- if (caseSet == "main") "file" else NA_character_
allow <- Sys.getenv("GIT_ALLOW_PROTOCOL", unset = NA)
say("step0.GIT_ALLOW_PROTOCOL", if (is.na(allow)) "unset" else allow)
if (!identical(allow, expectAllow)) fault("GIT_ALLOW_PROTOCOL not as the case set requires")
if (!identical(Sys.getenv("GIT_TERMINAL_PROMPT"), "0")) fault("GIT_TERMINAL_PROMPT is not 0")
if (!identical(Sys.getenv("GIT_CONFIG_NOSYSTEM"), "1")) fault("GIT_CONFIG_NOSYSTEM is not 1")
say("step0.GIT_TERMINAL_PROMPT", Sys.getenv("GIT_TERMINAL_PROMPT"))
say("step0.GIT_CONFIG_NOSYSTEM", Sys.getenv("GIT_CONFIG_NOSYSTEM"))
redirectFile <- Sys.getenv("K245_GITCONFIG_REDIRECT")
plainFile <- Sys.getenv("K245_GITCONFIG_PLAIN")
inForce <- Sys.getenv("GIT_CONFIG_GLOBAL")
wanted <- if (caseSet == "main") redirectFile else plainFile
say("step0.GIT_CONFIG_GLOBAL", inForce)
say("step0.config_in_force", if (identical(inForce, redirectFile)) "gitconfig-redirect" else
    if (identical(inForce, plainFile)) "gitconfig-plain" else "OTHER")
if (!nzchar(wanted) || !identical(normalizePath(inForce, mustWork = FALSE), normalizePath(wanted, mustWork = FALSE)) ||
    !file.exists(inForce)) {
    fault("GIT_CONFIG_GLOBAL is not the file the case set requires")
}
configLines <- readLines(inForce, warn = FALSE)
say("step0.config_file_content", if (length(configLines)) configLines else "<empty>")
hasRedirect <- any(grepl(paste0("insteadOf = https://x-access-token:", DUMMY, "@github.com/"), configLines, fixed = TRUE))
hasProxy <- any(grepl("proxy = http://127.0.0.1:9", configLines, fixed = TRUE))
if (caseSet == "main" && !(hasRedirect && hasProxy)) fault("gitconfig-redirect lacks the redirect or the dead proxy")
if (caseSet == "p2n" && (hasRedirect || hasProxy || length(grep("\\S", configLines)) > 0L)) fault("gitconfig-plain is not plain")
expectHome <- Sys.getenv("K245_EXPECT_HOME")
say("step0.HOME", Sys.getenv("HOME"))
if (!nzchar(expectHome) || !identical(Sys.getenv("HOME"), expectHome)) fault("HOME is not the empty folder")
homeGit <- file.exists(file.path(expectHome, c(".gitconfig", ".git-credentials", file.path(".config", "git"))))
if (any(homeGit)) fault("HOME holds git configuration or credentials")
if (isWindows) {
    say("step0.USERPROFILE", Sys.getenv("USERPROFILE"))
    if (!identical(Sys.getenv("USERPROFILE"), expectHome)) fault("USERPROFILE is not the empty folder")
}
workspace <- Sys.getenv("K245_WORKSPACE")
say("step0.getwd", getwd())
if (nzchar(workspace) && startsWith(normalizePath(getwd(), winslash = "/"), normalizePath(workspace, winslash = "/"))) {
    fault("R started inside the checkout")
}
if (hasGit) {
    q <- gitText(c("config", "--global", "--get-regexp", "^(url|http)\\."))
    say("step0.git_config_global_url_http", if (length(q$text)) q$text else "<none>")
    if (!q$status %in% c(0L, 1L)) fault(paste("git config --get-regexp failed with status", q$status))
    if (caseSet == "main") {
        target <- paste0("url.file://", remoteRoot, "/.insteadof https://x-access-token:", DUMMY, "@github.com/")
        if (!any(q$text == target)) fault("git does not see the redirect to the remote root")
        if (!any(q$text == "http.proxy http://127.0.0.1:9")) fault("git does not see the dead proxy")
    } else if (length(q$text)) {
        fault("git sees url or http configuration under gitconfig-plain")
    }
    s <- gitText(c("config", "--system", "--list"))
    say("step0.git_config_system", if (length(s$text)) s$text else "<none>")
    a <- gitText(c("config", "--list", "--show-origin"))
    say("step0.git_config_list_status", a$status)
    say("step0.git_config_all_with_origin", if (length(a$text)) a$text else "<none>")
    if (any(grepl("credential", a$text, ignore.case = TRUE))) fault("a credential setting is in force for git")
}
say("step0", "passed")

# ---- (1) The environment.
say("R", R.version.string)
say("os", paste(Sys.info()[["sysname"]], Sys.info()[["release"]]))
say("interactive()", interactive())
say("Sys.which(git)", if (hasGit) gitPath else "<none>")
if (hasGit) say("git --version", gitText("--version")$text)
say("PATH", Sys.getenv("PATH"))
say("MSR_BUNDLE", if (nzchar(Sys.getenv("MSR_BUNDLE"))) "set" else "unset")
say("tempdir()", tempdir())
for (p in c("git2r", "gh", "gitcreds", "MultiScholaR")) {
    say(paste0("packageVersion(", p, ")"), tryCatch(as.character(utils::packageVersion(p)), error = conditionMessage))
    say(paste0("find.package(", p, ")"), tryCatch(find.package(p), error = conditionMessage))
}
say("MultiScholaR SourceCommit",
    utils::packageDescription("MultiScholaR", fields = "Config/MultiScholaR/SourceCommit"))
say("git2r::libgit2_features()", tryCatch({
    f <- git2r::libgit2_features()
    paste(names(f), unlist(f), sep = "=")
}, error = function(e) paste("error:", conditionMessage(e))))
suppressPackageStartupMessages(library(MultiScholaR))

# ---- The recorder in gh's namespace (never calls the original).
REC <- new.env(parent = emptyenv())
REC$calls <- character()
recorder <- function(endpoint, ..., .token = NULL) {
    tok <- if (is.null(.token) || !nzchar(as.character(.token))) "none" else
        if (identical(as.character(unclass(.token)), DUMMY)) "dummy" else "other"
    a <- list(...)
    f <- function(x) if (is.null(x)) "NULL" else paste(format(x), collapse = ",")
    REC$calls <- c(REC$calls, sprintf("%s [.token=%s; org=%s; name=%s; private=%s; auto_init=%s]", endpoint, tok,
        f(a$org), f(a$name), f(a$private), f(a$auto_init)))
    if (tok != "dummy") {
        say("RECORDER_FOREIGN_TOKEN", paste("endpoint", endpoint, "token", tok))
        fault("the recorder saw a token that is not the dummy (stop rule 7)")
    }
    list(name = a$name, full_name = paste0(a$org, "/", a$name), private = TRUE)
}
ghNs <- asNamespace("gh")
unlockBinding("gh", ghNs)
assign("gh", recorder, envir = ghNs)
lockBinding("gh", ghNs)
if (!identical(gh::gh, recorder)) fault("gh::gh is not the recorder")
say("recorder", "gh::gh replaced in gh's namespace")

# ---- Fixtures (outside tempdir(), so the scan of tempdir() sees only what the function leaves).
template <- system.file("reports", "proteomics", "DIANN_report.rmd", package = "MultiScholaR")
template2 <- system.file("reports", "proteomics", "DIANN_limpa_report.rmd", package = "MultiScholaR")
if (!nzchar(template) || !nzchar(template2)) fault("shipped templates not found")
makeProject <- function(id, rproj = TRUE, twoTemplates = FALSE, dataFiles = FALSE) {
    root <- file.path(outDir, paste0("project-", id))
    unlink(file.path(outDir, paste0("project-", id)), recursive = TRUE, force = TRUE, expand = FALSE)
    scripts <- file.path(root, "scripts", "proteomics")
    dir.create(scripts, recursive = TRUE)
    file.copy(template, file.path(scripts, "DIANN_report.rmd"))
    if (twoTemplates) file.copy(template2, file.path(scripts, "DIANN_limpa_report.rmd"))
    writeLines(c("Run\tgroup", "s1\tA", "s2\tB"), file.path(scripts, "design_matrix.tab"))
    writeLines("kit245 fixture parameters", file.path(scripts, "study_parameters.txt"))
    if (rproj) writeLines("Version: 1.0", file.path(root, paste0(id, ".Rproj")))
    if (dataFiles) {
        dir.create(file.path(scripts, "import_sources"))
        writeLines(c("id\tvalue", "x\t1"), file.path(scripts, "import_sources", "x.tsv"))
        writeLines(c("id\tvalue", "x\t1"), file.path(scripts, "data_cln_LCMS_Pos.tab"))
        writeLines("placeholder archive", file.path(scripts, "filtered_session_x.msdia.tar"))
        saveRDS(list(kit245 = TRUE), file.path(scripts, "session_state_x.RDS"))
    }
    list(base_dir = root, source_dir = scripts)
}
filesWithDummy <- function(root) {
    f <- list.files(root, recursive = TRUE, all.files = TRUE, full.names = TRUE, no.. = TRUE)
    hit <- vapply(f, function(p) {
        size <- file.size(p)
        if (is.na(size) || size == 0) return(FALSE)
        length(grepRaw(DUMMY, readBin(p, "raw", size), fixed = TRUE)) > 0L
    }, logical(1))
    rel <- substring(f[hit], nchar(root) + 2L)
    if (length(rel)) rel else "<none>"
}
runPush <- function(case, fx, id) {
    REC$calls <- character()
    said <- character()
    err <- NULL
    status <- withCallingHandlers(
        tryCatch({
            r <- pushProjectToGithub(base_dir = fx$base_dir, source_dir = fx$source_dir, project_id = id,
                github_org = ORG, github_user_email = "kit245@example.invalid", github_user_name = "kit245-probe")
            paste("returned", format(r))
        }, error = function(e) {
            err <<- e
            "error"
        }),
        message = function(m) {
            said <<- c(said, conditionMessage(m))
            invokeRestart("muffleMessage")
        },
        warning = function(w) {
            said <<- c(said, paste("warning:", conditionMessage(w)))
            invokeRestart("muffleWarning")
        })
    say(paste0(case, ".status"), status)
    say(paste0(case, ".error_class"), if (is.null(err)) "<none>" else class(err))
    say(paste0(case, ".error_message"), if (is.null(err)) "<none>" else conditionMessage(err))
    say(paste0(case, ".messages"), if (length(said)) said else "<none>")
    say(paste0(case, ".recorder_calls"), length(REC$calls))
    say(paste0(case, ".recorder"), if (length(REC$calls)) REC$calls else "<none>")
    staged <- file.path(tempdir(), id)
    say(paste0(case, ".staged_folder_remains"), yn(dir.exists(staged)))
    cfg <- file.path(staged, ".git", "config")
    say(paste0(case, ".staged_git_config_holds_dummy"), if (file.exists(cfg)) {
        yn(length(grepRaw(DUMMY, readBin(cfg, "raw", file.size(cfg)), fixed = TRUE)) > 0L)
    } else "no file")
    say(paste0(case, ".files_under_tempdir_holding_dummy"), filesWithDummy(tempdir()))
    invisible(list(status = status, err = err, staged = staged))
}
rerunGitPush <- function(case, staged) {
    # The exact command of :331-338, with git's own standard error captured into R (base sends it to R's stderr).
    if (!hasGit || !dir.exists(staged)) {
        say(paste0(case, ".git_rerun"), if (!hasGit) "skipped: no git" else "skipped: no staged folder")
        return(invisible())
    }
    r <- gitText(c("-C", normalizePath(staged, winslash = "/"), "push", "origin", "refs/heads/main"))
    say(paste0(case, ".git_rerun_status"), r$status)
    say(paste0(case, ".git_rerun_text"), if (length(r$text)) r$text else "<none>")
    say(paste0(case, ".git_rerun_text_holds_dummy"), yn(any(grepl(DUMMY, r$text, fixed = TRUE))))
}
runCase <- function(case, expr) {
    tryCatch(expr, error = function(e) say(paste0(case, ".probe_error"), conditionMessage(e)))
}

if (caseSet == "main") {
    # ---- (2) T0: no token variable.
    runCase("T0", {
        t0 <- tryCatch(gh::gh_token(), error = function(e) e)
        if (inherits(t0, "error")) {
            say("T0.gh_token", paste("error:", paste(class(t0), collapse = "/"), "-", conditionMessage(t0)))
        } else {
            say("T0.gh_token_class", class(t0))
            if (nzchar(as.character(t0))) fault("gh::gh_token() found a token with no token variable set")
            say("T0.gh_token", "empty")
        }
        runPush("T0.push", makeProject("k245-t0", rproj = FALSE), "k245-t0")
    })

    # From here on the dummy is the token.
    Sys.setenv(GITHUB_PAT = DUMMY)
    tk <- tryCatch(gh::gh_token(), error = function(e) e)
    if (inherits(tk, "error") || !identical(as.character(unclass(tk)), DUMMY)) {
        say("token.gh_token_with_dummy", if (inherits(tk, "error")) paste("error:", conditionMessage(tk)) else "not the dummy")
        unlockBinding("gh_token", ghNs)
        assign("gh_token", function(...) DUMMY, envir = ghNs)
        lockBinding("gh_token", ghNs)
        say("token.route", "gh::gh_token stood in to return the dummy")
    } else {
        say("token.route", "GITHUB_PAT = dummy; gh::gh_token() returns it")
    }

    # ---- (3) P1: one template, no .Rproj.
    runCase("P1", runPush("P1", makeProject("k245-p1", rproj = FALSE), "k245-p1"))

    # ---- (4) P2: one template and .Rproj; the remote is redirected to a path where no repository exists.
    runCase("P2", {
        p2 <- runPush("P2", makeProject("k245-p2"), "k245-p2")
        rerunGitPush("P2", p2$staged)
    })

    # ---- (6) P2g: git present; the redirect leads to a local bare repository; the pushed tree.
    runCase("P2g", {
        if (!hasGit) {
            say("P2g", "skipped: no git")
        } else {
            bare <- file.path(remoteRoot, ORG, "k245-p2g.git")
            dir.create(bare, recursive = TRUE, showWarnings = FALSE)
            git2r::init(bare, bare = TRUE)
            runPush("P2g", makeProject("k245-p2g", dataFiles = TRUE), "k245-p2g")
            tree <- gitText(c("--git-dir", bare, "ls-tree", "-r", "--name-only", "refs/heads/main"))
            say("P2g.pushed_tree_git_ls_tree", if (length(tree$text)) tree$text else "<none>")
            say("P2g.pushed_tree_git2r_ls_tree", tryCatch({
                b <- git2r::repository(bare)
                br <- git2r::branches(b)[["main"]]
                lt <- git2r::ls_tree(git2r::tree(git2r::lookup(b, git2r::branch_target(br))), recursive = TRUE)
                paste0(lt$path, lt$name)[lt$type == "blob"]
            }, error = function(e) paste("error:", conditionMessage(e))))
            say("P2g.remote_files_holding_dummy", filesWithDummy(remoteRoot))
        }
    })

    # ---- (7) P3: two templates in the project.
    runCase("P3", runPush("P3", makeProject("k245-p3", twoTemplates = TRUE), "k245-p3"))

    # ---- (8) P4: (a)'s feasibility, base code untouched: git2r::push to a local bare repository with credentials.
    runCase("P4", {
        work <- file.path(outDir, "p4-work")
        bare4 <- file.path(outDir, "p4-remote.git")
        unlink(work, recursive = TRUE, force = TRUE, expand = FALSE)
        unlink(bare4, recursive = TRUE, force = TRUE, expand = FALSE)
        dir.create(work)
        dir.create(bare4)
        repo <- git2r::init(work)
        git2r::config(repo, user.name = "kit245-probe", user.email = "kit245@example.invalid")
        writeLines("kit245", file.path(work, "a.txt"))
        git2r::add(repo, "a.txt")
        git2r::commit(repo, "kit245 P4")
        git2r::init(bare4, bare = TRUE)
        git2r::remote_add(repo, "origin", normalizePath(bare4, winslash = "/"))
        branch <- git2r::repository_head(repo)$name
        r <- tryCatch({
            git2r::push(repo, "origin", paste0("refs/heads/", branch),
                credentials = git2r::cred_user_pass("x-access-token", DUMMY))
            "pushed"
        }, error = function(e) paste("error:", conditionMessage(e)))
        say("P4.git2r_push", r)
        say("P4.remote_branches", names(git2r::branches(git2r::repository(bare4))))
    })

    # ---- (10) P6 (W1', audit I01, I02): Windows with git present. The runner image's own Git for Windows is put on
    # R's PATH for this case only (K245_WINDOWS_GIT_DIR, set by the step); the guards are the redirect of
    # https://github.com/ to a local bare repository, GIT_ALLOW_PROTOCOL=file and the dead proxy, checked first.
    runCase("P6", (function() {
        gitDir <- Sys.getenv("K245_WINDOWS_GIT_DIR")
        if (!isWindows || !nzchar(gitDir)) {
            say("P6", "skipped: not the Windows invocation")
        } else {
            oldPath <- Sys.getenv("PATH")
            Sys.setenv(PATH = paste(gitDir, oldPath, sep = ";"), GCM_INTERACTIVE = "never")
            on.exit(Sys.setenv(PATH = oldPath), add = TRUE)
            say("P6.Sys.which(git)", Sys.which("git"))
            say("P6.git --version", gitText("--version")$text)
            q <- gitText(c("config", "--global", "--get-regexp", "^(url|http)\\."))
            say("P6.git_config_global_url_http", q$text)
            if (!any(grepl("insteadof https://github.com/$", q$text)) ||
                !any(q$text == "http.proxy http://127.0.0.1:9")) {
                fault("P6: git does not see the redirect of https://github.com/ and the dead proxy")
            }
            # The push, through the function under test, to the local bare repository.
            bare <- file.path(remoteRoot, ORG, "k245-p6.git")
            dir.create(bare, recursive = TRUE, showWarnings = FALSE)
            git2r::init(bare, bare = TRUE)
            runPush("P6", makeProject("k245-p6", dataFiles = TRUE), "k245-p6")
            tree <- gitText(c("--git-dir", bare, "ls-tree", "-r", "--name-only", "refs/heads/main"))
            say("P6.pushed_tree", if (length(tree$text)) tree$text else "<none>")
            # The credential mechanism with no network: the function's own git options (fix) or the audited head's
            # environment block (base), under a global configuration of two decoy helpers.
            home <- file.path(outDir, "p6-home")
            dir.create(home, showWarnings = FALSE)
            store <- file.path(home, "decoy-store")
            decoy <- file.path(home, "gitconfig-decoy")
            writeLines(c("[credential]", "\thelper = \"!f() { echo username=decoy; echo password=decoy; }; f\"",
                sprintf("\thelper = \"store --file %s\"", normalizePath(store, winslash = "/", mustWork = FALSE))), decoy)
            oldGlobal <- Sys.getenv("GIT_CONFIG_GLOBAL")
            Sys.setenv(GIT_CONFIG_GLOBAL = decoy)
            on.exit(Sys.setenv(GIT_CONFIG_GLOBAL = oldGlobal), add = TRUE)
            ns <- asNamespace("MultiScholaR")
            fixed <- exists("githubPushCredentialHelper", envir = ns, inherits = FALSE)
            say("P6.mechanism", if (fixed) "git -c options (fix)" else "GIT_CONFIG_COUNT environment block (audited head)")
            fill <- function(host, operation = "fill", extra = character()) {
                input <- c("protocol=https", paste0("host=", host), extra, "")
                if (fixed) {
                    cmd <- sprintf("git -c credential.helper= -c %s credential %s",
                        shQuote(get("githubPushCredentialHelper", envir = ns)()), operation)
                    out <- withr::with_envvar(c(MSR_GITHUB_PUSH_TOKEN = DUMMY),
                        suppressWarnings(system(cmd, intern = TRUE, input = input)))
                } else {
                    out <- withr::with_envvar(c(GIT_CONFIG_COUNT = "2", GIT_CONFIG_KEY_0 = "credential.helper",
                        GIT_CONFIG_VALUE_0 = "", GIT_CONFIG_KEY_1 = "credential.helper",
                        GIT_CONFIG_VALUE_1 = "!f() { echo username=x-access-token; echo \"password=$MSR_GITHUB_PUSH_TOKEN\"; }; f",
                        MSR_GITHUB_PUSH_TOKEN = DUMMY),
                        suppressWarnings(system(paste("git credential", operation), intern = TRUE, input = input)))
                }
                list(status = if (is.null(attr(out, "status"))) 0L else attr(out, "status"), text = out)
            }
            gh <- fill("github.com")
            say("P6.credential_github_status", gh$status)
            say("P6.credential_github_gives_dummy", yn(paste0("password=", DUMMY) %in% gh$text))
            say("P6.credential_github_gives_decoy", yn("password=decoy" %in% gh$text))
            other <- fill("example.com")
            say("P6.credential_example_gives_dummy", yn(any(grepl(DUMMY, other$text, fixed = TRUE))))
            fill("github.com", "approve", c("username=x-access-token", paste0("password=", DUMMY)))
            say("P6.decoy_store_written", yn(file.exists(store)))
        }
    })())

    # ---- (9) P5: (b)'s feasibility: the Summary UI of the four omics.
    runCase("P5", {
        for (omic in c("proteomics", "lipidomics", "metabolomics", "transcriptomics")) {
            html <- as.character(MultiScholaR:::mod_omic_summary_ui("summary",
                MultiScholaR:::omicGuiConfig(omic)$summary))
            say(paste0("P5.", omic), sprintf("push_to_github=%s; enable_github=%s; github_org=%s; heading=%s; not_available_sentence=%s",
                yn(grepl("id=\"summary-push_to_github\"", html, fixed = TRUE)),
                yn(grepl("id=\"summary-enable_github\"", html, fixed = TRUE)),
                yn(grepl("id=\"summary-github_org\"", html, fixed = TRUE)),
                yn(grepl("Version Control (Optional)", html, fixed = TRUE)),
                yn(grepl("Pushing to GitHub is not available in the installed application.", html, fixed = TRUE))))
        }
    })
} else {
    # ---- (5) P2n: no network, no redirect, GIT_ALLOW_PROTOCOL unset: git attempts HTTPS to github.com.
    Sys.setenv(GITHUB_PAT = DUMMY)
    tk <- tryCatch(gh::gh_token(), error = function(e) e)
    if (inherits(tk, "error") || !identical(as.character(unclass(tk)), DUMMY)) {
        unlockBinding("gh_token", ghNs)
        assign("gh_token", function(...) DUMMY, envir = ghNs)
        lockBinding("gh_token", ghNs)
        say("token.route", "gh::gh_token stood in to return the dummy")
    } else {
        say("token.route", "GITHUB_PAT = dummy; gh::gh_token() returns it")
    }
    runCase("P2n", {
        p2n <- runPush("P2n", makeProject("k245-p2n"), "k245-p2n")
        rerunGitPush("P2n", p2n$staged)
    })
}

say("after.files_under_tempdir_holding_dummy", filesWithDummy(tempdir()))
cat("G0 END", host, "\n")
