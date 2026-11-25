# MultiScholaR test build: a guide for testers

This build is for testing MultiScholaR on public data before release. It is the normal application with two
pre-release features that you can switch on: the **RUV design panel** (new choices for removing unwanted variation, on
each omic's Normalisation tab) and the **Integration module** (combining several omics with MOFA2). Both are off unless
you switch on the test setting below.

## 1. Install

**Linux** (Ubuntu 22.04 or newer, Debian 12 or newer, Fedora 36 or newer; no administrator rights needed)

1. Download `MultiScholaR-<version>-linux-x86_64.tar.gz` and unpack it: `tar -xzf MultiScholaR-<version>-linux-x86_64.tar.gz`
2. In the unpacked folder run `./install.sh`. It installs into `~/.local/share/multischolar/<version>` and adds
   MultiScholaR to your applications menu.
3. If it says `/usr/bin/which` is missing, install it (Fedora: `sudo dnf install which`; Arch: `sudo pacman -S which`)
   and run `./install.sh` again.
4. To remove it: `~/.local/share/multischolar/<version>/uninstall.sh`.

**Windows** (Windows 10 or 11, 64-bit; no administrator rights needed)

1. Download `MultiScholaR-Setup-<version>.exe` and double-click it. It installs for your user only and adds a
   Start-menu entry.
2. The installer is not signed yet, so Windows may show "Windows protected your PC". Click **More info**, then
   **Run anyway**.
3. To remove it: Settings, Apps, MultiScholaR, Uninstall.

Nothing is downloaded when the application starts; everything it needs is inside the installation.

## 2. Switch on the test features

The test setting is the environment variable `MULTISCHOLAR_TEST_FEATURES=1`. When it is on, the window title ends with
**(test features on)**, the Integration tile on the home page can be clicked, and each omic's Normalisation tab shows
the RUV design choices. Any other value, or no value, means off.

- **Linux, once:** start from a terminal: `MULTISCHOLAR_TEST_FEATURES=1 ~/.local/share/multischolar/<version>/multischolar`
- **Linux, always:** add the line `export MULTISCHOLAR_TEST_FEATURES=1` to `~/.profile`, log out and in again, then
  start from the menu. To switch off, remove the line.
- **Windows:** open PowerShell and run `setx MULTISCHOLAR_TEST_FEATURES 1`, then start MultiScholaR from the Start menu
  (close it first if it is running). To switch off: `setx MULTISCHOLAR_TEST_FEATURES 0`.

On Linux computers that can run each step under a memory limit (systemd, kernel 6.12 or newer), the test setting runs
steps without that limit. This is expected for the test build.

## 3. Where your projects and results go

- **Your project:** on the first page you choose a *Project Output Directory* (your home folder unless you change it)
  and an *experiment label*. MultiScholaR creates a folder `<directory>/<label>` and writes everything for that project
  inside it: the results tables, figures, reports and the saved sessions.
- **On Windows keep this path short**, for example `C:\MSR\study1`. With Windows' default settings, project paths longer
  than about 75 characters can make a step stop with a message that a path is too long.
- **Logs:** Help, *Open the logs folder*. One `shell-<time>.log` and one `r-<time>.log` per start
  (Linux `~/.config/MultiScholaR/logs`, Windows `%APPDATA%\MultiScholaR\logs`).
- **Working files of a session:** `~/.local/state/multischolar/run/<time>` (Linux) or
  `%LOCALAPPDATA%\MultiScholaR\run\<time>` (Windows); the last five are kept.

## 4. The STRING key, KEGG and Reactome

- STRING needs a personal key (free from string-db.org). In the application, type it into **Integration, Enrichment,
  Your STRING key**. It is used for that one query, the field is emptied straight after, and it is never saved in the
  project, the record or the logs. Type it again for the next query.
- If you run the integration from R with the package's `integration*` functions instead of the application, set the
  environment variable `STRING_API_KEY` before you start R and pass it as
  `integrationEnrichFactor(..., api_key = Sys.getenv("STRING_API_KEY"))`; never write the key in a script. It is not
  saved either.
- STRING, KEGG and Reactome are online services: enrichment needs an internet connection. The rest of the application
  works offline.

## 5. When something goes wrong: what to send

Send: what you did (omic, tab, button, the settings you chose), what you expected and what happened, the exact message
on screen (a screenshot is best), the version (Help, *About MultiScholaR*), your operating system, and the two log files
of that start (Help, *Open the logs folder*; the newest `shell-` and `r-` files). If a step stopped, also the
`*.log` files from that session's working-files folder (section 3).

Do **not** send: your data files, result tables, reports, saved sessions (`.tar`, `.rds`, `.qs2`), or the project
folder. Never send your STRING key. Logs can contain file and sample names: look through them first and replace
anything you cannot share.

## 6. Known limits of this build

- The Integration module and the RUV design panel are **pre-release**: their wording, defaults and figures may still
  change. The performance figures shown beside the RUV choices come from simulations and will be updated with public
  data before release.
- **DIA proteomics** (DIA, DIA with limpa, and enrichment after a DIA analysis) runs on **Linux only**; on Windows these
  routes stop with a message.
- **Linux, DIA with limpa:** the fast engine is built for the processor of the build machine and can stop with
  "illegal instruction" on older processors. If that happens, choose the `limpa` engine instead and tell us.
- **Windows:** long project paths (section 3); the installer is unsigned (section 1).
- **Linux:** needs glibc 2.35 or newer and `/usr/bin/which`; RHEL, Rocky and Alma 9 and older, Ubuntu 20.04 and Debian
  11 are not supported.
- **Lipidomics and metabolomics:** cyclic loess normalisation is unstable on data with only two or three features;
  with four or more it behaves normally.
- No automatic updates: install each new test build as it comes.
