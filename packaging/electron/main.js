'use strict';
// KIT-176: ELECTRON_RUN_AS_NODE from the caller's environment (a shell an editor extension or a coding agent
// opened). The binary ignores it (the runAsNode fuse is off, package.json electronFuses); it is kept here for
// the start event and removed, so nothing the shell starts (R, the probe, the system's opener) inherits it.
const RUN_AS_NODE_VARIABLE = process.env.ELECTRON_RUN_AS_NODE;
delete process.env.ELECTRON_RUN_AS_NODE;
// MultiScholaR desktop shell (KIT-134 P3, revision 2.0, 2.0.5). From the KIT-130 spike's main.js.
//
// One window that shows only the application served on 127.0.0.1 by the bundled R. The shell:
// takes a single-instance lock (a second start focuses the first window and starts no R); probes
// the host once for the launcher the package can use; starts the bundled Rscript (--vanilla, with
// R's home, library and Python inside the bundle and nothing from the user's R, packages, profile
// or Python) on a free port; waits for HTTP 200; refuses navigation and new windows to any other
// origin (https links open in the system browser); and, when the window closes, stops the whole R
// process tree. Logs go to Electron's logs folder, one shell log and one R log per start; the run
// folder (the application's working files) to %LOCALAPPDATA%\MultiScholaR\run\<stamp> or
// $XDG_STATE_HOME/multischolar/run/<stamp>, the last five kept.

const { app, BrowserWindow, Menu, dialog, shell } = require('electron');
const { spawn, spawnSync } = require('child_process');
const fs = require('fs');
const http = require('http');
const net = require('net');
const os = require('os');
const path = require('path');

const isWindows = process.platform === 'win32';
const resources = app.isPackaged ? process.resourcesPath : __dirname;
const runtime = path.join(resources, 'runtime');
const entryScript = path.join(resources, 'entry.R');
const stamp = readJson(path.join(resources, 'VERSION.json')) || { release: app.getVersion() };
const startedAt = new Date().toISOString().replace(/[:.]/g, '-');
const KEEP_RUN_FOLDERS = 5;

// Lifecycle-test hooks, set only by packaging/electron/lifecycle-test.js: an event log to append
// to, and a file whose appearance closes the window (the same path as the user closing it).
const testLog = process.env.MSR_LIFECYCLE_LOG;
const testCloseFile = process.env.MSR_LIFECYCLE_CLOSE;
if (testLog) app.disableHardwareAcceleration(); // the runners have no GPU

let shellLog = null;
let rProcess = null;
let rExited = false;
// The Linux launcher (packaging/linux/Rscript) exits 69 when the host has no /usr/bin/which (decision M25).
let rExitCode = null;
const MISSING_WHICH = 'MultiScholaR needs the program "which" (/usr/bin/which), which this computer does not have. '
  + 'Install it (Fedora, RHEL: sudo dnf install which; Arch Linux: sudo pacman -S which; Debian, Ubuntu: '
  + 'sudo apt install debianutils; openSUSE: sudo zypper install which), then start MultiScholaR again.';
let mainWindow = null;
let origin = null;
let launcher = { launcher: 'unknown', reason: 'R has not reported its launcher' };
// KIT-156: set when R reports the test setting (MULTISCHOLAR_TEST_FEATURES=1, packaging/electron/entry.R).
let testFeatures = null;
let quitting = false;

function readJson(file) {
  try { return JSON.parse(fs.readFileSync(file, 'utf8')); } catch (error) { return null; }
}

function record(event, data = {}) {
  const line = JSON.stringify({ time: new Date().toISOString(), event, ...data });
  if (shellLog) shellLog.write(`${line}\n`);
  if (testLog) fs.appendFileSync(testLog, `${line}\n`);
  if (!app.isPackaged) console.log(line);
}

// ------------------------------------------------------------------- folders
function stateRoot() {
  if (isWindows) return path.join(process.env.LOCALAPPDATA || app.getPath('appData'), 'MultiScholaR');
  return path.join(process.env.XDG_STATE_HOME || path.join(os.homedir(), '.local', 'state'), 'multischolar');
}

function newRunFolder() {
  const parent = path.join(stateRoot(), 'run');
  fs.mkdirSync(parent, { recursive: true });
  const folder = path.join(parent, startedAt);
  fs.mkdirSync(folder, { recursive: true });
  const old = fs.readdirSync(parent).sort();
  for (const name of old.slice(0, Math.max(0, old.length - KEEP_RUN_FOLDERS))) {
    try { fs.rmSync(path.join(parent, name), { recursive: true, force: true }); } catch (error) { /* in use */ }
  }
  return folder;
}

// ----------------------------------------------------------- the start-up probe
// The host facts that decide whether the package's systemd worker can run here (revision 2.0,
// 2.0.5, step 1): Linux; systemctl, systemd-run and busctl; a user manager that answers; the
// unified cgroup v2 hierarchy; the memory controller delegated to the user manager; kernel 6.12 or
// later (R/utils_workflow_execution.R, openWorkflowPeakWindow); and a user D-Bus on which systemd
// answers. Otherwise the plain launcher. The package's entry takes the result (MSR_LAUNCHER_PROBE)
// and owns the rest of the order: the trial unit, and the in-process fallback.
function probeLauncher() {
  const started = Date.now();
  const facts = { platform: process.platform, kernel: os.release() };
  const verdict = (choice, reason) => ({ launcher: choice, reason, facts, ms: Date.now() - started });
  if (process.platform !== 'linux') return verdict('process', `${process.platform}: no systemd`);
  const onPath = (name) => ['/usr/bin', '/bin'].some((dir) => fs.existsSync(path.join(dir, name)));
  facts.tools = Object.fromEntries(['systemctl', 'systemd-run', 'busctl'].map((name) => [name, onPath(name)]));
  if (!Object.values(facts.tools).every(Boolean)) return verdict('process', 'systemctl, systemd-run or busctl is missing');
  const manager = spawnSync('/usr/bin/systemctl', ['--user', 'show-environment'], { timeout: 1500, stdio: 'ignore' });
  facts.userManager = manager.status === 0;
  if (!facts.userManager) return verdict('process', 'the user service manager does not answer');
  let cgroup = '';
  try { cgroup = fs.readFileSync('/proc/self/cgroup', 'utf8').trim(); } catch (error) { /* none */ }
  facts.cgroup = cgroup;
  if (!/^0::[^\n]*$/.test(cgroup)) return verdict('process', 'not the unified cgroup v2 hierarchy');
  const uid = process.getuid();
  const managerCgroup = `/sys/fs/cgroup/user.slice/user-${uid}.slice/user@${uid}.service`;
  let controllers = '';
  try { controllers = fs.readFileSync(path.join(managerCgroup, 'cgroup.controllers'), 'utf8'); } catch (error) { /* none */ }
  facts.controllers = controllers.trim();
  if (!controllers.split(/\s+/).includes('memory')) return verdict('process', 'the memory controller is not delegated to the user manager');
  const [major, minor] = os.release().split(/[.-]/).map((part) => Number.parseInt(part, 10));
  if (!(major > 6 || (major === 6 && minor >= 12))) return verdict('process', `kernel ${major}.${minor} is older than 6.12`);
  // Last, the user D-Bus itself, not only busctl: the package creates the application's slice with
  // busctl --user (createWorkflowApplicationSlice). systemctl reaches the user manager over its
  // private socket, so it answers on a host without a user bus (GitHub's ubuntu-24.04 runner, runs
  // 44055833237 and 44059705566), where the slice call then fails with "Failed to connect to bus".
  const bus = spawnSync('/usr/bin/busctl', ['--user', 'call', 'org.freedesktop.systemd1', '/org/freedesktop/systemd1',
    'org.freedesktop.DBus.Peer', 'Ping'], { timeout: 1500, stdio: 'ignore' });
  facts.userBus = bus.status === 0;
  if (!facts.userBus) return verdict('process', 'the user D-Bus does not answer');
  return verdict('systemd', 'systemd worker: every host requirement holds');
}

// ------------------------------------------------------------------------- R
// The user's environment minus everything that could steer R or Python to another installation.
// On Linux the launcher inside the bundle (packaging/linux/Rscript) sets R's environment; on
// Windows it is set here (KIT-130, "What the launcher must set").
function rEnvironment(probe, runFolder) {
  const steering = /^(PATH|RHOME|R_.*|RTOOLS.*|RETICULATE_.*|PYTHON.*|CONDA_.*|VIRTUAL_ENV|BASILISK_.*|LD_LIBRARY_PATH|MSR_LIFECYCLE_.*)$/i;
  const env = {};
  for (const [key, value] of Object.entries(process.env)) {
    if (!steering.test(key)) env[key] = value;
  }
  env.MSR_LAUNCHER_PROBE = JSON.stringify(probe);
  env.MSR_RUN_FOLDER = runFolder;
  if (isWindows) {
    const systemRoot = process.env.SystemRoot || 'C:\\Windows';
    // Windows PowerShell's folder too: shinyFiles lists the drives with powershell where WMIC is
    // gone (Windows Server 2025, Windows 11 24H2), and the session fails without it (run 43945108934).
    env.PATH = `${systemRoot}\\System32;${systemRoot};${systemRoot}\\System32\\WindowsPowerShell\\v1.0`;
    // Neither folder exists: R keeps only its own library on the search path.
    env.R_LIBS_USER = path.join(runtime, 'R', 'no-user-library');
    env.R_LIBS_SITE = path.join(runtime, 'R', 'no-site-library');
    // R for Windows puts Rtools' folders at the front of its own PATH; this one does not exist.
    env.RTOOLS45_HOME = path.join(runtime, 'no-rtools');
    env.RETICULATE_PYTHON = path.join(runtime, 'python', 'python.exe');
    env.PYTHONNOUSERSITE = '1';
    env.RSTUDIO_PANDOC = path.join(runtime, 'pandoc');
    // The saved-entity store must not take WSL's bash or other GNU tools found on the host (F6).
    env.MULTISCHOLAR_DISABLE_TOOLS = 'bash,dd,sha256sum,sync,tar,cp';
    env.MSR_BUNDLE = runtime;
  }
  return env;
}

function freePort() {
  return new Promise((resolve, reject) => {
    const server = net.createServer();
    server.once('error', reject);
    server.listen(0, '127.0.0.1', () => {
      const { port } = server.address();
      server.close(() => resolve(port));
    });
  });
}

function startR(port, runFolder, probe) {
  const rscript = isWindows ? path.join(runtime, 'R', 'bin', 'Rscript.exe') : path.join(runtime, 'bin', 'Rscript');
  const rLog = fs.createWriteStream(path.join(app.getPath('logs'), `r-${startedAt}.log`), { flags: 'a' });
  rProcess = spawn(rscript, ['--vanilla', entryScript, String(port), runFolder, String(process.pid)], {
    env: rEnvironment(probe, runFolder),
    cwd: runFolder,
    stdio: ['ignore', 'pipe', 'pipe'],
    windowsHide: true,
    detached: !isWindows, // its own process group, so the whole group can be signalled
  });
  rProcess.once('error', (error) => {
    rExited = true;
    record('r-spawn-error', { message: String(error) });
  });
  rProcess.once('exit', (code, signal) => {
    rExited = true;
    rExitCode = code;
    record('r-exited', { code, signal });
    if (!quitting && mainWindow) {
      dialog.showErrorBox('MultiScholaR stopped', `The analysis engine stopped (exit ${code ?? signal}). Details are in the log:\n${app.getPath('logs')}`);
      app.quit();
    }
  });
  let pending = '';
  rProcess.stdout.on('data', (chunk) => {
    rLog.write(chunk);
    pending += chunk.toString();
    const lines = pending.split(/\r?\n/);
    pending = lines.pop();
    for (const line of lines) {
      if (line.startsWith('MSR-LAUNCHER ')) {
        try { launcher = JSON.parse(line.slice('MSR-LAUNCHER '.length)); } catch (error) { launcher = { launcher: 'unreadable', reason: line }; }
        record('launcher', launcher);
        applyLauncher();
      } else if (line.startsWith('MSR-TEST-FEATURES ')) {
        try { testFeatures = JSON.parse(line.slice('MSR-TEST-FEATURES '.length)); } catch (error) { testFeatures = { unreadable: line }; }
        record('test-features', testFeatures);
        if (mainWindow) mainWindow.setTitle(title());
      }
    }
  });
  rProcess.stderr.on('data', (chunk) => rLog.write(chunk));
  record('r-started', { pid: rProcess.pid, rscript, port, runFolder });
}

// Stops the R process tree: Linux, the process group, then SIGKILL after 5 s; Windows,
// taskkill /T /F (Rscript.exe runs Rterm.exe as a child; ending the parent alone leaves R running).
// Under the systemd worker the application runs in units of the user manager, outside R's process
// group: their slice (named on the MSR-LAUNCHER line) is stopped first.
function stopR() {
  if (launcher.launcher === 'systemd' && /^multischolarwfx[a-f0-9]{16}\.slice$/.test(launcher.slice || '')) {
    const stopped = spawnSync('/usr/bin/systemctl', ['--user', 'stop', launcher.slice], { timeout: 15000, stdio: 'ignore' });
    record('slice-stopped', { slice: launcher.slice, status: stopped.status });
  }
  return new Promise((resolve) => {
    if (!rProcess || rExited) {
      resolve();
      return;
    }
    const hardStop = setTimeout(() => {
      if (!isWindows) {
        try { process.kill(-rProcess.pid, 'SIGKILL'); } catch (error) { /* already gone */ }
      }
      resolve();
    }, 5000);
    rProcess.once('exit', () => {
      clearTimeout(hardStop);
      resolve();
    });
    if (isWindows) {
      const taskkill = path.join(process.env.SystemRoot || 'C:\\Windows', 'System32', 'taskkill.exe');
      spawnSync(taskkill, ['/pid', String(rProcess.pid), '/T', '/F'], { windowsHide: true });
    } else {
      try { process.kill(-rProcess.pid, 'SIGTERM'); } catch (error) { /* already gone */ }
    }
  });
}

function waitForApplication(port, timeoutMs = 300000) {
  const started = Date.now();
  return new Promise((resolve, reject) => {
    const attempt = () => {
      if (rExited) {
        reject(new Error(!isWindows && rExitCode === 69 ? MISSING_WHICH : 'R exited before the application answered'));
        return;
      }
      const retry = () => {
        if (Date.now() - started > timeoutMs) reject(new Error('the application did not answer in time'));
        else setTimeout(attempt, 500);
      };
      http.get({ host: '127.0.0.1', port, path: '/', timeout: 2000 }, (response) => {
        response.resume();
        if (response.statusCode === 200) resolve(Date.now() - started);
        else retry();
      }).on('error', retry).on('timeout', function onTimeout() { this.destroy(); });
    };
    attempt();
  });
}

// ------------------------------------------------------------------- window
function title() {
  const release = `MultiScholaR ${stamp.release}`;
  const base = testFeatures ? `${release} (test features on)` : release;
  return launcher.launcher === 'in-process' ? `${base} (steps run inside the application)` : base;
}

const IN_PROCESS_NOTICE = 'Steps run inside the application on this computer: memory used by a step is kept until the application is restarted. Details are in the log.';

// The in-process notice: once per start, in the page, dismissible, never blocking.
function applyLauncher() {
  if (!mainWindow) return;
  mainWindow.setTitle(title());
  if (launcher.launcher !== 'in-process' || !origin || !mainWindow.webContents.getURL().startsWith(origin)) return;
  const script = `(() => {
    if (document.getElementById('msr-launcher-notice')) return true;
    const bar = document.createElement('div');
    bar.id = 'msr-launcher-notice';
    bar.setAttribute('role', 'status');
    bar.style.cssText = 'position:fixed;bottom:12px;right:12px;max-width:420px;z-index:100000;background:#fff8e1;color:#4e3b00;border:1px solid #e0c060;border-radius:6px;padding:10px 32px 10px 12px;font:13px/1.4 sans-serif;box-shadow:0 2px 8px rgba(0,0,0,.2)';
    bar.textContent = ${JSON.stringify(IN_PROCESS_NOTICE)};
    const close = document.createElement('button');
    close.textContent = '\\u00d7';
    close.setAttribute('aria-label', 'Dismiss');
    close.style.cssText = 'position:absolute;top:4px;right:6px;border:0;background:none;font-size:18px;cursor:pointer;color:inherit';
    close.onclick = () => bar.remove();
    bar.appendChild(close);
    document.body.appendChild(bar);
    return true;
  })()`;
  mainWindow.webContents.executeJavaScript(script).then(() => record('notice-shown')).catch(() => {});
}

// The release, commit, lockfile, snapshot and build lines of About (KIT-176, D5): the commit lines are chosen by
// package_on_release_line first; a stamp without it (an unpackaged start) keeps the build commit's line.
function stampLines() {
  let commitLines;
  if (stamp.package_on_release_line === false) {
    commitLines = [`Analysis code: commit ${stamp.commit} (test build, not on the release line)`];
  } else if (stamp.package_on_release_line === true && stamp.package_commit === stamp.commit) {
    commitLines = [`Commit: ${stamp.package_commit}`];
  } else if (stamp.package_on_release_line === true) {
    commitLines = [`Analysis code: commit ${stamp.package_commit}`, `Test build from: commit ${stamp.commit}`];
  } else {
    commitLines = [`Commit: ${stamp.commit || 'unknown'}`];
  }
  return [
    `MultiScholaR ${stamp.release}`,
    ...commitLines,
    `Lockfile sha256: ${stamp.lockfile_sha256 || 'unknown'}`,
    `Package snapshot: ${stamp.snapshot || 'unknown'}`,
    `Built: ${stamp.build_date || 'unknown'} on ${stamp.runner_image || 'unknown'}`,
  ];
}

function showAbout() {
  const lines = [
    ...stampLines(),
    '',
    `How steps run: ${launcher.launcher}`,
    `Why: ${launcher.reason}`,
    '',
    `Logs: ${app.getPath('logs')}`,
    `Licences: ${path.join(resources, 'THIRD-PARTY-NOTICES.txt')}`,
  ];
  dialog.showMessageBox(mainWindow, { type: 'info', title: 'About MultiScholaR', message: `MultiScholaR ${stamp.release}`, detail: lines.join('\n') });
}

function buildMenu() {
  const template = [
    { label: 'File', submenu: [{ role: 'quit' }] },
    { label: 'Edit', submenu: [{ role: 'undo' }, { role: 'redo' }, { type: 'separator' }, { role: 'cut' }, { role: 'copy' }, { role: 'paste' }, { role: 'selectAll' }] },
    { label: 'View', submenu: [{ role: 'resetZoom' }, { role: 'zoomIn' }, { role: 'zoomOut' }, { type: 'separator' }, { role: 'togglefullscreen' }] },
    {
      label: 'Help',
      submenu: [
        { label: 'About MultiScholaR', click: showAbout },
        { label: 'Open the logs folder', click: () => shell.openPath(app.getPath('logs')) },
        { label: 'Third-party notices', click: () => shell.openPath(path.join(resources, 'THIRD-PARTY-NOTICES.txt')) },
      ],
    },
  ];
  Menu.setApplicationMenu(Menu.buildFromTemplate(template));
}

// Resolves with what the page shows once a live Shiny session has rendered the first tab: the
// session connected, the home tab's tiles present, and output rendered by the server.
const FIRST_TAB = `new Promise((resolve) => {
  const started = Date.now();
  const tick = () => {
    const connected = Boolean(window.Shiny && Shiny.shinyapp && Shiny.shinyapp.isConnected());
    const home = Boolean(document.getElementById('proteomics_box'));
    const output = document.getElementById('log_terminal_ui');
    const rendered = Boolean(output && output.children.length > 0);
    if ((connected && home && rendered) || Date.now() - started > 120000) {
      resolve({ connected, home, rendered, title: document.title });
    } else setTimeout(tick, 250);
  };
  tick();
})`;

function createWindow() {
  mainWindow = new BrowserWindow({
    width: 1400,
    height: 900,
    title: title(),
    show: true,
    webPreferences: { contextIsolation: true, nodeIntegration: false, sandbox: true },
  });
  mainWindow.on('page-title-updated', (event) => event.preventDefault());
  const contents = mainWindow.webContents;
  const external = (url) => {
    if (/^https:\/\//i.test(url)) shell.openExternal(url);
  };
  contents.on('will-navigate', (event, url) => {
    if (origin && url.startsWith(`${origin}/`)) return;
    event.preventDefault();
    external(url);
  });
  contents.setWindowOpenHandler(({ url }) => {
    if (origin && url.startsWith(`${origin}/`)) contents.downloadURL(url); // Shiny's download links
    else external(url);
    return { action: 'deny' };
  });
  contents.on('did-fail-load', (event, code, description, url) => record('did-fail-load', { code, description, url }));
  contents.on('did-finish-load', async () => {
    const url = contents.getURL();
    record('did-finish-load', { url });
    if (!origin || !url.startsWith(origin)) return;
    applyLauncher();
    const firstTab = await contents.executeJavaScript(FIRST_TAB).catch((error) => ({ error: String(error) }));
    record('first-tab', firstTab);
    if (testCloseFile) {
      const watcher = setInterval(() => {
        if (fs.existsSync(testCloseFile)) {
          clearInterval(watcher);
          if (mainWindow) mainWindow.close();
        }
      }, 300);
    }
  });
  mainWindow.on('closed', () => { mainWindow = null; });
  const waiting = `<!doctype html><meta charset="utf-8"><title>MultiScholaR</title><body style="font:15px sans-serif;display:flex;align-items:center;justify-content:center;height:90vh;color:#444">Starting MultiScholaR ${stamp.release}&hellip;</body>`;
  mainWindow.loadURL(`data:text/html;charset=utf-8,${encodeURIComponent(waiting)}`);
}

// ------------------------------------------------------------------ lifecycle
if (!app.requestSingleInstanceLock()) {
  // A second start joins the first: the first window comes to the front; no second R.
  if (testLog) fs.appendFileSync(testLog, `${JSON.stringify({ time: new Date().toISOString(), event: 'second-start-joined', pid: process.pid })}\n`);
  app.quit();
} else {
  app.on('second-instance', () => {
    record('second-instance');
    if (mainWindow) {
      if (mainWindow.isMinimized()) mainWindow.restore();
      mainWindow.focus();
    }
  });

  app.on('window-all-closed', async () => {
    quitting = true;
    record('window-closed');
    await stopR();
    record('quit', { rExited });
    app.quit();
  });

  app.on('before-quit', () => { quitting = true; });

  app.whenReady().then(async () => {
    fs.mkdirSync(app.getPath('logs'), { recursive: true });
    shellLog = fs.createWriteStream(path.join(app.getPath('logs'), `shell-${startedAt}.log`), { flags: 'a' });
    record('start', {
      release: stamp.release, resources, runtime, packaged: app.isPackaged, platform: process.platform,
      electron: process.versions.electron, pid: process.pid, noSandbox: app.commandLine.hasSwitch('no-sandbox'),
      logs: app.getPath('logs'), startedAt,
      runAsNodeVariable: RUN_AS_NODE_VARIABLE === undefined ? null : RUN_AS_NODE_VARIABLE, about: stampLines(),
    });
    buildMenu();
    try {
      const probe = probeLauncher();
      record('probe', probe);
      const runFolder = newRunFolder();
      const port = await freePort();
      origin = `http://127.0.0.1:${port}`;
      createWindow();
      startR(port, runFolder, probe);
      const waited = await waitForApplication(port);
      record('application-answered', { port, ms: waited });
      if (mainWindow) mainWindow.loadURL(`${origin}/`);
    } catch (error) {
      record('error', { message: String(error) });
      quitting = true;
      await stopR();
      if (!testLog) dialog.showErrorBox('MultiScholaR could not start', `${error.message}\n\nDetails are in the log:\n${app.getPath('logs')}`);
      app.exit(1);
    }
  });
}
