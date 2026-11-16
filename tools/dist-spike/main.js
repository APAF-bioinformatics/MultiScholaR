'use strict';
// KIT-130 packaging spike: an Electron shell that starts the hello-world Shiny app with the
// bundled Rscript on a free port, shows it in a window, and stops R when the window closes.

const { app, BrowserWindow } = require('electron');
const { spawn, spawnSync } = require('child_process');
const fs = require('fs');
const http = require('http');
const net = require('net');
const path = require('path');

const isWindows = process.platform === 'win32';
const resources = app.isPackaged ? process.resourcesPath : __dirname;
const runtime = path.join(resources, 'runtime');
const appScript = path.join(resources, 'app', 'run.R');

// Smoke-test hooks, set only by the workflow: an event log to write and a file whose
// appearance closes the window (the same path as the user closing it).
const smokeLog = process.env.SPIKE_LOG;
const closeFile = process.env.SPIKE_CLOSE_FILE;
if (smokeLog) {
  app.disableHardwareAcceleration(); // the runners have no GPU
}

let rProcess = null;
let rExited = false;

function record(event, data = {}) {
  const line = JSON.stringify({ event, ...data });
  console.log(line);
  if (smokeLog) fs.appendFileSync(smokeLog, `${line}\n`);
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

// The user's environment minus everything that could steer R or Python to another installation.
function rEnvironment() {
  const steering = /^(PATH|RHOME|R_.*|RETICULATE_.*|PYTHON.*|CONDA_.*|VIRTUAL_ENV|BASILISK_.*|LD_LIBRARY_PATH|SPIKE_.*)$/i;
  const env = {};
  for (const [key, value] of Object.entries(process.env)) {
    if (!steering.test(key)) env[key] = value;
  }
  if (isWindows) {
    // R for Windows finds its own home from the executable; only the library and Python need saying.
    env.PATH = `${process.env.SystemRoot}\\System32;${process.env.SystemRoot}`;
    env.R_LIBS_USER = path.join(runtime, 'R', 'no-user-library');
    env.R_LIBS_SITE = path.join(runtime, 'R', 'no-site-library');
    env.RETICULATE_PYTHON = path.join(runtime, 'python', 'python.exe');
    env.PYTHONNOUSERSITE = '1';
    env.MSR_BUNDLE = runtime;
  } else {
    // On Linux the launcher inside the bundle sets the same variables (tools/dist-spike/linux/Rscript).
    env.PATH = '/usr/bin:/bin';
  }
  return env;
}

function startR(port) {
  const rscript = isWindows
    ? path.join(runtime, 'R', 'bin', 'Rscript.exe')
    : path.join(runtime, 'bin', 'Rscript');
  rProcess = spawn(rscript, ['--vanilla', appScript, String(port)], {
    env: rEnvironment(),
    cwd: app.getPath('temp'),
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
    record('r-exited', { code, signal });
  });
  for (const stream of [rProcess.stdout, rProcess.stderr]) {
    stream.on('data', (chunk) => record('r-output', { text: chunk.toString().trimEnd() }));
  }
  record('r-started', { pid: rProcess.pid, rscript, port });
}

function stopR() {
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
      // Rscript.exe runs Rterm.exe as a child process: end the tree, not only the parent.
      const taskkill = path.join(process.env.SystemRoot, 'System32', 'taskkill.exe');
      spawnSync(taskkill, ['/pid', String(rProcess.pid), '/T', '/F'], { windowsHide: true });
    } else {
      try { process.kill(-rProcess.pid, 'SIGTERM'); } catch (error) { /* already gone */ }
    }
  });
}

function waitForShiny(port, timeoutMs = 120000) {
  const started = Date.now();
  return new Promise((resolve, reject) => {
    const attempt = () => {
      if (rExited) {
        reject(new Error('R exited before the app answered'));
        return;
      }
      const retry = () => {
        if (Date.now() - started > timeoutMs) reject(new Error('the app did not answer in time'));
        else setTimeout(attempt, 250);
      };
      http.get({ host: '127.0.0.1', port, path: '/', timeout: 2000 }, (response) => {
        response.resume();
        if (response.statusCode === 200) resolve();
        else retry();
      }).on('error', retry).on('timeout', function onTimeout() { this.destroy(); });
    };
    attempt();
  });
}

// Resolves with the text Shiny's server rendered into the page (empty after 30 s without one).
const RENDERED_TEXT = `new Promise((resolve) => {
  const started = Date.now();
  const tick = () => {
    const node = document.getElementById('runtime');
    const text = node ? node.textContent.trim() : '';
    if (text || Date.now() - started > 30000) resolve(text);
    else setTimeout(tick, 200);
  };
  tick();
})`;

function createWindow(port) {
  const window = new BrowserWindow({ width: 1000, height: 700 });
  window.webContents.on('did-fail-load', (event, code, description, url) => {
    record('did-fail-load', { code, description, url });
  });
  window.webContents.on('did-finish-load', async () => {
    record('did-finish-load', { url: window.webContents.getURL(), title: window.webContents.getTitle() });
    const text = await window.webContents.executeJavaScript(RENDERED_TEXT);
    record('shiny-session', { text });
    if (closeFile) {
      const watcher = setInterval(() => {
        if (fs.existsSync(closeFile)) {
          clearInterval(watcher);
          window.close();
        }
      }, 300);
    }
  });
  window.loadURL(`http://127.0.0.1:${port}/`);
}

app.on('window-all-closed', async () => {
  record('window-closed');
  await stopR();
  record('quit', { rExited });
  app.quit();
});

app.whenReady().then(async () => {
  record('start', {
    resources,
    runtime,
    packaged: app.isPackaged,
    platform: process.platform,
    electron: process.versions.electron,
  });
  try {
    const port = await freePort();
    startR(port);
    await waitForShiny(port);
    createWindow(port);
  } catch (error) {
    record('error', { message: String(error) });
    await stopR();
    app.exit(1);
  }
});
