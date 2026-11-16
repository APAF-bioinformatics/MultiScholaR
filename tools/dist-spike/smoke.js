'use strict';
// KIT-130 packaging spike: smoke test of the Electron application, run with the runner's Node
// (test harness only; the application itself uses nothing from the runner).
//
// Usage: node smoke.js <label> <application executable> [arguments...]
//
// It starts the application, waits until the window has loaded the Shiny page and a live
// session has rendered, fetches the page over HTTP, checks that the R serving it is the bundled
// 4.6.1 inside the application's own resources, closes the window, and checks that the
// application exits and that R is gone with it.

const { spawn } = require('child_process');
const fs = require('fs');
const http = require('http');
const os = require('os');
const path = require('path');

const [label, executable, ...appArguments] = process.argv.slice(2);
if (!label || !executable) {
  console.error('usage: node smoke.js <label> <application executable> [arguments...]');
  process.exit(2);
}

const workDir = fs.mkdtempSync(path.join(os.tmpdir(), 'kit130-smoke-'));
const logFile = path.join(workDir, 'events.jsonl');
const closeFile = path.join(workDir, 'close');
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));

function events() {
  if (!fs.existsSync(logFile)) return [];
  return fs.readFileSync(logFile, 'utf8').split('\n').filter(Boolean).map((line) => JSON.parse(line));
}

function fetchPage(url) {
  return new Promise((resolve) => {
    http.get(url, { timeout: 5000 }, (response) => {
      let body = '';
      response.on('data', (chunk) => { body += chunk; });
      response.on('end', () => resolve({ status: response.statusCode, body }));
    }).on('error', (error) => resolve({ status: null, error: String(error) }))
      .on('timeout', function onTimeout() { this.destroy(); });
  });
}

function isAlive(pid) {
  try {
    process.kill(pid, 0);
    return true;
  } catch (error) {
    return error.code === 'EPERM';
  }
}

// Same folder, whatever the spelling: forward slashes, and case-insensitive on Windows.
function canonical(p) {
  const slashed = p.replace(/\\/g, '/');
  return process.platform === 'win32' ? slashed.toLowerCase() : slashed;
}

async function main() {
  const failures = [];
  const check = (name, passed, detail = '') => {
    console.log(`${passed ? 'ok  ' : 'FAIL'} ${name}${detail ? `: ${detail}` : ''}`);
    if (!passed) failures.push(name);
  };

  const started = Date.now();
  let exitCode;
  const application = spawn(executable, appArguments, {
    env: { ...process.env, SPIKE_LOG: logFile, SPIKE_CLOSE_FILE: closeFile },
    stdio: 'inherit',
  });
  application.once('exit', (code, signal) => { exitCode = code === null ? `signal ${signal}` : code; });
  application.once('error', (error) => { exitCode = `spawn error: ${error}`; });

  // 1. The window loads the page and a live Shiny session renders.
  const find = (name) => events().find((event) => event.event === name);
  while (!find('shiny-session') && !find('error') && exitCode === undefined && Date.now() - started < 300000) {
    await sleep(500);
  }
  const startUp = find('start') || {};
  const rStarted = find('r-started') || {};
  const loaded = find('did-finish-load');
  const session = find('shiny-session');
  check('the application started R', Boolean(rStarted.pid), JSON.stringify(rStarted));
  check('the window finished loading the Shiny page', Boolean(loaded), JSON.stringify(loaded || find('error') || find('did-fail-load') || { exitCode }));
  check('a live Shiny session rendered in the window', Boolean(session && session.text), `${((Date.now() - started) / 1000).toFixed(1)} s after launch`);

  const reported = Object.fromEntries(
    ((session && session.text) || '').split('\n').filter((line) => line.includes('=')).map((line) => {
      const at = line.indexOf('=');
      return [line.slice(0, at), line.slice(at + 1)];
    }),
  );
  console.log('reported by the page:', JSON.stringify(reported));
  console.log('application resources:', startUp.resources);

  // 2. The page answers over HTTP on the port the application chose.
  const url = loaded ? loaded.url : `http://127.0.0.1:${rStarted.port}/`;
  const page = await fetchPage(url);
  check('the port serves the page', page.status === 200 && page.body.includes('hello from the bundled R'), `${url} -> ${page.status}`);

  // 3. The R behind it is the bundled one.
  const runtimeRoot = `${canonical(startUp.runtime || 'missing')}/`;
  const libPaths = (reported.lib_paths || '').split(';').filter(Boolean);
  check('R is 4.6.1', reported.r_version === '4.6.1', reported.r_version);
  check('R home is inside the application', canonical(reported.r_home || '').startsWith(runtimeRoot), reported.r_home);
  check('every library path is inside the application', libPaths.length > 0 && libPaths.every((p) => canonical(p).startsWith(runtimeRoot)), reported.lib_paths);
  const rPids = [...new Set([rStarted.pid, Number(reported.r_pid)].filter(Boolean))];
  check('R is running while the window is open', rPids.length > 0 && rPids.every(isAlive), `pids ${rPids.join(', ')}`);

  // 4. Closing the window ends the application and R.
  fs.writeFileSync(closeFile, 'close');
  const closing = Date.now();
  while (exitCode === undefined && Date.now() - closing < 60000) {
    await sleep(250);
  }
  check('the application exited after the window closed', exitCode === 0, `exit ${exitCode}`);
  if (exitCode === undefined) application.kill();
  let survivors = rPids.filter(isAlive);
  for (let waited = 0; survivors.length > 0 && waited < 5000; waited += 250) {
    await sleep(250);
    survivors = rPids.filter(isAlive);
  }
  check('R is gone after the window closed', survivors.length === 0, survivors.length ? `still alive: ${survivors.join(', ')}` : `pids ${rPids.join(', ')}`);
  const after = await fetchPage(url);
  check('the port no longer answers', after.status === null, String(after.error || after.status));

  console.log('--- application events');
  for (const event of events()) console.log(JSON.stringify(event));
  if (failures.length > 0) {
    console.log(`KIT-130 CHECK FAILED: ${label}: ${failures.join('; ')}`);
    process.exit(1);
  }
  console.log(`KIT-130 CHECK PASSED: ${label}`);
}

main().catch((error) => {
  console.error(error);
  process.exit(1);
});
