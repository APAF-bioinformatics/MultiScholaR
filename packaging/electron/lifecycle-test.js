'use strict';
// KIT-134 P3: the lifecycle test of the packaged application (revision 2.0, 2.0.6). From the
// KIT-130 spike's smoke.js. Run with the runner's Node (test harness only; the application itself
// uses nothing from the runner).
//
// Usage: node lifecycle-test.js <label> <application executable> <application folder>
//            --expect-launcher <in-process|process|systemd> [--expect-reason <text>] [--network-closed]
//            [--expect-test-features <on|off>] [--expect-package-commit <sha>] [--expect-run-as-node-variable <value>]
//
// <application folder> is the installed or unpacked application; every process whose executable
// lies under it counts as "from the bundle", and those under its resources/runtime as R.
// The test: (1) start; the window loads the application and a live Shiny session renders the first
// tab; the launcher line says what was expected (and, with --expect-reason, its reason contains
// that text); the port answers 200. (2) A second start joins
// the first: it exits, the first window gets the event, and no second R starts. (3) The window is
// closed: the application exits 0, then no process from the bundle is left and the port refuses.
// (4) Started again and the shell killed hard: R is gone within 10 s. With --network-closed it
// first checks that this process cannot reach the network either (the caller closed it).
// KIT-176: with --expect-package-commit <sha>, the VERSION.json the application reads (beside the runtime folder)
// names that package commit, its package_on_release_line is a JSON boolean, and the About lines the application
// recorded at its start are the ones that VERSION.json gives. With --expect-run-as-node-variable <value>, the
// application's start shows that ELECTRON_RUN_AS_NODE reached it with that value.

const { spawn, execFileSync } = require('child_process');
const fs = require('fs');
const http = require('http');
const net = require('net');
const os = require('os');
const path = require('path');

const argv = process.argv.slice(2);
const flag = (name) => argv.includes(name);
const option = (name) => (argv.includes(name) ? argv[argv.indexOf(name) + 1] : undefined);
const [label, executable, folder] = argv;
const expectLauncher = option('--expect-launcher');
const expectReason = option('--expect-reason');
// KIT-156: --expect-test-features on|off: whether the test setting (MULTISCHOLAR_TEST_FEATURES=1) is expected, as R
// reports it (on: the RUV design panel only). KIT-240: the served page's Integration tile is live either way (the
// package's default).
const expectTestFeatures = option('--expect-test-features');
const expectPackageCommit = option('--expect-package-commit');
const expectRunAsNode = option('--expect-run-as-node-variable');
if (!label || !executable || !folder || !expectLauncher) {
  console.error('usage: node lifecycle-test.js <label> <executable> <application folder> --expect-launcher <launcher> [--expect-reason <text>] [--network-closed]');
  process.exit(2);
}
if (flag('--expect-package-commit') && !/^[0-9a-f]{40}$/.test(expectPackageCommit || '')) {
  console.error(`--expect-package-commit needs a 40-character hexadecimal commit, not "${expectPackageCommit}"`);
  process.exit(2);
}
const isWindows = process.platform === 'win32';
const appFolder = canonical(fs.realpathSync(folder));
// The unpacked or Windows-installed application holds resources/ directly; the Linux install keeps
// the application under app/ beside its launcher script.
const runtimeFolder = [`${appFolder}/resources/runtime`, `${appFolder}/app/resources/runtime`]
  .find((candidate) => fs.existsSync(candidate)) || `${appFolder}/resources/runtime`;
const sleep = (ms) => new Promise((resolve) => setTimeout(resolve, ms));
const failures = [];
const check = (name, passed, detail = '') => {
  console.log(`${passed ? 'ok  ' : 'FAIL'} ${name}${detail ? `: ${detail}` : ''}`);
  if (!passed) failures.push(name);
};

function canonical(p) {
  const slashed = p.replace(/\\/g, '/').replace(/\/+$/, '');
  return isWindows ? slashed.toLowerCase() : slashed;
}

// Every running process with its executable path.
function processes() {
  if (isWindows) {
    const out = execFileSync('powershell.exe', ['-NoProfile', '-Command',
      'Get-CimInstance Win32_Process | Select-Object ProcessId, ExecutablePath | ConvertTo-Json -Compress'],
    { encoding: 'utf8', maxBuffer: 64 * 1024 * 1024 });
    return JSON.parse(out).filter((p) => p.ExecutablePath).map((p) => ({ pid: p.ProcessId, exe: canonical(p.ExecutablePath) }));
  }
  const list = [];
  for (const entry of fs.readdirSync('/proc')) {
    if (!/^\d+$/.test(entry)) continue;
    try { list.push({ pid: Number(entry), exe: canonical(fs.readlinkSync(`/proc/${entry}/exe`)) }); } catch (error) { /* gone or not ours */ }
  }
  return list;
}
const fromBundle = () => processes().filter((p) => p.exe.startsWith(`${appFolder}/`));
const rProcesses = () => processes().filter((p) => p.exe.startsWith(`${runtimeFolder}/`));

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

function reachable(host, port) {
  return new Promise((resolve) => {
    const socket = net.connect({ host, port, timeout: 5000 });
    socket.once('connect', () => { socket.destroy(); resolve(true); });
    socket.once('error', () => resolve(false));
    socket.once('timeout', () => { socket.destroy(); resolve(false); });
  });
}

function launch(env) {
  const child = spawn(executable, [], { env: { ...process.env, ...env }, stdio: 'ignore' });
  const state = { child, exit: undefined };
  child.once('exit', (code, signal) => { state.exit = code === null ? `signal ${signal}` : code; });
  child.once('error', (error) => { state.exit = `spawn error: ${error}`; });
  return state;
}

// KIT-176 (D5): the About lines for a VERSION.json, the form chosen by package_on_release_line first.
function aboutLines(stamp) {
  let commitLines;
  if (stamp.package_on_release_line === false) {
    commitLines = [`Analysis code: commit ${stamp.commit} (test build, not on the release line)`];
  } else if (stamp.package_commit === stamp.commit) {
    commitLines = [`Commit: ${stamp.package_commit}`];
  } else {
    commitLines = [`Analysis code: commit ${stamp.package_commit}`, `Test build from: commit ${stamp.commit}`];
  }
  return [`MultiScholaR ${stamp.release}`, ...commitLines, `Lockfile sha256: ${stamp.lockfile_sha256 || 'unknown'}`,
    `Package snapshot: ${stamp.snapshot || 'unknown'}`, `Built: ${stamp.build_date || 'unknown'} on ${stamp.runner_image || 'unknown'}`];
}

function reader(logFile) {
  const events = () => (fs.existsSync(logFile)
    ? fs.readFileSync(logFile, 'utf8').split('\n').filter(Boolean).map((line) => JSON.parse(line)) : []);
  return { events, find: (name) => events().find((event) => event.event === name) };
}

async function waitFor(predicate, ms, step = 500) {
  const started = Date.now();
  while (Date.now() - started < ms) {
    const value = predicate();
    if (value) return value;
    await sleep(step);
  }
  return predicate();
}

function killHard(pid) {
  if (isWindows) execFileSync('taskkill.exe', ['/PID', String(pid), '/F'], { stdio: 'ignore' });
  else process.kill(pid, 'SIGKILL');
}

async function main() {
  const workDir = fs.mkdtempSync(path.join(os.tmpdir(), 'kit134-lifecycle-'));
  console.log(`${label}: ${executable}\napplication folder: ${appFolder}`);
  check('no process from the bundle before the start', fromBundle().length === 0, `${fromBundle().length} found`);

  if (flag('--network-closed')) {
    const out = await reachable('1.1.1.1', 443);
    check('the network is closed', !out, out ? 'a connection to 1.1.1.1:443 opened' : 'a connection to 1.1.1.1:443 failed');
  }

  // (1) Start: the first tab rendered by a live session.
  const first = { log: path.join(workDir, 'first.jsonl'), close: path.join(workDir, 'close-first') };
  const started = Date.now();
  const app1 = launch({ MSR_LIFECYCLE_LOG: first.log, MSR_LIFECYCLE_CLOSE: first.close });
  const log1 = reader(first.log);
  await waitFor(() => log1.find('first-tab') || log1.find('error') || app1.exit !== undefined, 360000);
  const startUp = log1.find('start') || {};
  const firstTab = log1.find('first-tab');
  const launcher = log1.find('launcher') || {};
  const answered = log1.find('application-answered') || {};
  check('the shell started', Boolean(startUp.pid), JSON.stringify({ electron: startUp.electron, noSandbox: startUp.noSandbox }));
  if (expectRunAsNode !== undefined) {
    check(`ELECTRON_RUN_AS_NODE=${expectRunAsNode} reached the application`, startUp.runAsNodeVariable === expectRunAsNode,
      JSON.stringify({ runAsNodeVariable: startUp.runAsNodeVariable === undefined ? '<absent>' : startUp.runAsNodeVariable }));
  }
  if (expectPackageCommit !== undefined) {
    const stampFile = path.join(path.dirname(runtimeFolder), 'VERSION.json');
    let stamp = null;
    try { stamp = JSON.parse(fs.readFileSync(stampFile, 'utf8')); } catch (error) { stamp = null; }
    check(`VERSION.json names the package commit ${expectPackageCommit}`, Boolean(stamp) && stamp.package_commit === expectPackageCommit,
      `${stampFile}: ${JSON.stringify(stamp ? { commit: stamp.commit, package_commit: stamp.package_commit, package_on_release_line: stamp.package_on_release_line } : null)}`);
    const onLine = Boolean(stamp) && typeof stamp.package_on_release_line === 'boolean';
    check('VERSION.json\'s package_on_release_line is a JSON boolean', onLine, JSON.stringify(stamp ? stamp.package_on_release_line : null));
    const expected = onLine ? aboutLines(stamp) : null;
    check('the About lines the application recorded are the ones VERSION.json gives',
      Boolean(expected) && JSON.stringify(startUp.about) === JSON.stringify(expected),
      `recorded ${JSON.stringify(startUp.about === undefined ? '<absent>' : startUp.about)}; expected ${JSON.stringify(expected)}`);
  }
  check('the probe ran', Boolean(log1.find('probe')), JSON.stringify(log1.find('probe') || {}));
  check('the window loaded the application and a live session rendered the first tab',
    Boolean(firstTab && firstTab.connected && firstTab.home && firstTab.rendered),
    `${JSON.stringify(firstTab || log1.find('error') || { exit: app1.exit })}; ${((Date.now() - started) / 1000).toFixed(1)} s after launch`);
  check(`the launcher is ${expectLauncher}`, launcher.launcher === expectLauncher, JSON.stringify(launcher));
  if (expectReason) {
    check(`the launcher's reason names "${expectReason}"`, String(launcher.reason || '').includes(expectReason), JSON.stringify(launcher.reason));
  }
  const url = `http://127.0.0.1:${answered.port}/`;
  const page = await fetchPage(url);
  check('the port serves the application', page.status === 200 && /MultiScholaR/.test(page.body || ''), `${url} -> ${page.status}`);
  if (expectTestFeatures) {
    const on = expectTestFeatures === 'on';
    const reported = log1.find('test-features');
    check(`the test setting is ${expectTestFeatures} as R reports it`, on
      ? Boolean(reported && reported.ruv_design_panel === true && !('integration_gui' in reported)) : !reported,
    JSON.stringify(reported || null));
    const tile = ((page.body || '').match(/<div[^>]*id="integration_box"[^>]*>/) || [''])[0];
    check('the Integration tile on the served page is live',
      Boolean(tile) && !/disabled/.test(tile), tile || 'no integration_box');
  }
  const r1 = rProcesses();
  check('R from the bundle is running', r1.length > 0, r1.map((p) => `${p.pid} ${p.exe}`).join('; '));

  // (2) A second start joins the first.
  const second = launch({ MSR_LIFECYCLE_LOG: first.log });
  await waitFor(() => second.exit !== undefined, 60000);
  check('a second start exits', second.exit === 0, `exit ${second.exit}`);
  check('the first window received the second start', Boolean(await waitFor(() => log1.find('second-instance'), 10000)));
  await sleep(2000);
  const r2 = rProcesses();
  check('no second R started', r2.length === r1.length && r2.every((p) => r1.some((q) => q.pid === p.pid)),
    `before ${r1.map((p) => p.pid).join(',')}; after ${r2.map((p) => p.pid).join(',')}`);

  // (3) Close the window.
  fs.writeFileSync(first.close, 'close');
  await waitFor(() => app1.exit !== undefined, 60000, 250);
  check('the application exited 0 after the window closed', app1.exit === 0, `exit ${app1.exit}`);
  const left = await waitFor(() => fromBundle().length === 0, 10000, 250) ? [] : fromBundle();
  check('no process from the bundle is left', left.length === 0, left.map((p) => `${p.pid} ${p.exe}`).join('; ') || 'none');
  const after = await fetchPage(url);
  check('the port refuses', after.status === null, String(after.error || after.status));

  // (4) Start again, then kill the shell hard: R must go within 10 s.
  const third = { log: path.join(workDir, 'third.jsonl') };
  const app3 = launch({ MSR_LIFECYCLE_LOG: third.log });
  const log3 = reader(third.log);
  await waitFor(() => log3.find('first-tab') || log3.find('error') || app3.exit !== undefined, 360000);
  const r3 = rProcesses();
  check('started again, R running', Boolean(log3.find('first-tab')) && r3.length > 0, `${r3.length} R processes`);
  killHard(app3.child.pid);
  const killedAt = Date.now();
  const gone = await waitFor(() => rProcesses().length === 0, 10000, 250);
  check('R is gone within 10 s of the shell being killed', Boolean(gone),
    gone ? `${((Date.now() - killedAt) / 1000).toFixed(1)} s` : rProcesses().map((p) => `${p.pid} ${p.exe}`).join('; '));
  for (const p of fromBundle()) { try { killHard(p.pid); } catch (error) { /* gone */ } }

  console.log('--- shell events (first start)');
  for (const event of log1.events()) console.log(JSON.stringify(event));
  if (failures.length > 0 && startUp.logs) {
    // The application's own logs of the first start, for the diagnosis.
    for (const name of [`r-${startUp.startedAt}.log`, `shell-${startUp.startedAt}.log`]) {
      const file = path.join(startUp.logs, name);
      const lines = fs.existsSync(file) ? fs.readFileSync(file, 'utf8').split(/\r?\n/) : ['(missing)'];
      console.log(`--- ${name}, last 60 lines`);
      console.log(lines.slice(-60).join('\n'));
    }
  }
  if (failures.length > 0) {
    console.log(`KIT-134 CHECK FAILED: ${label}: ${failures.join('; ')}`);
    process.exit(1);
  }
  console.log(`KIT-134 CHECK PASSED: ${label}`);
}

main().catch((error) => {
  console.error(error);
  process.exit(1);
});
