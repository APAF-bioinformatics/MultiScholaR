const assert = require('node:assert/strict');
const fs = require('fs');
const os = require('os');
const path = require('path');
const test = require('node:test');
const { execFileSync } = require('child_process');

const LOOP_SCRIPT = path.resolve(process.cwd(), 'tools/refactor/stabilization-loop.py');
const STATUS_SCRIPT = path.resolve(process.cwd(), 'tools/refactor/stabilization-status.py');
const REVIEWER_SCRIPT = path.resolve(process.cwd(), 'tools/refactor/stabilization-reviewer.py');
const SUPERVISOR_SCRIPT = path.resolve(process.cwd(), 'tools/refactor/stabilization-supervisor.py');
const NODE_BIN = process.execPath;
const PYTHON_BIN = process.env.PYTHON || '/usr/bin/python3';

function makeSandbox(prefix) {
  return fs.mkdtempSync(path.join(os.tmpdir(), prefix));
}

function writeFile(filePath, content, mode) {
  fs.mkdirSync(path.dirname(filePath), { recursive: true });
  fs.writeFileSync(filePath, content, 'utf8');
  if (mode) {
    fs.chmodSync(filePath, mode);
  }
}

function runJson(command, args, options = {}) {
  const executable = command === 'node'
    ? NODE_BIN
    : command === 'python3'
      ? PYTHON_BIN
      : command;
  let stdout;
  try {
    stdout = execFileSync(executable, args, {
      encoding: 'utf8',
      cwd: options.cwd || process.cwd(),
      env: options.env || process.env
    });
  } catch (error) {
    if (error && error.status === 0 && typeof error.stdout === 'string') {
      stdout = error.stdout;
    } else {
      throw error;
    }
  }
  return JSON.parse(stdout);
}

function makeBacklog(rootDir, targetPath) {
  const backlogPath = path.join(rootDir, 'tools', 'refactor', 'GOD_MODULE_STABILIZATION_BACKLOG.md');
  writeFile(backlogPath, `# Backlog

### 1. Completed Item

Current state:

- completed in live \`R/\`

### 2. Active Target

- Files:
  - [R/mod_example.R](${targetPath}:1) \`100\`

Current state:

- next active target
`);
  return backlogPath;
}

function makeTargetFiles(rootDir) {
  const target = path.join(rootDir, 'R', 'mod_example.R');
  const helper = path.join(rootDir, 'R', 'mod_example_helpers.R');
  const builder = path.join(rootDir, 'R', 'mod_example_builder.R');
  writeFile(target, `mod_example_server <- function(id) {\n  shiny::moduleServer(id, function(input, output, session) {\n    shiny::observeEvent(input$go, {\n      invisible(NULL)\n    })\n  })\n}\n`);
  writeFile(helper, `helperOne <- function(x) {\n  x\n}\n\nhelperTwo <- function(x) {\n  x\n}\n`);
  writeFile(builder, `mod_example_builder_server <- function(id) {\n  shiny::moduleServer(id, function(input, output, session) {\n    invisible(NULL)\n  })\n}\n`);
  return { target, helper, builder };
}

function makeTargetOverrides(rootDir, itemId, override) {
  const overridesPath = path.join(rootDir, 'tools', 'refactor', 'stabilization-target-overrides.json');
  writeFile(overridesPath, JSON.stringify({
    items: {
      [itemId]: override
    }
  }, null, 2));
  return overridesPath;
}

function makeHandover(rootDir, lines) {
  const handoverPath = path.join(rootDir, 'tools', 'refactor', 'HANDOVER-example.md');
  writeFile(handoverPath, lines.join('\n') + '\n');
  return handoverPath;
}

function makeExecutorRunner(rootDir) {
  const runner = path.join(rootDir, 'fake-executor.cjs');
  writeFile(runner, `#!/usr/bin/env node
const fs = require('fs');
const payload = JSON.parse(fs.readFileSync(process.argv[2], 'utf8'));
const state = JSON.parse(fs.readFileSync(payload.statePath, 'utf8'));
const first = state.iteration === 0;
process.stdout.write(JSON.stringify({
  status: 'completed',
  target: payload.target,
  checkpoint: first ? 'seam-1' : 'seam-2',
  targetStatusAfterIteration: first ? 'in_progress' : 'done',
  summary: first ? 'first checkpoint' : 'second checkpoint',
  verification: {
    replayCommands: [
      { argv: ['python3', '-c', 'print(123)'], label: 'smoke replay', phase: 'post' }
    ],
    display: ['python3 -c "print(123)"']
  },
  filesChanged: ['R/mod_example.R', 'tools/refactor/GOD_MODULE_STABILIZATION_BACKLOG.md'],
  notes: []
}));
`, 0o755);
  return runner;
}

function makeFailingRunner(rootDir, stderrText) {
  const runner = path.join(rootDir, 'fake-failing-executor.cjs');
  writeFile(runner, `#!/usr/bin/env node
process.stderr.write(${JSON.stringify(stderrText)});
process.exit(1);
`, 0o755);
  return runner;
}

function makeInvalidSchema(rootDir) {
  const schemaPath = path.join(rootDir, 'invalid-executor-schema.json');
  writeFile(schemaPath, JSON.stringify({
    type: 'object',
    properties: {
      status: { type: 'string' },
      verification: {
        type: 'object',
        properties: {
          replayCommands: {
            type: 'array',
            items: {
              type: 'object',
              properties: {
                argv: { type: 'array', items: { type: 'string' } },
                replayMode: { type: ['string', 'null'] }
              },
              required: ['argv'],
              additionalProperties: false
            }
          }
        },
        required: ['replayCommands'],
        additionalProperties: false
      }
    },
    required: ['status', 'verification'],
    additionalProperties: false
  }, null, 2));
  return schemaPath;
}

function makeDuplicateJsonRunner(rootDir) {
  const runner = path.join(rootDir, 'fake-duplicate-json-executor.cjs');
  writeFile(runner, `#!/usr/bin/env node
const payload = JSON.parse(require('fs').readFileSync(process.argv[2], 'utf8'));
const result = {
  status: 'completed',
  target: payload.target,
  checkpoint: 'seam-dup',
  targetStatusAfterIteration: 'in_progress',
  summary: 'duplicate json payload',
  verification: {
    replayCommands: [
      { argv: ['python3', '-c', 'print(123)'], label: 'smoke replay', phase: 'post' }
    ],
    display: ['python3 -c "print(123)"']
  },
  filesChanged: ['R/mod_example.R'],
  notes: []
};
process.stdout.write(JSON.stringify(result) + '\\n' + JSON.stringify(result) + '\\n');
`, 0o755);
  return runner;
}

function makeSleepingRunner(rootDir, sleepSeconds) {
  const runner = path.join(rootDir, 'fake-sleeping-executor.py');
  writeFile(runner, `#!/usr/bin/env python3
import time
time.sleep(${JSON.stringify(sleepSeconds)})
`, 0o755);
  return runner;
}

function makeFakeCodex(rootDir, logPath) {
  const binDir = path.join(rootDir, 'bin');
  fs.mkdirSync(binDir, { recursive: true });
  const script = path.join(binDir, 'codex');
  writeFile(script, `#!/usr/bin/env node
const fs = require('fs');
const path = require('path');
const args = process.argv.slice(2);
const outputIndex = args.indexOf('--output-last-message');
const schemaIndex = args.indexOf('--output-schema');
const outputPath = outputIndex >= 0 ? args[outputIndex + 1] : null;
const schemaPath = schemaIndex >= 0 ? args[schemaIndex + 1] : null;
const codexHome = process.env.CODEX_HOME || null;
const skillsPath = codexHome ? path.join(codexHome, 'skills') : null;
const sessionsPath = codexHome ? path.join(codexHome, 'sessions') : null;
const authPath = codexHome ? path.join(codexHome, 'auth.json') : null;
fs.writeFileSync(${JSON.stringify(logPath)}, JSON.stringify({
  cwd: process.cwd(),
  tmpdir: process.env.TMPDIR,
  codexHome,
  skillsPath,
  skillsExists: skillsPath ? fs.existsSync(skillsPath) : false,
  skillsIsSymlink: skillsPath ? fs.lstatSync(skillsPath).isSymbolicLink() : false,
  sessionsPath,
  sessionsExists: sessionsPath ? fs.existsSync(sessionsPath) : false,
  authPath,
  authExists: authPath ? fs.existsSync(authPath) : false,
  authIsSymlink: authPath ? fs.lstatSync(authPath).isSymbolicLink() : false,
  outputPath,
  schemaPath,
  args
}, null, 2));
if (outputPath) {
  fs.writeFileSync(outputPath, JSON.stringify({ ok: true }));
}
`, 0o755);
  return script;
}

function makeStatusScript(rootDir, snapshot) {
  const script = path.join(rootDir, 'fake-status.py');
  writeFile(script, `#!/usr/bin/env python3
import json
print(json.dumps(${JSON.stringify(snapshot)}))
`, 0o755);
  return script;
}

function makeBlockedLoopScript(rootDir, target) {
  const script = path.join(rootDir, 'tools', 'refactor', 'stabilization-loop.py');
  writeFile(script, `#!/usr/bin/env python3
import json
import sys

cmd = sys.argv[1]
if cmd == "status":
    print(json.dumps({
        "ok": True,
        "summary": {
            "schemaVersion": 2,
            "status": "blocked",
            "runtimeStatus": "idle",
            "iteration": 5,
            "maxIterations": 50,
            "currentItemId": "blocked-item",
            "itemStatus": {"pending": 0, "in_progress": 0, "done": 0, "blocked": 1, "failed": 0},
            "totalItems": 1
        },
        "runtime": {
            "status": "idle",
            "reason": "no_active_run",
            "loopPid": None,
            "childPid": None,
            "phase": None,
            "heartbeatAt": None,
            "heartbeatAgeSeconds": None
        },
        "currentItem": {
            "id": "blocked-item",
            "bucketNumber": 5,
            "title": "Blocked Item",
            "targetPath": ${JSON.stringify(target)},
            "status": "blocked"
        }
    }))
elif cmd == "run":
    print(json.dumps({"ok": True, "returncode": 0}))
else:
    raise SystemExit(1)
`, 0o755);
  return script;
}

test('stabilization status reports function distribution estimates', () => {
  const sandbox = makeSandbox('godmod-status-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);

  const result = runJson('python3', [
    STATUS_SCRIPT,
    '--target',
    target,
    '--backlog',
    backlogPath,
    '--json'
  ]);

  assert.equal(result.target.top_level_functions, 1);
  assert.equal(result.family.refactored_top_level_functions, 3);
  assert.equal(result.progress.target_refactored_pct_estimate, 75.0);
  assert.equal(result.backlog.repo.done, 1);
  assert.equal(result.backlog.repo.in_progress, 1);
});

test('stabilization loop resumes the same active target until reviewer-approved done', () => {
  const sandbox = makeSandbox('godmod-loop-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  const first = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  assert.equal(first.ok, true);
  let state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  let activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.iteration, 1);
  assert.ok(activeItem);
  assert.equal(activeItem.status, 'in_progress');
  assert.equal(activeItem.lastCheckpoint, 'seam-1');
  assert.equal(activeItem.lastReview.status, 'approved');

  const status = runJson('python3', [
    LOOP_SCRIPT,
    'status',
    '--state-path',
    statePath
  ]);

  assert.equal(status.ok, true);
  assert.equal(status.summary.iteration, 1);
  assert.equal(status.currentItem.targetPath, target);
  assert.equal(status.currentItem.lastCheckpoint, 'seam-1');

  const second = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  assert.equal(second.ok, true);
  state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.iteration, 2);
  assert.ok(activeItem);
  assert.equal(activeItem.status, 'done');
  assert.equal(activeItem.lastCheckpoint, 'seam-2');
});

test('stabilization loop target override can reopen a done backlog item', () => {
  const sandbox = makeSandbox('godmod-loop-override-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '2'
  ]);

  let state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  let activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(activeItem.status, 'done');

  const rerun = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1',
    '--target',
    target
  ]);

  assert.equal(rerun.ok, true);
  state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(activeItem.status, 'done');
  assert.equal(activeItem.lastCheckpoint, 'seam-2');
  assert.equal(state.iteration, 3);
  assert.equal(state.history[state.history.length - 1].targetPath, target);
});

test('stabilization loop target override can queue a manual target outside the backlog', () => {
  const sandbox = makeSandbox('godmod-loop-manual-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');
  const manualTarget = path.join(sandbox, 'R', 'manual_example.R');
  writeFile(manualTarget, 'manual_example <- function() NULL\n');

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1',
    '--target',
    manualTarget
  ]);

  assert.equal(result.ok, true);
  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.id === 'manual-manual-example');
  assert.ok(activeItem);
  assert.equal(activeItem.id, 'manual-manual-example');
  assert.equal(path.basename(activeItem.targetPath), 'manual_example.R');
  assert.equal(activeItem.status, 'in_progress');
});

test('stabilization loop applies target overrides and emits a compact work-target prompt', () => {
  const sandbox = makeSandbox('godmod-loop-overrides-');
  const { target, builder } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');
  const handoverPath = makeHandover(sandbox, [
    '# Example Handover',
    '',
    '## Current Position',
    '- target: mod_example.R',
    '- next stop point: builder helper seam in mod_example_builder.R',
    '',
    '## History',
    '- old checkpoint 1',
    '- old checkpoint 2'
  ]);

  makeTargetOverrides(sandbox, '2-active-target', {
    handoverPath,
    workTargetPath: builder,
    focusedGateCommands: [['Rscript', 'tools/test_with_renv.R', 'tests/testthat/test-prot-04-design.R']],
    promptHints: ['Prefer structural seams in mod_example_builder.R over wrapper-only review checkpoints.']
  });

  const first = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  assert.equal(first.ok, true);
  const second = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  assert.equal(second.ok, true);
  const runtimeDir = path.join(sandbox, '.god-module-stabilization');
  const payloadCandidates = fs.readdirSync(runtimeDir)
    .filter((entry) => /^2-active-target-iter-\d+-executor-payload\.json$/.test(entry))
    .sort();
  assert.ok(payloadCandidates.length > 0);
  const payloadPath = path.join(runtimeDir, payloadCandidates[payloadCandidates.length - 1]);
  const payload = JSON.parse(fs.readFileSync(payloadPath, 'utf8'));
  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const item = state.items.find((entry) => entry.id === '2-active-target');

  assert.match(item.handoverPath, /HANDOVER-example\.md$/);
  assert.match(item.workTargetPath, /mod_example_builder\.R$/);
  assert.match(payload.target, /mod_example_builder\.R$/);
  assert.match(payload.bucketTarget, /mod_example\.R$/);
  assert.match(payload.handoverPath, /HANDOVER-example\.md$/);
  assert.match(payload.prompt, /Primary work file: .*mod_example_builder\.R/);
  assert.match(payload.prompt, /Target handover: .*HANDOVER-example\.md/);
  assert.match(payload.prompt, /Latest checkpoint:/);
  assert.match(payload.prompt, /Prefer structural seams in mod_example_builder\.R/);
  assert.match(payload.prompt, /Current handover excerpt:/);
});

test('stabilization loop status returns a safe error when no state file exists', () => {
  const sandbox = makeSandbox('godmod-loop-status-');
  const statePath = path.join(sandbox, '.god-module-stabilization', 'missing-state.json');

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'status',
    '--state-path',
    statePath
  ]);

  assert.equal(result.ok, false);
  assert.match(result.error, /state file not found/);
});

test('stabilization loop status distinguishes idle runtime from in-progress workflow state', () => {
  const sandbox = makeSandbox('godmod-loop-idle-runtime-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  const status = runJson('python3', [
    LOOP_SCRIPT,
    'status',
    '--state-path',
    statePath
  ]);

  assert.equal(status.ok, true);
  assert.equal(status.summary.status, 'in_progress');
  assert.equal(status.runtime.status, 'idle');
  assert.equal(status.runtime.reason, 'no_active_run');
  assert.equal(status.currentItem.targetPath, target);
});

test('stabilization loop times out an overlong executor turn into a blocked hold', () => {
  const sandbox = makeSandbox('godmod-loop-timeout-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeSleepingRunner(sandbox, 5);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1',
    '--executor-timeout-ms',
    '200'
  ]);

  assert.equal(result.ok, true);
  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.status, 'blocked');
  assert.equal(activeItem.status, 'blocked');
  assert.equal(activeItem.lastReasonCode, 'blocked_executor_timeout');
  assert.match(activeItem.lastSummary, /executor timed out/i);
  assert.equal(state.activeRun, null);
});

test('stabilization loop status reconciles a stale active run', () => {
  const sandbox = makeSandbox('godmod-loop-stale-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  const staleStdoutPath = path.join(sandbox, '.god-module-stabilization', 'stale-executor.stdout.log');
  const staleStderrPath = path.join(sandbox, '.god-module-stabilization', 'stale-executor.stderr.log');
  const staleReviewerStdoutPath = path.join(sandbox, '.god-module-stabilization', 'stale-reviewer.stdout.log');
  const staleReviewerStderrPath = path.join(sandbox, '.god-module-stabilization', 'stale-reviewer.stderr.log');
  writeFile(staleStdoutPath, '');
  writeFile(staleStderrPath, 'executor disappeared unexpectedly');
  writeFile(staleReviewerStdoutPath, '');
  writeFile(staleReviewerStderrPath, '');

  state.status = 'in_progress';
  activeItem.status = 'in_progress';
  state.activeRun = {
    runId: 'stale-run',
    loopPid: 999999,
    itemId: activeItem.id,
    targetPath: target,
    startedAt: '2026-04-12T00:00:00+00:00',
    heartbeatAt: '2026-04-12T00:00:00+00:00',
    phase: 'executor',
    executor: {
      payloadPath: path.join(sandbox, '.god-module-stabilization', 'stale-executor-payload.json'),
      stdoutPath: staleStdoutPath,
      stderrPath: staleStderrPath,
      pid: 999998,
      startedAt: '2026-04-12T00:00:00+00:00',
      completedAt: null,
      returncode: null
    },
    reviewer: {
      payloadPath: path.join(sandbox, '.god-module-stabilization', 'stale-reviewer-payload.json'),
      stdoutPath: staleReviewerStdoutPath,
      stderrPath: staleReviewerStderrPath,
      pid: null,
      startedAt: null,
      completedAt: null,
      returncode: null
    }
  };
  fs.writeFileSync(statePath, JSON.stringify(state, null, 2));

  const status = runJson('python3', [
    LOOP_SCRIPT,
    'status',
    '--state-path',
    statePath
  ]);

  assert.equal(status.ok, true);
  assert.equal(status.summary.status, 'blocked');
  assert.equal(status.runtime.status, 'idle');
  assert.ok(status.reconciled);
  assert.equal(status.reconciled.event, 'stale_active_run_reconciled');

  const repaired = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const repairedItem = repaired.items.find((item) => item.targetPath === target);
  assert.equal(repaired.activeRun, null);
  assert.equal(repaired.status, 'blocked');
  assert.equal(repairedItem.status, 'blocked');
  assert.match(repairedItem.lastSummary, /Recovered stale active run/);
});

test('stabilization loop status trusts a recent heartbeat even when pid checks are inconclusive', () => {
  const sandbox = makeSandbox('godmod-loop-heartbeat-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  const staleStdoutPath = path.join(sandbox, '.god-module-stabilization', 'heartbeat-executor.stdout.log');
  const staleStderrPath = path.join(sandbox, '.god-module-stabilization', 'heartbeat-executor.stderr.log');
  writeFile(staleStdoutPath, '');
  writeFile(staleStderrPath, '');

  state.status = 'in_progress';
  activeItem.status = 'in_progress';
  state.activeRun = {
    runId: 'heartbeat-run',
    loopPid: 999999,
    itemId: activeItem.id,
    targetPath: target,
    startedAt: '2026-04-12T00:00:00+00:00',
    heartbeatAt: new Date().toISOString(),
    phase: 'executor',
    executor: {
      payloadPath: path.join(sandbox, '.god-module-stabilization', 'heartbeat-executor-payload.json'),
      stdoutPath: staleStdoutPath,
      stderrPath: staleStderrPath,
      pid: 999998,
      startedAt: '2026-04-12T00:00:00+00:00',
      completedAt: null,
      returncode: null
    },
    reviewer: {
      payloadPath: path.join(sandbox, '.god-module-stabilization', 'heartbeat-reviewer-payload.json'),
      stdoutPath: path.join(sandbox, '.god-module-stabilization', 'heartbeat-reviewer.stdout.log'),
      stderrPath: path.join(sandbox, '.god-module-stabilization', 'heartbeat-reviewer.stderr.log'),
      pid: null,
      startedAt: null,
      completedAt: null,
      returncode: null
    }
  };
  fs.writeFileSync(statePath, JSON.stringify(state, null, 2));

  const status = runJson('python3', [
    LOOP_SCRIPT,
    'status',
    '--state-path',
    statePath
  ]);

  assert.equal(status.ok, true);
  assert.equal(status.summary.status, 'in_progress');
  assert.equal(status.runtime.status, 'running');
  assert.equal(status.runtime.reason, 'heartbeat_recent_pid_unconfirmed');
  assert.equal(status.reconciled, null);
});

test('stabilization loop status reconciles legacy idle error state without activeRun metadata', () => {
  const sandbox = makeSandbox('godmod-loop-legacy-error-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeFailingRunner(sandbox, 'failed to lookup address information');
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  state.activeRun = null;
  state.status = 'in_progress';
  activeItem.status = 'in_progress';
  fs.writeFileSync(statePath, JSON.stringify(state, null, 2));

  const status = runJson('python3', [
    LOOP_SCRIPT,
    'status',
    '--state-path',
    statePath
  ]);

  assert.equal(status.ok, true);
  assert.equal(status.summary.status, 'blocked');
  assert.ok(status.reconciled);
  assert.equal(status.reconciled.event, 'legacy_idle_error_reconciled');

  const repaired = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const repairedItem = repaired.items.find((item) => item.targetPath === target);
  assert.equal(repaired.activeRun, null);
  assert.equal(repaired.status, 'blocked');
  assert.equal(repairedItem.status, 'blocked');
  assert.match(repairedItem.lastSummary, /Recovered legacy ambiguous loop state/);

  const recover = runJson('python3', [
    LOOP_SCRIPT,
    'recover',
    '--state-path',
    statePath
  ]);

  assert.equal(recover.ok, true);
  assert.equal(recover.recovered, false);
});

test('stabilization loop classifies network executor failures as blocked', () => {
  const sandbox = makeSandbox('godmod-loop-network-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeFailingRunner(sandbox, 'failed to lookup address information');
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  assert.equal(result.ok, true);
  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.status, 'blocked');
  assert.equal(activeItem.status, 'blocked');
  assert.match(activeItem.lastSummary, /lookup address information/);
});

test('stabilization loop classifies invalid output schema failures as blocked_protocol_schema', () => {
  const sandbox = makeSandbox('godmod-loop-bad-schema-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeFailingRunner(
    sandbox,
    'protocol schema validation failed: properties.verification.properties.replayCommands.items: required must exactly match properties (missing required keys [\'replayMode\'])'
  );
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  assert.equal(result.ok, true);
  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.status, 'blocked');
  assert.equal(activeItem.status, 'blocked');
  assert.equal(activeItem.lastReasonCode, 'blocked_protocol_schema');
});

test('stabilization codex runner rejects invalid output schema before launching codex exec', () => {
  const sandbox = makeSandbox('godmod-runner-schema-');
  const schemaPath = makeInvalidSchema(sandbox);
  const payloadPath = path.join(sandbox, 'payload.json');
  writeFile(payloadPath, JSON.stringify({
    repoRoot: sandbox,
    prompt: 'test prompt',
    sandbox: 'workspace-write',
    approvalPolicy: 'never',
    networkAccessEnabled: false
  }, null, 2));

  let error;
  try {
    execFileSync(PYTHON_BIN, [
      path.resolve(process.cwd(), 'tools/refactor/stabilization-codex-runner.py'),
      payloadPath
    ], {
      cwd: process.cwd(),
      env: {
        ...process.env,
        STABILIZATION_EXECUTOR_SCHEMA_PATH: schemaPath
      },
      encoding: 'utf8',
      stdio: ['ignore', 'pipe', 'pipe']
    });
  } catch (caught) {
    error = caught;
  }

  assert.ok(error);
  assert.notEqual(error.status, 0);
  assert.match(error.stderr, /protocol schema validation failed/i);
});

test('stabilization codex runner uses repo-local runtime roots instead of inherited TMPDIR', () => {
  const sandbox = makeSandbox('godmod-runner-tmp-root-');
  const logPath = path.join(sandbox, 'fake-codex-log.json');
  const fakeCodex = makeFakeCodex(sandbox, logPath);
  const payloadPath = path.join(sandbox, 'payload.json');
  const realCodexHome = path.join(sandbox, 'real-codex-home');
  const poisonedTmp = path.join(sandbox, 'poisoned-tmp');
  writeFile(path.join(realCodexHome, 'skills', 'god-module-stabilization', 'SKILL.md'), '# test skill\n');
  writeFile(path.join(realCodexHome, 'auth.json'), '{"token":"test"}\n');
  writeFile(path.join(realCodexHome, 'config.toml'), 'model = "gpt-5.4"\n');
  fs.mkdirSync(poisonedTmp, { recursive: true });
  writeFile(payloadPath, JSON.stringify({
    repoRoot: sandbox,
    prompt: 'test prompt',
    sandbox: 'workspace-write',
    approvalPolicy: 'never',
    networkAccessEnabled: false
  }, null, 2));

  const stdout = execFileSync(PYTHON_BIN, [
    path.resolve(process.cwd(), 'tools/refactor/stabilization-codex-runner.py'),
    payloadPath
  ], {
    cwd: process.cwd(),
    env: {
      ...process.env,
      CODEX_HOME: realCodexHome,
      PATH: `${path.dirname(fakeCodex)}${path.delimiter}${process.env.PATH}`,
      TMPDIR: poisonedTmp
    },
    encoding: 'utf8',
    stdio: ['ignore', 'pipe', 'pipe']
  });

  const log = JSON.parse(fs.readFileSync(logPath, 'utf8'));
  assert.match(log.cwd, new RegExp(`${path.basename(sandbox)}$`));
  assert.match(log.tmpdir, /\/\.god-module-stabilization\/runner-runtime\/godmod-codex-runner-[^/]+$/);
  assert.match(log.codexHome, /\/\.god-module-stabilization\/runner-runtime\/godmod-codex-runner-[^/]+\/codex-home$/);
  assert.notEqual(log.codexHome, realCodexHome);
  assert.equal(log.skillsExists, true);
  assert.equal(log.skillsIsSymlink, true);
  assert.equal(log.sessionsExists, true);
  assert.equal(log.authExists, true);
  assert.equal(log.authIsSymlink, true);
  assert.ok(!log.outputPath.startsWith(poisonedTmp));
  assert.ok(!log.schemaPath.startsWith(poisonedTmp));
  assert.deepEqual(JSON.parse(stdout), { ok: true });
});

test('stabilization loop accepts duplicate identical JSON lines from executor stdout', () => {
  const sandbox = makeSandbox('godmod-loop-dup-json-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeDuplicateJsonRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  assert.equal(result.ok, true);
  const state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.iteration, 1);
  assert.equal(activeItem.status, 'in_progress');
  assert.equal(activeItem.lastCheckpoint, 'seam-dup');
  assert.equal(activeItem.lastSummary, 'duplicate json payload');
});

test('stabilization loop target override does not clobber paused max-iteration state', () => {
  const sandbox = makeSandbox('godmod-loop-paused-override-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1'
  ]);

  let state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const activeItem = state.items.find((item) => item.targetPath === target);
  state.iteration = 1;
  state.maxIterations = 1;
  state.status = 'paused_max_iterations';
  activeItem.status = 'blocked';
  fs.writeFileSync(statePath, JSON.stringify(state, null, 2));

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1',
    '--max-iterations',
    '1',
    '--target',
    target
  ]);

  assert.equal(result.ok, true);
  state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const repairedItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.status, 'paused_max_iterations');
  assert.equal(state.currentItemId, repairedItem.id);
  assert.equal(repairedItem.status, 'blocked');
});

test('stabilization loop can raise max-iterations on an existing state file', () => {
  const sandbox = makeSandbox('godmod-loop-max-raise-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runner = makeExecutorRunner(sandbox);
  const statePath = path.join(sandbox, '.god-module-stabilization', 'loop-state.json');

  runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1',
    '--max-iterations',
    '1'
  ]);

  let state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  state.iteration = 1;
  state.maxIterations = 1;
  state.status = 'paused_max_iterations';
  fs.writeFileSync(statePath, JSON.stringify(state, null, 2));

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'run',
    '--repo-root',
    sandbox,
    '--backlog-path',
    backlogPath,
    '--state-path',
    statePath,
    '--runner-script',
    runner,
    '--iteration-limit',
    '1',
    '--max-iterations',
    '5',
    '--target',
    target
  ]);

  assert.equal(result.ok, true);
  state = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const updatedItem = state.items.find((item) => item.targetPath === target);
  assert.equal(state.maxIterations, 5);
  assert.equal(state.iteration, 2);
  assert.equal(updatedItem.status, 'done');
});

test('stabilization reviewer sanitizes annotated replay commands and reruns them directly', () => {
  const sandbox = makeSandbox('godmod-reviewer-sanitize-');
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const markerPath = path.join(sandbox, 'ran.txt');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: path.join(sandbox, 'R', 'mod_example.R'), labels: [], next_step: 'none' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: path.join(sandbox, 'backlog.md'), repo: { total: 1, done: 0, in_progress: 1, pending: 0, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 1 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });
  const payload = {
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: path.join(sandbox, 'R', 'mod_example.R'),
    backlogPath: path.join(sandbox, 'backlog.md'),
    timeoutMs: 10000,
    executorResult: {
      targetStatusAfterIteration: 'in_progress',
      filesChanged: [],
      testsRun: [
        `python3 -c 'from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("ok")' '(pre-check)'`,
        `python3 -c 'from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("ok")' '(pre-check)'`
      ]
    }
  };
  writeFile(payloadPath, JSON.stringify(payload, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(fs.readFileSync(markerPath, 'utf8'), 'ok');
  assert.equal(result.testsRun.length, 1);
  assert.match(result.testsRun[0], /python3 -c/);
  assert.ok(result.notes.some((note) => /Skipped duplicate replay command/.test(note)));
  assert.equal(result.verification.legacyFallbackUsed, true);
});

test('stabilization reviewer prefers structured replay commands over legacy strings', () => {
  const sandbox = makeSandbox('godmod-reviewer-structured-');
  const markerPath = path.join(sandbox, 'structured-ran.txt');
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: path.join(sandbox, 'R', 'mod_example.R'), labels: [], next_step: 'none' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: path.join(sandbox, 'backlog.md'), repo: { total: 1, done: 0, in_progress: 1, pending: 0, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 1 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });
  const payload = {
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: path.join(sandbox, 'R', 'mod_example.R'),
    backlogPath: path.join(sandbox, 'backlog.md'),
    timeoutMs: 10000,
    executorResult: {
      status: 'completed',
      checkpoint: 'structured-replay',
      target: path.join(sandbox, 'R', 'mod_example.R'),
      targetStatusAfterIteration: 'in_progress',
      summary: 'structured replay test',
      filesChanged: [],
      notes: [],
      verification: {
        replayCommands: [
          {
            argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("structured")`],
            label: 'structured replay',
            phase: 'post'
          }
        ],
        display: ['python3 -c <structured>']
      },
      testsRun: ['python3 -c "raise SystemExit(99)"']
    }
  };
  writeFile(payloadPath, JSON.stringify(payload, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(fs.readFileSync(markerPath, 'utf8'), 'structured');
  assert.equal(result.verification.legacyFallbackUsed, false);
  assert.equal(result.verification.replayedCommands.length, 1);
  assert.equal(result.verification.replayedCommands[0].label, 'structured replay');
});

test('stabilization reviewer reruns supplemental loop-state review commands in addition to executor replay', () => {
  const sandbox = makeSandbox('godmod-reviewer-supplemental-');
  const primaryMarkerPath = path.join(sandbox, 'primary-ran.txt');
  const supplementalMarkerPath = path.join(sandbox, 'supplemental-ran.txt');
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: path.join(sandbox, 'R', 'mod_example.R'), labels: [], next_step: 'none' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: path.join(sandbox, 'backlog.md'), repo: { total: 1, done: 0, in_progress: 1, pending: 0, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 1 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });
  const payload = {
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: path.join(sandbox, 'R', 'mod_example.R'),
    backlogPath: path.join(sandbox, 'backlog.md'),
    timeoutMs: 10000,
    supplementalReviewCommands: [
      {
        argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(supplementalMarkerPath)}).write_text("supplemental")`],
        label: 'focused gate',
        phase: 'post',
        replayMode: 'verify'
      }
    ],
    executorResult: {
      status: 'completed',
      checkpoint: 'supplemental-review',
      target: path.join(sandbox, 'R', 'mod_example.R'),
      targetStatusAfterIteration: 'in_progress',
      summary: 'supplemental reviewer command test',
      filesChanged: [],
      notes: [],
      verification: {
        replayCommands: [
          {
            argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(primaryMarkerPath)}).write_text("primary")`],
            label: 'reported gate',
            phase: 'post',
            replayMode: 'verify'
          }
        ],
        display: ['reported gate passed']
      }
    }
  };
  writeFile(payloadPath, JSON.stringify(payload, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(fs.readFileSync(primaryMarkerPath, 'utf8'), 'primary');
  assert.equal(fs.readFileSync(supplementalMarkerPath, 'utf8'), 'supplemental');
  assert.equal(result.verification.replayedCommands.length, 2);
  assert.ok(result.notes.some((note) => /supplemental reviewer command/.test(note)));
});

test('stabilization reviewer skips non-replay-safe structured apply commands', () => {
  const sandbox = makeSandbox('godmod-reviewer-skip-operate-');
  const preMarkerPath = path.join(sandbox, 'pre-ran.txt');
  const postMarkerPath = path.join(sandbox, 'post-ran.txt');
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: path.join(sandbox, 'R', 'mod_example.R'), labels: [], next_step: 'continue' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: path.join(sandbox, 'backlog.md'), repo: { total: 1, done: 0, in_progress: 1, pending: 0, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 1 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });

  const payload = {
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: path.join(sandbox, 'R', 'mod_example.R'),
    backlogPath: path.join(sandbox, 'backlog.md'),
    timeoutMs: 10000,
    executorResult: {
      status: 'completed',
      checkpoint: 'wave-apply',
      target: path.join(sandbox, 'R', 'mod_example.R'),
      targetStatusAfterIteration: 'in_progress',
      summary: 'wave apply checkpoint',
      reasonCode: 'completed_checkpoint',
      filesChanged: [],
      notes: [],
      verification: {
        replayCommands: [
          {
            argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(preMarkerPath)}).write_text("pre")`],
            label: 'pre gate',
            phase: 'pre'
          },
          {
            argv: ['python3', '-c', 'raise SystemExit(99)'],
            label: 'apply wave',
            phase: 'apply'
          },
          {
            argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(postMarkerPath)}).write_text("post")`],
            label: 'post gate',
            phase: 'post'
          }
        ],
        display: ['pre gate passed', 'wave apply succeeded', 'post gate passed']
      }
    }
  };
  writeFile(payloadPath, JSON.stringify(payload, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(fs.existsSync(preMarkerPath), false);
  assert.equal(fs.readFileSync(postMarkerPath, 'utf8'), 'post');
  assert.equal(result.verification.replayedCommands.length, 1);
  assert.ok(result.verification.replayedCommands.every((entry) => entry.replayMode === 'verify'));
  assert.ok(result.notes.some((note) => /pre-phase verification is not replay-safe/.test(note)));
  assert.ok(result.notes.some((note) => /Skipped non-replay-safe command/.test(note)));
});

test('stabilization reviewer skips verify_refactor replay after a mutating checkpoint even if labeled verify', () => {
  const sandbox = makeSandbox('godmod-reviewer-skip-verify-refactor-');
  const markerPath = path.join(sandbox, 'post-ran.txt');
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: path.join(sandbox, 'R', 'mod_example.R'), labels: [], next_step: 'continue' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: path.join(sandbox, 'backlog.md'), repo: { total: 1, done: 0, in_progress: 1, pending: 0, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 1 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });

  const payload = {
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: path.join(sandbox, 'R', 'mod_example.R'),
    backlogPath: path.join(sandbox, 'backlog.md'),
    timeoutMs: 10000,
    executorResult: {
      status: 'completed',
      checkpoint: 'wave-apply',
      target: path.join(sandbox, 'R', 'mod_example.R'),
      targetStatusAfterIteration: 'done',
      summary: 'wave apply checkpoint',
      reasonCode: 'completed_checkpoint',
      filesChanged: [],
      notes: [],
      verification: {
        replayCommands: [
          {
            argv: ['Rscript', 'tools/refactor/verify_refactor.R', '--manifest', 'tools/refactor/manifest-wave.yml'],
            label: 'manifest verification',
            phase: 'post',
            replayMode: 'verify'
          },
          {
            argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("post")`],
            label: 'focused gate',
            phase: 'post',
            replayMode: 'verify'
          }
        ],
        display: ['manifest verification passed', 'focused gate passed']
      }
    }
  };
  writeFile(payloadPath, JSON.stringify(payload, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(fs.readFileSync(markerPath, 'utf8'), 'post');
  assert.equal(result.verification.replayedCommands.length, 1);
  assert.equal(result.verification.replayedCommands[0].label, 'focused gate');
  assert.ok(result.notes.some((note) => /verify_refactor\.r/.test(note)));
});

test('stabilization reviewer ignores missing generated staging collate artifacts in filesChanged', () => {
  const sandbox = makeSandbox('godmod-reviewer-ignore-staging-collate-');
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const existingPath = path.join(sandbox, 'tools', 'refactor', 'manifest-wave.yml');
  writeFile(existingPath, 'version: 1\nentries: []\n');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: path.join(sandbox, 'R', 'mod_example.R'), labels: [], next_step: 'continue' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: path.join(sandbox, 'backlog.md'), repo: { total: 1, done: 0, in_progress: 1, pending: 0, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 1 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });

  const payload = {
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: path.join(sandbox, 'R', 'mod_example.R'),
    backlogPath: path.join(sandbox, 'backlog.md'),
    timeoutMs: 10000,
    executorResult: {
      status: 'completed',
      checkpoint: 'stage-wave',
      target: path.join(sandbox, 'R', 'mod_example.R'),
      targetStatusAfterIteration: 'in_progress',
      summary: 'stage wave checkpoint',
      reasonCode: 'completed_checkpoint',
      filesChanged: [
        existingPath,
        path.join(sandbox, 'tools', 'refactor', 'staging', 'wave-1', 'collate-wave-1.txt')
      ],
      notes: [],
      verification: {
        replayCommands: [],
        display: []
      },
      testsRun: []
    }
  };
  writeFile(payloadPath, JSON.stringify(payload, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.ok(result.notes.some((note) => /Ignored missing generated staging artifact/.test(note)));
});

test('stabilization supervisor treats blocked loop state as a hold and does not select it for relaunch', () => {
  const status = {
    summary: { status: 'blocked' },
    currentItem: {
      id: 'blocked-item',
      targetPath: '/tmp/blocked.R',
      status: 'blocked'
    }
  };

  const parsed = runJson('python3', [
    '-c',
    [
      'import importlib.util, json',
      `spec = importlib.util.spec_from_file_location("stabilization_supervisor", ${JSON.stringify(SUPERVISOR_SCRIPT)})`,
      'module = importlib.util.module_from_spec(spec)',
      'spec.loader.exec_module(module)',
      `status = json.loads(${JSON.stringify(JSON.stringify(status))})`,
      'result = {"is_blocked_hold_state": module.is_blocked_hold_state(status), "next_target_override": module.next_target_override(status)}',
      'print(json.dumps(result))'
    ].join('\n')
  ]);
  assert.equal(parsed.is_blocked_hold_state, true);
  assert.equal(parsed.next_target_override, null);
});

test('stabilization loop recover can rerun and clear a false blocked reviewer result', () => {
  const sandbox = makeSandbox('godmod-loop-recover-review-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runtimeDir = path.join(sandbox, '.god-module-stabilization');
  const statePath = path.join(runtimeDir, 'loop-state.json');
  const logPath = path.join(runtimeDir, 'loop.jsonl');
  const reviewerPayloadPath = path.join(runtimeDir, 'blocked-reviewer-payload.json');
  const markerPath = path.join(sandbox, 'recover-ran.txt');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: target, labels: [], next_step: 'continue' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: backlogPath, repo: { total: 2, done: 0, in_progress: 1, pending: 1, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 2 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });
  const executorResult = {
    status: 'completed',
    target,
    checkpoint: 'recover-checkpoint',
    targetStatusAfterIteration: 'in_progress',
    summary: 'recover blocked review',
    testsRun: [`python3 -c 'from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("recovered")' '(post-check)'`],
    filesChanged: ['R/mod_example.R'],
    notes: []
  };
  writeFile(reviewerPayloadPath, JSON.stringify({
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: target,
    backlogPath,
    timeoutMs: 10000,
    executorResult
  }, null, 2));
  writeFile(logPath, '');
  writeFile(statePath, JSON.stringify({
    schemaVersion: 2,
    createdAt: '2026-04-12T00:00:00+00:00',
    updatedAt: '2026-04-12T00:00:00+00:00',
    repoRoot: sandbox,
    backlogPath,
    statePath,
    logPath,
    stopFile: path.join(runtimeDir, 'loop.stop'),
    runtimeDir,
    status: 'blocked',
    iteration: 1,
    maxIterations: 25,
    currentItemId: 'recover-item',
    items: [
      {
        id: 'recover-item',
        bucketNumber: 2,
        title: 'Recover Item',
        targetPath: target,
        handoverPath: null,
        status: 'blocked',
        startedAt: '2026-04-12T00:00:00+00:00',
        completedAt: null,
        lastCheckpoint: 'recover-checkpoint',
        lastSummary: 'blocked by reviewer replay',
        lastProgress: null,
        lastReview: {
          status: 'blocked',
          summary: 'false blocked state',
          reasonCode: 'blocked_reviewer_replay',
          testsRun: [],
          issues: ['One or more reported gate commands failed when rerun by the reviewer.'],
          notes: [],
          testFailures: [
            {
              command: `python3 -c 'from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("recovered")' '(post-check)'`,
              returncode: 1,
              stdout: '',
              stderr: 'path does not exist'
            }
          ]
        }
      }
    ],
    history: [
      {
        iteration: 1,
        itemId: 'recover-item',
        targetPath: target,
        executor: executorResult,
        reviewer: {
          status: 'blocked',
          summary: 'false blocked state',
          reasonCode: 'blocked_reviewer_replay',
          testsRun: [],
          issues: ['One or more reported gate commands failed when rerun by the reviewer.'],
          notes: [],
          testFailures: []
        },
        progress: null,
        artifacts: {
          reviewer: {
            payloadPath: reviewerPayloadPath
          }
        },
        timestamp: '2026-04-12T00:00:00+00:00'
      }
    ],
    activeRun: null
  }, null, 2));

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'recover',
    '--state-path',
    statePath,
    '--reviewer-script',
    REVIEWER_SCRIPT
  ], { cwd: sandbox });

  assert.equal(result.ok, true);
  assert.equal(result.recovered, true);
  assert.equal(result.reviewRecovered.event, 'blocked_review_recovered');
  assert.equal(fs.readFileSync(markerPath, 'utf8'), 'recovered');

  const recoveredState = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const recoveredItem = recoveredState.items.find((item) => item.id === 'recover-item');
  assert.equal(recoveredState.status, 'in_progress');
  assert.equal(recoveredItem.status, 'in_progress');
  assert.equal(recoveredItem.lastReview.status, 'approved');
  assert.equal(recoveredState.history[recoveredState.history.length - 1].runtimeEvent, 'blocked_review_recovered');
});

test('stabilization loop recover can clear a false needs_changes reviewer result', () => {
  const sandbox = makeSandbox('godmod-loop-recover-needs-changes-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runtimeDir = path.join(sandbox, '.god-module-stabilization');
  const statePath = path.join(runtimeDir, 'loop-state.json');
  const logPath = path.join(runtimeDir, 'loop.jsonl');
  const reviewerPayloadPath = path.join(runtimeDir, 'needs-changes-reviewer-payload.json');
  const markerPath = path.join(sandbox, 'recover-needs-changes-ran.txt');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: target, labels: [], next_step: 'continue' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 1, legacy_lines: 10, family_total_top_level_functions: 1, family_total_lines: 10 },
    backlog: { path: backlogPath, repo: { total: 2, done: 0, in_progress: 1, pending: 1, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 2 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });
  const executorResult = {
    status: 'completed',
    target,
    checkpoint: 'recover-needs-changes',
    targetStatusAfterIteration: 'in_progress',
    summary: 'recover false needs_changes review',
    verification: {
      replayCommands: [
        {
          argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("approved")`],
          label: 'focused gate',
          phase: 'post',
          replayMode: 'verify'
        }
      ]
    },
    filesChanged: [path.join(sandbox, 'tools', 'refactor', 'staging', 'wave-1', 'collate-wave-1.txt')],
    notes: []
  };
  writeFile(reviewerPayloadPath, JSON.stringify({
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: target,
    backlogPath,
    timeoutMs: 10000,
    executorResult
  }, null, 2));
  writeFile(logPath, '');
  writeFile(statePath, JSON.stringify({
    schemaVersion: 2,
    createdAt: '2026-04-12T00:00:00+00:00',
    updatedAt: '2026-04-12T00:00:00+00:00',
    repoRoot: sandbox,
    backlogPath,
    statePath,
    logPath,
    stopFile: path.join(runtimeDir, 'loop.stop'),
    runtimeDir,
    status: 'blocked',
    iteration: 1,
    maxIterations: 25,
    currentItemId: 'recover-needs-changes-item',
    items: [
      {
        id: 'recover-needs-changes-item',
        bucketNumber: 2,
        title: 'Recover Needs Changes Item',
        targetPath: target,
        handoverPath: null,
        status: 'failed',
        startedAt: '2026-04-12T00:00:00+00:00',
        completedAt: null,
        lastCheckpoint: 'recover-needs-changes',
        lastSummary: 'failed by reviewer inconsistency',
        lastProgress: null,
        lastReview: {
          status: 'needs_changes',
          summary: 'false changed file mismatch',
          reasonCode: 'needs_changes_missing_changed_file',
          testsRun: [],
          issues: ['Executor reported changed files that do not exist in the repo root.'],
          notes: [],
          testFailures: []
        }
      }
    ],
    history: [
      {
        iteration: 1,
        itemId: 'recover-needs-changes-item',
        targetPath: target,
        executor: executorResult,
        reviewer: {
          status: 'needs_changes',
          summary: 'false changed file mismatch',
          reasonCode: 'needs_changes_missing_changed_file',
          testsRun: [],
          issues: ['Executor reported changed files that do not exist in the repo root.'],
          notes: [],
          testFailures: []
        },
        progress: null,
        artifacts: {
          reviewer: {
            payloadPath: reviewerPayloadPath
          }
        },
        timestamp: '2026-04-12T00:00:00+00:00'
      }
    ],
    activeRun: null
  }, null, 2));

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'recover',
    '--state-path',
    statePath,
    '--reviewer-script',
    REVIEWER_SCRIPT
  ], { cwd: sandbox });

  assert.equal(result.ok, true);
  assert.equal(result.recovered, true);
  assert.equal(result.reviewRecovered.event, 'blocked_review_recovered');
  assert.equal(fs.readFileSync(markerPath, 'utf8'), 'approved');

  const recoveredState = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const recoveredItem = recoveredState.items.find((item) => item.id === 'recover-needs-changes-item');
  assert.equal(recoveredState.status, 'in_progress');
  assert.equal(recoveredItem.status, 'in_progress');
  assert.equal(recoveredItem.lastReview.status, 'approved');
  assert.equal(recoveredState.history[recoveredState.history.length - 1].runtimeEvent, 'blocked_review_recovered');
});

test('reviewer downgrades an over-optimistic done claim to in_progress when seam work remains', () => {
  const sandbox = makeSandbox('godmod-reviewer-state-override-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const markerPath = path.join(sandbox, 'reviewer-override-ran.txt');
  const statusScript = makeStatusScript(sandbox, {
    target: {
      path: target,
      labels: ['high-risk-wrapper', 'needs-seam-introduction'],
      next_step: 'Introduce one top-level helper seam in the same file and rerun the focused gate.'
    },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 2, legacy_lines: 100, family_total_top_level_functions: 2, family_total_lines: 100 },
    backlog: { path: backlogPath, repo: { total: 2, done: 0, in_progress: 1, pending: 1, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 2 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });

  writeFile(payloadPath, JSON.stringify({
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: target,
    backlogPath,
    timeoutMs: 10000,
    executorResult: {
      status: 'completed',
      target,
      checkpoint: 'checkpoint-1',
      targetStatusAfterIteration: 'done',
      summary: 'executor thinks target is done',
      verification: {
        replayCommands: [
          {
            argv: ['python3', '-c', `from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("ok")`],
            label: 'focused gate',
            phase: 'post',
            replayMode: 'verify'
          }
        ],
        display: ['focused gate passed']
      },
      filesChanged: ['R/mod_example.R'],
      notes: []
    }
  }, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(result.reasonCode, 'approved_target_state_override');
  assert.equal(result.targetStateOverride, 'in_progress');
  assert.match(result.summary, /downgraded the target state to in_progress/i);
  assert.equal(fs.readFileSync(markerPath, 'utf8'), 'ok');
});

test('reviewer repairs malformed structured replay argv tails before rerun', () => {
  const sandbox = makeSandbox('godmod-reviewer-argv-tail-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const markerPath = path.join(sandbox, 'reviewer-argv-tail-ran.txt');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: target, labels: [], next_step: 'continue' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 2, legacy_lines: 100, family_total_top_level_functions: 2, family_total_lines: 100 },
    backlog: { path: backlogPath, repo: { total: 2, done: 0, in_progress: 1, pending: 1, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 2 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });

  writeFile(payloadPath, JSON.stringify({
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: target,
    backlogPath,
    timeoutMs: 10000,
    executorResult: {
      status: 'completed',
      target,
      checkpoint: 'checkpoint-argv-tail',
      targetStatusAfterIteration: 'in_progress',
      summary: 'executor leaked metadata fields into structured replay argv',
      verification: {
        replayCommands: [
          {
            argv: [
              'python3',
              '-c',
              `from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("ok")`,
              'label',
              'focused gate',
              'phase',
              'post',
              'replayMode',
              'verify',
              'argv_fixme'
            ],
            label: 'focused gate',
            phase: 'post',
            replayMode: 'verify'
          }
        ],
        display: ['focused gate passed']
      },
      filesChanged: ['R/mod_example.R'],
      notes: []
    }
  }, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(fs.readFileSync(markerPath, 'utf8'), 'ok');
  assert.ok(result.notes.some((note) => /Trimmed malformed structured replay metadata tail/.test(note)));
  assert.deepEqual(result.verification.replayedCommands[0].argv, [
    'python3',
    '-c',
    `from pathlib import Path; Path(${JSON.stringify(markerPath)}).write_text("ok")`
  ]);
});

test('reviewer repairs stale R top-level count replay commands before rerun', () => {
  const sandbox = makeSandbox('godmod-reviewer-r-top-level-count-fix-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const payloadPath = path.join(sandbox, 'reviewer-payload.json');
  const statusScript = makeStatusScript(sandbox, {
    target: { path: target, labels: [], next_step: 'continue' },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 2, legacy_lines: 100, family_total_top_level_functions: 2, family_total_lines: 100 },
    backlog: { path: backlogPath, repo: { total: 2, done: 0, in_progress: 1, pending: 1, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 2 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });

  writeFile(payloadPath, JSON.stringify({
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: target,
    backlogPath,
    timeoutMs: 10000,
    executorResult: {
      status: 'completed',
      target,
      checkpoint: 'checkpoint-r-top-level-count-fix',
      targetStatusAfterIteration: 'in_progress',
      summary: 'executor emitted a stale top-level count replay command',
      verification: {
        replayCommands: [
          {
            argv: [
              'Rscript',
              '-e',
              `exprs <- parse(file=${JSON.stringify(target)}, keep.source=TRUE); n <- sum(vapply(expr, function(expr) is.call(expr) && length(expr) >= 3 && as.character(expr[[1]]) %in% c('<-','=') && is.symbol(expr[[2]]), logical(1)), USE.NAMES = FALSE)`
            ],
            label: 'top-level count check',
            phase: 'post',
            replayMode: 'verify'
          }
        ],
        display: ['top-level count passed']
      }
    }
  }, null, 2));

  const result = runJson('python3', [REVIEWER_SCRIPT, payloadPath], { cwd: sandbox });

  assert.equal(result.status, 'approved');
  assert.equal(result.verification.replayedCommands.length, 1);
  assert.ok(result.notes.some((note) => /top-level count replay command/.test(note)));
  assert.deepEqual(result.verification.replayedCommands[0].argv, [
    'Rscript',
    '-e',
    `exprs <- parse(file=${JSON.stringify(target)}, keep.source=TRUE); n <- sum(vapply(exprs, function(expr) is.call(expr) && length(expr) >= 3 && as.character(expr[[1]]) %in% c('<-','=') && is.symbol(expr[[2]]), logical(1)), USE.NAMES = FALSE)`
  ]);
});

test('stabilization loop canonicalizes malformed structured replay argv tails', () => {
  const result = runJson('python3', [
    '-c',
    [
      'import importlib.util, json, pathlib, sys',
      `sys.path.insert(0, str(pathlib.Path(${JSON.stringify(path.dirname(LOOP_SCRIPT))}).resolve()))`,
      `spec = importlib.util.spec_from_file_location("stabilization_loop", ${JSON.stringify(LOOP_SCRIPT)})`,
      'module = importlib.util.module_from_spec(spec)',
      'sys.modules["stabilization_loop"] = module',
      'spec.loader.exec_module(module)',
      'payload = {"verification": {"replayCommands": [{"argv": ["python3", "-c", "print(123)", "label", "focused gate", "phase", "post", "replayMode", "verify", "argv_fixme"], "label": "focused gate", "phase": "post", "replayMode": "verify"}]}}',
      'verification, used_legacy = module.canonicalize_replay_commands(payload)',
      'print(json.dumps({"verification": verification, "used_legacy": used_legacy}))'
    ].join('\n')
  ]);

  assert.equal(result.used_legacy, false);
  assert.deepEqual(result.verification.replayCommands[0].argv, ['python3', '-c', 'print(123)']);
});

test('stabilization loop recover can reopen a false status mismatch block as in_progress', () => {
  const sandbox = makeSandbox('godmod-loop-recover-status-mismatch-');
  const { target } = makeTargetFiles(sandbox);
  const backlogPath = makeBacklog(sandbox, target);
  const runtimeDir = path.join(sandbox, '.god-module-stabilization');
  const statePath = path.join(runtimeDir, 'loop-state.json');
  const logPath = path.join(runtimeDir, 'loop.jsonl');
  const reviewerPayloadPath = path.join(runtimeDir, 'status-mismatch-reviewer-payload.json');
  const statusScript = makeStatusScript(sandbox, {
    target: {
      path: target,
      labels: ['high-risk-wrapper', 'needs-seam-introduction'],
      next_step: 'Introduce one top-level helper seam in the same file and rerun the focused gate.'
    },
    family: { helper_files: [], helper_file_count: 0, refactored_top_level_functions: 0, refactored_lines: 0, legacy_top_level_functions: 2, legacy_lines: 100, family_total_top_level_functions: 2, family_total_lines: 100 },
    backlog: { path: backlogPath, repo: { total: 2, done: 0, in_progress: 1, pending: 1, blocked: 0 }, proteomics: { total: 0, done: 0, in_progress: 0, pending: 0, blocked: 0 }, active_bucket: 2 },
    progress: { target_refactored_pct_estimate: 0.0, target_legacy_pct_estimate: 100.0, repo_completed_pct_estimate: 0.0, proteomics_completed_pct_estimate: 0.0 },
    estimation_basis: { target: 'test', backlog: 'test' },
    progress_line: 'GODMOD_PROGRESS test'
  });
  const executorResult = {
    status: 'completed',
    target,
    checkpoint: 'checkpoint-status-mismatch',
    targetStatusAfterIteration: 'done',
    summary: 'executor marked target done too early',
    verification: {
      replayCommands: [
        {
          argv: ['python3', '-c', 'print(123)'],
          label: 'focused gate',
          phase: 'post',
          replayMode: 'verify'
        }
      ],
      display: ['focused gate passed']
    },
    filesChanged: ['R/mod_example.R'],
    notes: []
  };
  writeFile(reviewerPayloadPath, JSON.stringify({
    repoRoot: sandbox,
    statusScriptPath: statusScript,
    targetPath: target,
    backlogPath,
    timeoutMs: 10000,
    executorResult
  }, null, 2));
  writeFile(logPath, '');
  writeFile(statePath, JSON.stringify({
    schemaVersion: 2,
    createdAt: '2026-04-12T00:00:00+00:00',
    updatedAt: '2026-04-12T00:00:00+00:00',
    repoRoot: sandbox,
    backlogPath,
    statePath,
    logPath,
    stopFile: path.join(runtimeDir, 'loop.stop'),
    runtimeDir,
    status: 'blocked',
    iteration: 1,
    maxIterations: 25,
    currentItemId: 'recover-status-mismatch-item',
    items: [
      {
        id: 'recover-status-mismatch-item',
        bucketNumber: 2,
        title: 'Recover Status Mismatch Item',
        targetPath: target,
        handoverPath: null,
        status: 'failed',
        startedAt: '2026-04-12T00:00:00+00:00',
        completedAt: null,
        lastCheckpoint: 'checkpoint-status-mismatch',
        lastSummary: 'blocked by status mismatch',
        lastProgress: null,
        lastReview: {
          status: 'needs_changes',
          summary: 'false status mismatch',
          reasonCode: 'needs_changes_status_mismatch',
          testsRun: [],
          issues: ['Executor reported target as done, but the fresh classifier still marks it as needs-seam-introduction.'],
          notes: [],
          testFailures: []
        }
      }
    ],
    history: [
      {
        iteration: 1,
        itemId: 'recover-status-mismatch-item',
        targetPath: target,
        executor: executorResult,
        reviewer: {
          status: 'needs_changes',
          summary: 'false status mismatch',
          reasonCode: 'needs_changes_status_mismatch',
          testsRun: [],
          issues: ['Executor reported target as done, but the fresh classifier still marks it as needs-seam-introduction.'],
          notes: [],
          testFailures: []
        },
        progress: null,
        artifacts: {
          reviewer: {
            payloadPath: reviewerPayloadPath
          }
        },
        timestamp: '2026-04-12T00:00:00+00:00'
      }
    ],
    activeRun: null
  }, null, 2));

  const result = runJson('python3', [
    LOOP_SCRIPT,
    'recover',
    '--state-path',
    statePath,
    '--reviewer-script',
    REVIEWER_SCRIPT
  ], { cwd: sandbox });

  assert.equal(result.ok, true);
  assert.equal(result.recovered, true);
  assert.equal(result.reviewRecovered.event, 'blocked_review_recovered');

  const recoveredState = JSON.parse(fs.readFileSync(statePath, 'utf8'));
  const recoveredItem = recoveredState.items.find((item) => item.id === 'recover-status-mismatch-item');
  assert.equal(recoveredState.status, 'in_progress');
  assert.equal(recoveredItem.status, 'in_progress');
  assert.equal(recoveredItem.lastReview.status, 'approved');
  assert.equal(recoveredItem.lastReview.targetStateOverride, 'in_progress');
  assert.equal(recoveredItem.lastReasonCode, 'approved_target_state_override');
});
