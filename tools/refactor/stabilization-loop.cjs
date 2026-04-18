#!/usr/bin/env node

const fs = require('fs');
const path = require('path');
const { execFileSync } = require('child_process');

const STATE_SCHEMA_VERSION = 1;
const EXECUTOR_SCHEMA_PATH = path.join(__dirname, 'stabilization-executor-output.schema.json');
const NODE_BIN = process.execPath;
const PYTHON_BIN = process.env.PYTHON || '/usr/bin/python3';

function usage() {
  console.error('Usage:');
  console.error('  stabilization-loop.cjs run [--repo-root <path>] [--backlog-path <path>] [--state-path <path>] [--log-path <path>] [--stop-file <path>] [--runner-script <path>] [--reviewer-script <path>] [--iteration-limit <n>] [--max-iterations <n>] [--model <name>] [--sandbox <mode>]');
  console.error('  stabilization-loop.cjs status [--state-path <path>]');
}

function parseArgs(argv) {
  const [command, ...rest] = argv;
  const options = {};
  for (let index = 0; index < rest.length; index += 1) {
    const token = rest[index];
    if (!token.startsWith('--')) {
      continue;
    }
    const key = token.slice(2);
    const next = rest[index + 1];
    if (next === undefined || next.startsWith('--')) {
      options[key] = 'true';
      continue;
    }
    options[key] = next;
    index += 1;
  }
  return { command, options };
}

function nowIso() {
  return new Date().toISOString();
}

function toPositiveInt(value, fallback) {
  const parsed = Number.parseInt(String(value || ''), 10);
  if (!Number.isInteger(parsed) || parsed <= 0) {
    return fallback;
  }
  return parsed;
}

function ensureParent(filePath) {
  fs.mkdirSync(path.dirname(filePath), { recursive: true });
}

function writeJson(filePath, value) {
  ensureParent(filePath);
  fs.writeFileSync(filePath, `${JSON.stringify(value, null, 2)}\n`, 'utf8');
}

function runCommand(executable, args, options) {
  try {
    return execFileSync(executable, args, options);
  } catch (error) {
    if (error && error.status === 0 && typeof error.stdout === 'string') {
      return error.stdout;
    }
    throw error;
  }
}

function readJson(filePath) {
  return JSON.parse(fs.readFileSync(filePath, 'utf8'));
}

function appendJsonl(filePath, value) {
  ensureParent(filePath);
  fs.appendFileSync(filePath, `${JSON.stringify(value)}\n`, 'utf8');
}

function resolveRepoRoot(options) {
  return path.resolve(options['repo-root'] || process.cwd());
}

function resolveBacklogPath(options) {
  return path.resolve(options['backlog-path'] || path.join(resolveRepoRoot(options), 'tools', 'refactor', 'GOD_MODULE_STABILIZATION_BACKLOG.md'));
}

function resolveStatePath(options) {
  return path.resolve(options['state-path'] || path.join(resolveRepoRoot(options), '.god-module-stabilization', 'loop-state.json'));
}

function resolveLogPath(options) {
  return path.resolve(options['log-path'] || path.join(resolveRepoRoot(options), '.god-module-stabilization', 'loop.jsonl'));
}

function resolveStopFile(options) {
  return path.resolve(options['stop-file'] || path.join(resolveRepoRoot(options), '.god-module-stabilization', 'loop.stop'));
}

function parseBacklog(backlogPath) {
  const content = fs.readFileSync(backlogPath, 'utf8');
  const headingRe = /^###\s+(\d+)\.\s+(.+?)\s*$/gm;
  const matches = [...content.matchAll(headingRe)];
  const buckets = [];
  for (let index = 0; index < matches.length; index += 1) {
    const match = matches[index];
    const start = match.index;
    const end = index + 1 < matches.length ? matches[index + 1].index : content.length;
    const block = content.slice(start, end);
    const filePaths = [...block.matchAll(/\((\/[^)]+\/R\/[^):]+\.R):\d+\)/g)].map((item) => item[1]);
    const handoverMatches = [...block.matchAll(/\((\/[^)]+\/tools\/refactor\/HANDOVER[^):]+\.md):\d+\)/g)].map((item) => item[1]);
    let status = 'pending';
    const lower = block.toLowerCase();
    if (lower.includes('completed in live `r/`') || lower.includes('no longer a blocker')) {
      status = 'done';
    } else if (lower.includes('next active target')) {
      status = 'in_progress';
    }
    buckets.push({
      bucketNumber: Number.parseInt(match[1], 10),
      title: match[2].trim(),
      status,
      targetPath: filePaths[0] || null,
      handoverPath: handoverMatches[0] || null
    });
  }
  return buckets;
}

function buildInitialState(options) {
  const repoRoot = resolveRepoRoot(options);
  const backlogPath = resolveBacklogPath(options);
  const buckets = parseBacklog(backlogPath).filter((bucket) => bucket.targetPath);
  const items = buckets.map((bucket) => ({
    id: `${bucket.bucketNumber}-${bucket.title.toLowerCase().replace(/[^a-z0-9]+/g, '-')}`,
    bucketNumber: bucket.bucketNumber,
    title: bucket.title,
    targetPath: bucket.targetPath,
    handoverPath: bucket.handoverPath,
    status: bucket.status,
    startedAt: bucket.status === 'in_progress' ? nowIso() : null,
    completedAt: bucket.status === 'done' ? nowIso() : null,
    lastCheckpoint: null,
    lastSummary: null,
    lastProgress: null,
    lastReview: null
  }));
  const firstOpenItem = items.find((item) => item.status === 'in_progress') || items.find((item) => item.status === 'pending') || null;
  return {
    schemaVersion: STATE_SCHEMA_VERSION,
    createdAt: nowIso(),
    updatedAt: nowIso(),
    repoRoot,
    backlogPath,
    statePath: resolveStatePath(options),
    logPath: resolveLogPath(options),
    stopFile: resolveStopFile(options),
    status: 'in_progress',
    iteration: 0,
    maxIterations: toPositiveInt(options['max-iterations'], 25),
    currentItemId: firstOpenItem ? firstOpenItem.id : null,
    items,
    history: []
  };
}

function loadOrCreateState(options) {
  const statePath = resolveStatePath(options);
  if (fs.existsSync(statePath)) {
    return { statePath, state: readJson(statePath) };
  }
  const state = buildInitialState(options);
  writeJson(statePath, state);
  return { statePath, state };
}

function saveState(statePath, state) {
  state.updatedAt = nowIso();
  writeJson(statePath, state);
}

function summarizeState(state) {
  const counts = { pending: 0, in_progress: 0, done: 0, blocked: 0, failed: 0 };
  for (const item of state.items) {
    if (counts[item.status] !== undefined) {
      counts[item.status] += 1;
    }
  }
  return {
    schemaVersion: state.schemaVersion,
    status: state.status,
    iteration: state.iteration,
    maxIterations: state.maxIterations,
    currentItemId: state.currentItemId,
    itemStatus: counts,
    totalItems: state.items.length
  };
}

function nextOpenItem(state) {
  return state.items.find((item) => item.status === 'in_progress') || state.items.find((item) => item.status === 'pending') || null;
}

function runStatusSnapshot(state, item) {
  const statusScript = path.join(__dirname, 'stabilization-status.py');
  const stdout = runCommand(PYTHON_BIN, [
    statusScript,
    '--target',
    item.targetPath,
    '--backlog',
    state.backlogPath,
    '--json'
  ], {
    cwd: state.repoRoot,
    encoding: 'utf8'
  });
  return JSON.parse(stdout);
}

function buildExecutorSchema() {
  return {
    type: 'object',
    properties: {
      status: { type: 'string', enum: ['completed', 'blocked', 'failed'] },
      target: { type: 'string' },
      checkpoint: { type: 'string' },
      targetStatusAfterIteration: { type: 'string', enum: ['in_progress', 'done', 'blocked'] },
      summary: { type: 'string' },
      testsRun: { type: 'array', items: { type: 'string' } },
      filesChanged: { type: 'array', items: { type: 'string' } },
      notes: { type: 'array', items: { type: 'string' } }
    },
    required: ['status', 'target', 'checkpoint', 'targetStatusAfterIteration', 'summary', 'testsRun', 'filesChanged', 'notes'],
    additionalProperties: false
  };
}

function ensureExecutorSchema() {
  if (!fs.existsSync(EXECUTOR_SCHEMA_PATH)) {
    writeJson(EXECUTOR_SCHEMA_PATH, buildExecutorSchema());
  }
}

function buildExecutorPrompt(state, item, statePath) {
  return [
    '$god-module-stabilization',
    '',
    'Continue the installed god-module-stabilization skill as one loop-driven stabilization iteration.',
    'Do exactly one bounded checkpoint on the target below and then stop.',
    '',
    `Repository root: ${state.repoRoot}`,
    `Loop state path: ${statePath}`,
    `Loop log path: ${state.logPath}`,
    `Backlog path: ${state.backlogPath}`,
    `Target bucket: ${item.bucketNumber}. ${item.title}`,
    `Target file: ${item.targetPath}`,
    `Target handover: ${item.handoverPath || 'none'}`,
    '',
    'Required behavior:',
    '- Use $god-module-stabilization in stabilize mode.',
    '- Read the backlog and handover before editing.',
    '- Perform exactly one bounded seam, staged wave/apply, or equivalent clean checkpoint.',
    '- Rerun the focused gate for the target.',
    '- Update the handover and backlog if the checkpoint changes the stop point.',
    '- Stop after the checkpoint. Do not continue into a second seam.',
    '- Return structured JSON only matching the provided schema.',
    '',
    'When deciding targetStatusAfterIteration:',
    '- use "in_progress" if the target still needs more stabilization work',
    '- use "done" only if this backlog target is genuinely complete',
    '- use "blocked" if safe progress is blocked'
  ].join('\n');
}

function runExecutor(options, state, item, statePath) {
  const runnerScript = options['runner-script'] || process.env.GODMOD_LOOP_RUNNER;
  if (!runnerScript) {
    throw new Error('No executor backend configured. Pass --runner-script or set GODMOD_LOOP_RUNNER.');
  }

  const payload = {
    prompt: buildExecutorPrompt(state, item, statePath),
    repoRoot: state.repoRoot,
    statePath,
    target: item.targetPath,
    item
  };
  const payloadPath = path.join(path.dirname(statePath), `${item.id}-executor-payload.json`);
  writeJson(payloadPath, payload);
  const stdout = runCommand(NODE_BIN, [path.resolve(runnerScript), payloadPath], {
    cwd: state.repoRoot,
    encoding: 'utf8'
  });
  return JSON.parse(stdout);
}

function runReviewer(options, state, item, executorResult) {
  const reviewerScript = options['reviewer-script'] || path.join(__dirname, 'stabilization-reviewer.py');
  const payload = {
    repoRoot: state.repoRoot,
    targetPath: item.targetPath,
    backlogPath: state.backlogPath,
    statusScriptPath: path.join(__dirname, 'stabilization-status.py'),
    executorResult,
    timeoutMs: toPositiveInt(options['review-timeout-ms'], 900000)
  };
  const payloadPath = path.join(path.dirname(state.statePath || resolveStatePath(options)), `${item.id}-reviewer-payload.json`);
  writeJson(payloadPath, payload);
  const stdout = runCommand(PYTHON_BIN, [path.resolve(reviewerScript), payloadPath], {
    cwd: state.repoRoot,
    encoding: 'utf8'
  });
  return JSON.parse(stdout);
}

function updateItemStatusFromExecutor(item, executorResult, reviewerResult, progressSnapshot) {
  item.lastCheckpoint = executorResult.checkpoint;
  item.lastSummary = executorResult.summary;
  item.lastProgress = progressSnapshot;
  item.lastReview = reviewerResult;

  if (reviewerResult.status === 'approved') {
    if (executorResult.targetStatusAfterIteration === 'done') {
      item.status = 'done';
      item.completedAt = nowIso();
      return 'completed';
    }
    if (executorResult.targetStatusAfterIteration === 'blocked') {
      item.status = 'blocked';
      return 'blocked';
    }
    item.status = 'in_progress';
    if (!item.startedAt) {
      item.startedAt = nowIso();
    }
    return 'in_progress';
  }

  item.status = reviewerResult.status === 'blocked' ? 'blocked' : 'failed';
  return 'blocked';
}

function printProgress(progressSnapshot, reviewerStatus, checkpoint) {
  const line = `${progressSnapshot.progress_line} review=${reviewerStatus} checkpoint=${checkpoint}`;
  console.log(line);
}

function runLoop(options) {
  ensureExecutorSchema();
  const { statePath, state } = loadOrCreateState(options);
  state.statePath = statePath;
  const iterationLimit = toPositiveInt(options['iteration-limit'], Number.MAX_SAFE_INTEGER);
  let ranIterations = 0;

  while (ranIterations < iterationLimit) {
    if (fs.existsSync(state.stopFile)) {
      state.status = 'stopped';
      saveState(statePath, state);
      return { ok: true, statePath, summary: summarizeState(state) };
    }

    const item = nextOpenItem(state);
    if (!item) {
      state.status = 'completed';
      state.currentItemId = null;
      saveState(statePath, state);
      return { ok: true, statePath, summary: summarizeState(state) };
    }

    if (state.iteration >= state.maxIterations) {
      state.status = 'paused_max_iterations';
      state.currentItemId = item.id;
      saveState(statePath, state);
      return { ok: true, statePath, summary: summarizeState(state) };
    }

    if (!item.startedAt) {
      item.startedAt = nowIso();
    }
    item.status = 'in_progress';
    state.currentItemId = item.id;
    saveState(statePath, state);

    try {
      const executorResult = runExecutor(options, state, item, statePath);
      state.iteration += 1;
      ranIterations += 1;

      let reviewerResult = {
        status: 'blocked',
        summary: 'Reviewer did not run.',
        testsRun: [],
        issues: ['Reviewer did not run.'],
        notes: []
      };
      let progressSnapshot = runStatusSnapshot(state, item);

      if (executorResult.status === 'completed') {
        reviewerResult = runReviewer(options, state, item, executorResult);
        progressSnapshot = reviewerResult.statusSnapshot || progressSnapshot;
      }

      const loopOutcome = updateItemStatusFromExecutor(item, executorResult, reviewerResult, progressSnapshot);
      state.history.push({
        iteration: state.iteration,
        itemId: item.id,
        targetPath: item.targetPath,
        executor: executorResult,
        reviewer: reviewerResult,
        progress: progressSnapshot,
        timestamp: nowIso()
      });
      appendJsonl(state.logPath, state.history[state.history.length - 1]);

      if (executorResult.status !== 'completed' || reviewerResult.status !== 'approved') {
        state.status = 'blocked';
        saveState(statePath, state);
        printProgress(progressSnapshot, reviewerResult.status, executorResult.checkpoint);
        return { ok: true, statePath, summary: summarizeState(state) };
      }

      state.status = loopOutcome === 'completed' && nextOpenItem(state) === null ? 'completed' : 'in_progress';
      saveState(statePath, state);
      printProgress(progressSnapshot, reviewerResult.status, executorResult.checkpoint);
    } catch (error) {
      state.iteration += 1;
      ranIterations += 1;
      state.status = 'blocked';
      const failedItem = nextOpenItem(state) || state.items.find((entry) => entry.id === state.currentItemId);
      if (failedItem) {
        failedItem.status = 'failed';
        failedItem.lastSummary = error.message;
      }
      state.history.push({
        iteration: state.iteration,
        itemId: failedItem ? failedItem.id : null,
        targetPath: failedItem ? failedItem.targetPath : null,
        executor: null,
        reviewer: null,
        progress: null,
        error: error.message,
        timestamp: nowIso()
      });
      appendJsonl(state.logPath, state.history[state.history.length - 1]);
      saveState(statePath, state);
      return { ok: true, statePath, summary: summarizeState(state) };
    }
  }

  saveState(statePath, state);
  return { ok: true, statePath, summary: summarizeState(state) };
}

function main() {
  const { command, options } = parseArgs(process.argv.slice(2));
  if (!command) {
    usage();
    process.exit(1);
  }

  if (command === 'status') {
    const statePath = resolveStatePath(options);
    const state = readJson(statePath);
    console.log(JSON.stringify({
      ok: true,
      statePath,
      summary: summarizeState(state),
      currentItem: state.items.find((item) => item.id === state.currentItemId) || null
    }, null, 2));
    return;
  }

  if (command === 'run') {
    console.log(JSON.stringify(runLoop(options), null, 2));
    return;
  }

  usage();
  process.exit(1);
}

main();
