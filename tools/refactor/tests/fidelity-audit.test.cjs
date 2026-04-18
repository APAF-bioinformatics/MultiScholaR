const assert = require('node:assert/strict');
const fs = require('node:fs');
const os = require('node:os');
const path = require('node:path');
const { execFileSync } = require('node:child_process');

const AUDIT_SCRIPT = path.resolve(process.cwd(), 'tools/refactor/fidelity_audit.py');
const INVENTORY_SCRIPT = path.resolve(process.cwd(), 'tools/refactor/fidelity_inventory_snapshot.R');
const PYTHON_BIN = process.env.PYTHON || '/usr/bin/python3';
const RSCRIPT_BIN = process.env.RSCRIPT || 'Rscript';
const queuedTests = [];
const serialTest = (name, fn) => {
  queuedTests.push({ name, fn });
};

function makeSandbox(prefix) {
  return fs.mkdtempSync(path.join(os.tmpdir(), prefix));
}

function writeFile(filePath, content) {
  fs.mkdirSync(path.dirname(filePath), { recursive: true });
  fs.writeFileSync(filePath, content, 'utf8');
}

function hasRecoverableStdout(error) {
  return Boolean(
    error &&
    typeof error.stdout === 'string' &&
    (error.status === 0 || error.code === 'EPERM')
  );
}

function runJson(command, args, options = {}) {
  const executable = command === 'python3'
    ? PYTHON_BIN
    : command === 'Rscript'
      ? RSCRIPT_BIN
      : command;
  const argv = command === 'Rscript'
    ? ['--vanilla', ...args]
    : args;
  let stdout;
  try {
    stdout = execFileSync(executable, argv, {
      cwd: options.cwd || process.cwd(),
      env: options.env || process.env,
      encoding: 'utf8'
    });
  } catch (error) {
    if (hasRecoverableStdout(error)) {
      stdout = error.stdout;
    } else {
      throw error;
    }
  }
  return JSON.parse(stdout);
}

function runText(command, args, options = {}) {
  const executable = command === 'python3'
    ? PYTHON_BIN
    : command === 'Rscript'
      ? RSCRIPT_BIN
      : command;
  const argv = command === 'Rscript'
    ? ['--vanilla', ...args]
    : args;
  try {
    return execFileSync(executable, argv, {
      cwd: options.cwd || process.cwd(),
      env: options.env || process.env,
      encoding: 'utf8'
    });
  } catch (error) {
    if (hasRecoverableStdout(error)) {
      return error.stdout;
    }
    throw error;
  }
}

function querySqlite(dbPath, sql) {
  let stdout;
  try {
    stdout = execFileSync(PYTHON_BIN, [
      '-c',
      'import json, sqlite3, sys; conn = sqlite3.connect(sys.argv[1]); rows = conn.execute(sys.argv[2]).fetchall(); print(json.dumps(rows))',
      dbPath,
      sql
    ], { encoding: 'utf8' });
  } catch (error) {
    if (hasRecoverableStdout(error)) {
      stdout = error.stdout;
    } else {
      throw error;
    }
  }
  return JSON.parse(stdout);
}

function makeFixtureRepo() {
  const root = makeSandbox('fidelity-audit-');
  writeFile(path.join(root, 'DESCRIPTION'), `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Demo package for fidelity audit tests.
License: MIT
Encoding: UTF-8
Collate:
    'a.R'
    'b.R'
`);
  writeFile(path.join(root, 'NAMESPACE'), `export(foo)
exportMethods(processThing)
exportClasses(DemoData)
`);
  writeFile(path.join(root, 'R', 'a.R'), `foo <- function(x) {
  x + 1
}

helper = function(y) {
  y
}

setGeneric("processThing", function(x) {
  standardGeneric("processThing")
})
`);
  writeFile(path.join(root, 'R', 'b.R'), `setClass("DemoData", representation(value = "numeric"))

setMethod("processThing", "DemoData", function(x) {
  x
})
`);
  return root;
}

function makeManifestFixtureRepo() {
  const root = makeFixtureRepo();
  writeFile(path.join(root, 'R', 'c.R'), `# BEGIN BLOCK
anchor_value <- 1
# END BLOCK
`);
  writeFile(path.join(root, 'tools', 'refactor', 'manifest-test.yml'), `version: 1
entries:
  - id: symbol_foo
    source: R/a.R
    selector: { kind: symbol, value: foo }
    target: R/extracted_symbols.R

  - id: expr_helper
    source: R/a.R
    selector: { kind: expr_index, value: 2 }
    target: R/extracted_symbols.R

  - id: method_process_thing
    source: R/b.R
    selector:
      kind: setMethod
      value: processThing
      signature: DemoData
    target: R/extracted_s4.R

  - id: anchor_block
    source: R/c.R
    selector:
      kind: anchor_range
      start: "^# BEGIN BLOCK$"
      end: "^# END BLOCK$"
      end_inclusive: true
    target: R/extracted_anchor.R
`);
  return root;
}

function makeSurfaceDiffPair() {
  const baselineRoot = makeSandbox('fidelity-surface-baseline-');
  const targetRoot = makeSandbox('fidelity-surface-target-');

  writeFile(path.join(baselineRoot, 'DESCRIPTION'), `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Baseline package for fidelity audit surface tests.
License: MIT
Encoding: UTF-8
Collate:
    'a.R'
    'b.R'
`);
  writeFile(path.join(baselineRoot, 'NAMESPACE'), `export(foo)
exportMethods(processThing)
exportClasses(DemoData)
`);
  writeFile(path.join(baselineRoot, 'R', 'a.R'), `foo <- function(x) {
  x + 1
}

helper <- function(y) {
  y
}

setGeneric("processThing", function(x) {
  standardGeneric("processThing")
})
`);
  writeFile(path.join(baselineRoot, 'R', 'b.R'), `setClass("DemoData", representation(value = "numeric"))

setMethod("processThing", "DemoData", function(x) {
  x
})
`);

  writeFile(path.join(targetRoot, 'DESCRIPTION'), `Package: DemoPkg
Title: Demo
Version: 0.0.2
Description: Target package for fidelity audit surface tests.
License: MIT
Encoding: UTF-8
Collate:
    'b.R'
    'c.R'
    'd.R'
`);
  writeFile(path.join(targetRoot, 'NAMESPACE'), `export(foo)
export(bar)
exportClasses(DemoData)
`);
  writeFile(path.join(targetRoot, 'R', 'b.R'), `foo <- function(x) {
  x + 1
}

setClass("DemoData", representation(value = "numeric"))
`);
  writeFile(path.join(targetRoot, 'R', 'c.R'), `setGeneric("processThing", function(x) {
  standardGeneric("processThing")
})

setMethod("processThing", "DemoData", function(x) {
  x
})

bar <- function(z) {
  z
}
`);
  writeFile(path.join(targetRoot, 'R', 'd.R'), `bar <- function(z) {
  z
}
`);

  return { baselineRoot, targetRoot };
}

function makeManifestComparisonPair() {
  const baselineRoot = makeSandbox('fidelity-manifest-baseline-');
  const targetRoot = makeSandbox('fidelity-manifest-target-');

  writeFile(path.join(baselineRoot, 'R', 'source_symbol.R'), `foo <- function(x) {
  x + 1
}
`);
  writeFile(path.join(baselineRoot, 'R', 'source_comment.R'), `# baseline comment
commentHelper <- function(z) {
  z
}
`);
  writeFile(path.join(baselineRoot, 'R', 'source_expr.R'), `unusedHelper <- function(a) {
  a
}

exprHelper <- function(b) {
  b
}
`);
  writeFile(path.join(baselineRoot, 'R', 'source_method.R'), `setMethod("processThing", "DemoData", function(x) {
  x
})
`);
  writeFile(path.join(baselineRoot, 'R', 'source_anchor.R'), `# BEGIN BLOCK
anchorHelper <- function(a) {
  a + 1
}
# END BLOCK
`);
  writeFile(path.join(baselineRoot, 'R', 'source_missing.R'), `missingHelper <- function(v) {
  v
}
`);
  writeFile(path.join(baselineRoot, 'R', 'source_manual.R'), `manualHelper <- function(m) {
  m
}
`);

  writeFile(path.join(targetRoot, 'R', 'target_symbol.R'), `foo <- function(x) {
  x + 1
}
`);
  writeFile(path.join(targetRoot, 'R', 'target_comment.R'), `# updated comment
commentHelper <- function(z) {
  z
}
`);
  writeFile(path.join(targetRoot, 'R', 'target_expr.R'), `exprHelper <- function(b) {
  b
}
`);
  writeFile(path.join(targetRoot, 'R', 'target_method.R'), `setMethod("processThing", "DemoData", function(x) {
  x
})
`);
  writeFile(path.join(targetRoot, 'R', 'target_anchor.R'), `# BEGIN BLOCK
anchorHelper <- function(a) {
    a + 1
}
# END BLOCK
`);
  writeFile(path.join(targetRoot, 'R', 'target_missing.R'), `otherHelper <- function(q) {
  q
}
`);
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-fidelity.yml'), `version: 1
entries:
  - id: symbol_exact
    source: R/source_symbol.R
    selector: { kind: symbol, value: foo }
    target: R/target_symbol.R
    group: fidelity.symbol

  - id: symbol_comment_only
    source: R/source_comment.R
    selector: { kind: symbol, value: commentHelper }
    target: R/target_comment.R
    group: fidelity.comment

  - id: expr_exact
    source: R/source_expr.R
    selector: { kind: expr_index, value: 2 }
    target: R/target_expr.R
    group: fidelity.expr

  - id: method_exact
    source: R/source_method.R
    selector:
      kind: setMethod
      value: processThing
      signature: DemoData
    target: R/target_method.R
    group: fidelity.method

  - id: anchor_normalized
    source: R/source_anchor.R
    selector:
      kind: anchor_range
      start: "^# BEGIN BLOCK$"
      end: "^# END BLOCK$"
      end_inclusive: true
    target: R/target_anchor.R
    group: fidelity.anchor

  - id: missing_target
    source: R/source_missing.R
    selector: { kind: symbol, value: missingHelper }
    target: R/target_missing.R
    group: fidelity.missing

  - id: manual_expected
    action: manual_merge
    source: R/source_manual.R
    selector: { kind: symbol, value: manualHelper }
    target: R/target_manual.R
    group: fidelity.manual
`);

  return { baselineRoot, targetRoot };
}

function makeManifestSourceGapPair() {
  const baselineRoot = makeSandbox('fidelity-manifest-source-gap-baseline-');
  const targetRoot = makeSandbox('fidelity-manifest-source-gap-target-');

  writeFile(path.join(baselineRoot, 'R', 'source_gap.R'), `otherHelper <- function(x) {
  x
}
`);
  writeFile(path.join(targetRoot, 'R', 'target_gap.R'), `gapHelper <- function(x) {
  x
}
`);
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-source-gap.yml'), `version: 1
entries:
  - id: source_gap
    source: R/source_gap.R
    selector: { kind: symbol, value: gapHelper }
    target: R/target_gap.R
    group: fidelity.source_gap
`);

  return { baselineRoot, targetRoot };
}

function makeBehaviorComparisonPair() {
  const baselineRoot = makeSandbox('fidelity-behavior-baseline-');
  const targetRoot = makeSandbox('fidelity-behavior-target-');

  writeFile(path.join(baselineRoot, 'R', 'exact.R'), `exactValue <- function(x) {
  x + 1L
}
`);
  writeFile(path.join(targetRoot, 'R', 'exact.R'), `exactValue <- function(x) {
  x + 1L
}
`);

  writeFile(path.join(baselineRoot, 'R', 'normalize.R'), `normalizedTable <- function() {
  data.frame(key = c("b", "a"), value = c(2L, 1L), stringsAsFactors = FALSE)
}
`);
  writeFile(path.join(targetRoot, 'R', 'normalize.R'), `normalizedTable <- function() {
  data.frame(key = c("a", "b"), value = c(1L, 2L), stringsAsFactors = FALSE)
}
`);

  writeFile(path.join(baselineRoot, 'R', 'warn.R'), `warningValue <- function() {
  warning("baseline warning")
  1L
}
`);
  writeFile(path.join(targetRoot, 'R', 'warn.R'), `warningValue <- function() {
  warning("target warning")
  1L
}
`);

  writeFile(path.join(baselineRoot, 'R', 'message.R'), `messageValue <- function() {
  message("baseline message")
  1L
}
`);
  writeFile(path.join(targetRoot, 'R', 'message.R'), `messageValue <- function() {
  message("target message")
  1L
}
`);

  writeFile(path.join(baselineRoot, 'R', 'error.R'), `errorValue <- function() {
  stop("shared error")
}
`);
  writeFile(path.join(targetRoot, 'R', 'error.R'), `errorValue <- function() {
  stop("shared error")
}
`);

  writeFile(path.join(baselineRoot, 'R', 'grid.R'), `setClass("DemoGrid", slots = list(items = "list"))

setGeneric("InitialiseDemoGrid", function(dummy = NULL) {
  standardGeneric("InitialiseDemoGrid")
})

setMethod("InitialiseDemoGrid", signature(dummy = "ANY"), function(dummy = NULL) {
  new("DemoGrid", items = list())
})
`);
  writeFile(path.join(targetRoot, 'R', 'grid.R'), `setClass("DemoGrid", slots = list(items = "list"))

setGeneric("InitialiseDemoGrid", function(dummy = NULL) {
  standardGeneric("InitialiseDemoGrid")
})

setMethod("InitialiseDemoGrid", signature(dummy = "ANY"), function(dummy = NULL) {
  new("DemoGrid", items = list())
})
`);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'behavior-cases-fixture.json'), JSON.stringify({
    cases: [
      {
        id: 'exact_scalar',
        family: 'pure_helper',
        entity_key: 'symbol::exactValue',
        fixture_kind: 'inline_setup',
        normalizer_key: 'scalar_integer',
        risk_level: 'low',
        files: ['R/exact.R'],
        setup_code: ['x <- 2L'],
        call_expr: 'exactValue(x)',
        normalize_expr: 'as.integer(result)'
      },
      {
        id: 'normalized_table',
        family: 'pure_helper',
        entity_key: 'symbol::normalizedTable',
        fixture_kind: 'inline_setup',
        normalizer_key: 'sorted_table',
        risk_level: 'medium',
        files: ['R/normalize.R'],
        call_expr: 'normalizedTable()',
        normalize_expr: 'local({ ordered <- result[order(result[[\"key\"]]), , drop = FALSE]; list(key = unname(as.character(ordered[[\"key\"]])), value = unname(as.integer(ordered[[\"value\"]]))) })'
      },
      {
        id: 'warning_diff',
        family: 'pure_helper',
        entity_key: 'symbol::warningValue',
        fixture_kind: 'inline_setup',
        normalizer_key: 'scalar_integer',
        risk_level: 'medium',
        files: ['R/warn.R'],
        call_expr: 'warningValue()',
        normalize_expr: 'as.integer(result)'
      },
      {
        id: 'message_diff',
        family: 'pure_helper',
        entity_key: 'symbol::messageValue',
        fixture_kind: 'inline_setup',
        normalizer_key: 'scalar_integer',
        risk_level: 'low',
        compare_messages: true,
        files: ['R/message.R'],
        call_expr: 'messageValue()',
        normalize_expr: 'as.integer(result)'
      },
      {
        id: 'error_match',
        family: 'pure_helper',
        entity_key: 'symbol::errorValue',
        fixture_kind: 'inline_setup',
        normalizer_key: 'identity',
        risk_level: 'high',
        files: ['R/error.R'],
        call_expr: 'errorValue()'
      },
      {
        id: 's4_grid',
        family: 's4_isolated',
        entity_key: 'setMethod::InitialiseDemoGrid::ANY',
        fixture_kind: 'inline_setup',
        normalizer_key: 'slot_lengths',
        risk_level: 'high',
        compare_warnings: false,
        files: ['R/grid.R'],
        call_expr: 'InitialiseDemoGrid()',
        normalize_expr: 'list(class = as.character(class(result)), slot_names = methods::slotNames(result), slot_lengths = as.integer(vapply(methods::slotNames(result), function(name) length(methods::slot(result, name)), integer(1))))'
      }
    ]
  }, null, 2));

  return { baselineRoot, targetRoot };
}

function makeContractCatalogRepo() {
  const targetRoot = makeSandbox('fidelity-contracts-target-');

  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-prot-01c-import-module-contracts.R'), `library(testthat)

test_that("mod_prot_import_server initializes default outputs", {
  expect_true(TRUE)
})

test_that("mod_prot_import_server reports invalid input errors", {
  expect_true(TRUE)
})

test_that("mod_prot_import_server exports reports through download seam", {
  expect_true(TRUE)
})
`);

  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-prot-07c-da-results-characterization.R'), `library(testthat)

test_that("extractResults preserves names and pulls nested results", {
  expect_equal(1, 1)
})
`);

  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-prot-07b-da-handlers-compat.R'), `library(testthat)

test_that("writeInteractiveVolcanoPlotProteomicsMain accepts legacy aliases", {
  expect_true(FALSE)
})
`);

  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-protein-da-golden-master.R'), `library(testthat)

test_that("protein_deAnalysisWrapperFunction and its alias are exported", {
  skip("golden master fixture")
})
`);

  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-prot-03-rollup.R'), `library(testthat)

test_that("rollup helper returns expected grouping", {
  expect_true(TRUE)
})
`);

  return targetRoot;
}

function makeCloseoutComparisonTriplet() {
  const baselineRoot = makeSandbox('fidelity-closeout-baseline-');
  const mainlineRoot = makeSandbox('fidelity-closeout-mainline-');
  const targetRoot = makeSandbox('fidelity-closeout-target-');

  const description = `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Closeout fixture package.
License: MIT
Encoding: UTF-8
Collate:
    'core.R'
`;
  const namespace = `export(commentHelper)
export(exactValue)
`;
  const baselineCore = `# baseline comment
commentHelper <- function(x) {
  x
}

exactValue <- function(x) {
  x + 1L
}
`;
  const targetCore = `# updated comment
commentHelper <- function(x) {
  x
}

exactValue <- function(x) {
  x + 1L
}
`;

  for (const root of [baselineRoot, mainlineRoot, targetRoot]) {
    writeFile(path.join(root, 'DESCRIPTION'), description);
    writeFile(path.join(root, 'NAMESPACE'), namespace);
  }
  writeFile(path.join(baselineRoot, 'R', 'core.R'), baselineCore);
  writeFile(path.join(mainlineRoot, 'R', 'core.R'), targetCore);
  writeFile(path.join(targetRoot, 'R', 'core.R'), targetCore);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-closeout.yml'), `version: 1
entries:
  - id: comment_helper
    source: R/core.R
    selector: { kind: symbol, value: commentHelper }
    target: R/core.R
    group: fidelity.comment
`);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'behavior-cases-closeout.json'), JSON.stringify({
    cases: [
      {
        id: 'exact_value',
        family: 'pure_helper',
        entity_key: 'symbol::exactValue',
        fixture_kind: 'inline_setup',
        normalizer_key: 'scalar_integer',
        risk_level: 'low',
        files: ['R/core.R'],
        setup_code: ['x <- 2L'],
        call_expr: 'exactValue(x)',
        normalize_expr: 'as.integer(result)'
      }
    ]
  }, null, 2));

  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-prot-03-rollup.R'), `library(testthat)

test_that("rollup helper returns expected grouping", {
  expect_true(TRUE)
})
`);

  return { baselineRoot, mainlineRoot, targetRoot };
}

serialTest('inventory snapshot captures symbols, S4 entities, exports, and collate order', () => {
  const repoRoot = makeFixtureRepo();
  const snapshot = runJson('Rscript', [INVENTORY_SCRIPT, '--repo-root', repoRoot]);

  assert.equal(snapshot.counts.r_files, 2);
  assert.equal(snapshot.counts.collate, 2);
  assert.equal(snapshot.counts.exports.functions, 1);
  assert.equal(snapshot.counts.exports.methods, 1);
  assert.equal(snapshot.counts.exports.classes, 1);
  assert.equal(snapshot.counts.parse_failures, 0);

  const byKey = new Map(snapshot.entities.map((entity) => [entity.entity_key, entity]));
  assert.equal(byKey.get('symbol::foo').exported, true);
  assert.equal(byKey.get('symbol::foo').collate_index, 1);
  assert.equal(byKey.get('symbol::helper').exported, false);
  assert.equal(byKey.get('setGeneric::processThing').exported, false);
  assert.equal(byKey.get('setClass::DemoData').exported, true);
  assert.equal(byKey.get('setMethod::processThing::DemoData').exported, true);
  assert.equal(byKey.get('setMethod::processThing::DemoData').collate_index, 2);

  assert.deepEqual(
    snapshot.collate.map((entry) => entry.file_path),
    ['R/a.R', 'R/b.R']
  );
});

serialTest('audit bootstrap and inventory import create durable state and persist entities', () => {
  const repoRoot = makeFixtureRepo();
  const bootstrap = runJson('python3', [AUDIT_SCRIPT, 'bootstrap', '--repo-root', repoRoot]);

  assert.equal(bootstrap.status, 'bootstrap_ready');
  assert.equal(fs.existsSync(path.join(repoRoot, '.refactor-fidelity-audit', 'audit.db')), true);
  assert.equal(fs.existsSync(path.join(repoRoot, '.refactor-fidelity-audit', 'events.jsonl')), true);
  assert.equal(fs.existsSync(path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-summary.json')), true);
  assert.equal(fs.existsSync(path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json')), true);

  const tables = querySqlite(
    path.join(repoRoot, '.refactor-fidelity-audit', 'audit.db'),
    "select name from sqlite_master where type = 'table' order by name"
  ).map((row) => row[0]);
  assert.deepEqual(tables, [
    'audit_runs',
    'behavior_cases',
    'behavior_results',
    'contract_results',
    'contract_test_cases',
    'contract_test_files',
    'exceptions',
    'inventory_diffs',
    'inventory_entities',
    'inventory_exports',
    'inventory_files',
    'manifest_comparisons',
    'manifest_entries',
    'testing_matrix_entries'
  ]);

  const inventory = runJson('python3', [
    AUDIT_SCRIPT,
    'inventory',
    '--repo-root',
    repoRoot,
    '--side',
    'target',
    '--target-ref',
    'HEAD'
  ]);

  assert.equal(inventory.status, 'completed');
  assert.equal(inventory.counts.r_files, 2);

  const auditDb = path.join(repoRoot, '.refactor-fidelity-audit', 'audit.db');
  const runCount = querySqlite(auditDb, 'select count(*) from audit_runs')[0][0];
  const entityCount = querySqlite(auditDb, 'select count(*) from inventory_entities')[0][0];
  assert.equal(runCount, 1);
  assert.equal(entityCount, 5);

  const events = fs
    .readFileSync(path.join(repoRoot, '.refactor-fidelity-audit', 'events.jsonl'), 'utf8')
    .trim()
    .split('\n')
    .map((line) => JSON.parse(line));
  assert.equal(events.length, 3);
  assert.equal(events[0].phase, 'bootstrap');
  assert.equal(events[1].phase, 'inventory');
  assert.equal(events[2].status, 'completed');

  const summaryMarkdown = fs.readFileSync(
    path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-summary.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /R files: `2`/);
  assert.match(summaryMarkdown, /Exported methods: `1`/);
});

serialTest('surface command persists file and export surfaces and reports drift classes', () => {
  const { baselineRoot, targetRoot } = makeSurfaceDiffPair();
  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'surface',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot
  ]);

  assert.equal(summary.mode, 'surface');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.diff_counts.missing_definition, 1);
  assert.equal(summary.diff_counts.extra_definition, 1);
  assert.equal(summary.diff_counts.duplicate_definition, 1);
  assert.equal(summary.diff_counts.moved_definition, 3);
  assert.equal(summary.diff_counts.missing_export, 1);
  assert.equal(summary.diff_counts.extra_export, 1);
  assert.equal(summary.diff_counts.missing_file, 1);
  assert.equal(summary.diff_counts.extra_file, 2);
  assert.equal(summary.diff_counts.collate_drift, 1);
  assert.equal(summary.exception_count, 9);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const fileSurfaceCount = querySqlite(auditDb, 'select count(*) from inventory_files')[0][0];
  const exportSurfaceCount = querySqlite(auditDb, 'select count(*) from inventory_exports')[0][0];
  const diffCount = querySqlite(auditDb, 'select count(*) from inventory_diffs')[0][0];
  const exceptionCount = querySqlite(auditDb, 'select count(*) from exceptions')[0][0];
  assert.equal(fileSurfaceCount, 5);
  assert.equal(exportSurfaceCount, 6);
  assert.equal(diffCount, 12);
  assert.equal(exceptionCount, 9);

  const diffClasses = querySqlite(
    auditDb,
    'select diff_class, entity_key from inventory_diffs order by diff_class, entity_key'
  );
  assert.deepEqual(diffClasses, [
    ['collate_drift', 'R/b.R'],
    ['duplicate_definition', 'symbol::bar'],
    ['extra_definition', 'symbol::bar'],
    ['extra_export', 'functions::bar'],
    ['extra_file', 'R/c.R'],
    ['extra_file', 'R/d.R'],
    ['missing_definition', 'symbol::helper'],
    ['missing_export', 'methods::processThing'],
    ['missing_file', 'R/a.R'],
    ['moved_definition', 'setGeneric::processThing'],
    ['moved_definition', 'setMethod::processThing::DemoData'],
    ['moved_definition', 'symbol::foo']
  ]);

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-summary.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Total diffs: `12`/);
  assert.match(summaryMarkdown, /Exception candidates: `9`/);
  assert.match(summaryMarkdown, /`missing_definition`: `1`/);
  assert.match(summaryMarkdown, /`moved_definition`: `3`/);
});

serialTest('manifest command catalogs entries, classifies fidelity tiers, and emits exception candidates', () => {
  const { baselineRoot, targetRoot } = makeManifestComparisonPair();
  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'manifest',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--manifest-path',
    'tools/refactor/manifest-fidelity.yml'
  ]);

  assert.equal(summary.mode, 'manifest');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.manifest_count, 1);
  assert.equal(summary.entry_count, 7);
  assert.equal(summary.exception_count, 5);
  assert.equal(summary.status_counts.raw_exact, 3);
  assert.equal(summary.status_counts.normalized_only, 1);
  assert.equal(summary.status_counts.ast_only, 1);
  assert.equal(summary.status_counts.content_missing, 1);
  assert.equal(summary.status_counts.manual_merge_expected, 1);
  assert.equal(summary.resolver_counts.selector, 3);
  assert.equal(summary.resolver_counts.raw_text_search, 1);
  assert.equal(summary.resolver_counts.file_normalized_text_search, 1);
  assert.equal(summary.resolver_counts.content_missing, 1);
  assert.equal(summary.resolver_counts.missing_target_file, 1);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const manifestEntryCount = querySqlite(auditDb, 'select count(*) from manifest_entries')[0][0];
  const comparisonCount = querySqlite(auditDb, 'select count(*) from manifest_comparisons')[0][0];
  const exceptionCount = querySqlite(auditDb, 'select count(*) from exceptions')[0][0];
  assert.equal(manifestEntryCount, 7);
  assert.equal(comparisonCount, 7);
  assert.equal(exceptionCount, 5);

  const comparisons = querySqlite(
    auditDb,
    'select entry_id, status, target_resolver from manifest_entries join manifest_comparisons using (manifest_entry_id) order by entry_id'
  );
  assert.deepEqual(comparisons, [
    ['anchor_normalized', 'normalized_only', 'file_normalized_text_search'],
    ['expr_exact', 'raw_exact', 'raw_text_search'],
    ['manual_expected', 'manual_merge_expected', 'missing_target_file'],
    ['method_exact', 'raw_exact', 'selector'],
    ['missing_target', 'content_missing', 'content_missing'],
    ['symbol_comment_only', 'ast_only', 'selector'],
    ['symbol_exact', 'raw_exact', 'selector']
  ]);

  const exceptions = querySqlite(
    auditDb,
    'select exception_type, severity from exceptions order by exception_type'
  );
  assert.deepEqual(exceptions, [
    ['comment_only_drift', 'low'],
    ['manual_merge_expected', 'medium'],
    ['missing_target_block', 'high'],
    ['target_resolution_asymmetry', 'low'],
    ['whitespace_only_drift', 'low']
  ]);

  const exceptionsJson = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json'),
    'utf8'
  ));
  assert.equal(exceptionsJson.exceptions.length, 5);
  assert.equal(exceptionsJson.summary.open_count, 5);

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-summary.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Manifest entries: `7`/);
  assert.match(summaryMarkdown, /`normalized_only`: `1`/);
  assert.match(summaryMarkdown, /`manual_merge_expected`: `1`/);
});

serialTest('manifest command records source lineage gaps instead of aborting when baseline selectors no longer resolve', () => {
  const { baselineRoot, targetRoot } = makeManifestSourceGapPair();
  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'manifest',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--manifest-path',
    'tools/refactor/manifest-source-gap.yml'
  ]);

  assert.equal(summary.mode, 'manifest');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.entry_count, 1);
  assert.equal(summary.exception_count, 1);
  assert.equal(summary.status_counts.source_missing, 1);
  assert.equal(summary.resolver_counts.source_missing, 1);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const comparisonRows = querySqlite(
    auditDb,
    'select entry_id, status, target_resolver from manifest_entries join manifest_comparisons using (manifest_entry_id)'
  );
  assert.deepEqual(comparisonRows, [[
    'source_gap',
    'source_missing',
    'source_missing'
  ]]);

  const exceptionRows = querySqlite(
    auditDb,
    'select exception_type, severity, reason from exceptions order by exception_type'
  );
  assert.deepEqual(exceptionRows, [[
    'source_lineage_gap',
    'high',
    'Baseline/source block could not be resolved via `selector:symbol`: Selector did not match any block for entry source_gap'
  ]]);

  const latestExceptions = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json'),
    'utf8'
  ));
  assert.equal(latestExceptions.summary.open_count, 1);
  assert.equal(latestExceptions.exceptions[0].exception_type, 'source_lineage_gap');
});

serialTest('behavior command catalogs replay cases, records parity tiers, and emits behavior exceptions', () => {
  const { baselineRoot, targetRoot } = makeBehaviorComparisonPair();
  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'behavior',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--case-path',
    'tools/refactor/behavior-cases-fixture.json'
  ]);

  assert.equal(summary.mode, 'behavior');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.case_file_count, 1);
  assert.equal(summary.case_count, 6);
  assert.equal(summary.exception_count, 2);
  assert.equal(summary.status_counts.exact_match, 2);
  assert.equal(summary.status_counts.normalized_match, 1);
  assert.equal(summary.status_counts.error_match, 1);
  assert.equal(summary.status_counts.mismatch, 2);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const behaviorCaseCount = querySqlite(auditDb, 'select count(*) from behavior_cases')[0][0];
  const behaviorResultCount = querySqlite(auditDb, 'select count(*) from behavior_results')[0][0];
  const exceptionCount = querySqlite(auditDb, 'select count(*) from exceptions')[0][0];
  assert.equal(behaviorCaseCount, 6);
  assert.equal(behaviorResultCount, 6);
  assert.equal(exceptionCount, 2);

  const behaviorRows = querySqlite(
    auditDb,
    'select behavior_case_id, baseline_status, target_status, normalized_equal, warning_equal, message_equal, error_equal, artifact_equal from behavior_results order by behavior_case_id'
  );
  assert.deepEqual(behaviorRows, [
    ['error_match', 'error', 'error', 0, 1, 1, 1, 0],
    ['exact_scalar', 'ok', 'ok', 1, 1, 1, 1, 1],
    ['message_diff', 'ok', 'ok', 1, 1, 0, 1, 1],
    ['normalized_table', 'ok', 'ok', 1, 1, 1, 1, 0],
    ['s4_grid', 'ok', 'ok', 1, 1, 1, 1, 1],
    ['warning_diff', 'ok', 'ok', 1, 0, 1, 1, 1]
  ]);

  const exceptionRows = querySqlite(
    auditDb,
    'select exception_type, severity, entity_key from exceptions order by exception_type'
  );
  assert.deepEqual(exceptionRows, [
    ['behavior_message_mismatch', 'low', 'symbol::messageValue'],
    ['behavior_warning_mismatch', 'medium', 'symbol::warningValue']
  ]);

  const exceptionsJson = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json'),
    'utf8'
  ));
  assert.equal(exceptionsJson.exceptions.length, 2);

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-summary.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Cases: `6`/);
  assert.match(summaryMarkdown, /`normalized_match`: `1`/);
  assert.match(summaryMarkdown, /`mismatch`: `2`/);
});

serialTest('contracts command materializes the testing matrix, catalogs testthat families, and records execution results', () => {
  const targetRoot = makeContractCatalogRepo();
  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'contracts',
    '--repo-root',
    targetRoot,
    '--target-path',
    targetRoot,
    '--execute'
  ]);

  assert.equal(summary.mode, 'contracts');
  assert.equal(summary.status, 'completed_with_failures');
  assert.equal(summary.test_file_count, 5);
  assert.equal(summary.test_case_count, 7);
  assert.equal(summary.executed_file_count, 5);
  assert.equal(summary.exception_count, 1);
  assert.equal(summary.family_counts.module_contract, 1);
  assert.equal(summary.family_counts.characterization, 1);
  assert.equal(summary.family_counts.compat, 1);
  assert.equal(summary.family_counts.golden_master, 1);
  assert.equal(summary.family_counts.general_unit, 1);
  assert.equal(summary.execution_status_counts.failed, 1);
  assert.equal(summary.execution_status_counts.passed, 4);
  assert.ok(summary.matrix_surface_count >= 7);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const matrixCount = querySqlite(auditDb, 'select count(*) from testing_matrix_entries')[0][0];
  const fileCount = querySqlite(auditDb, 'select count(*) from contract_test_files')[0][0];
  const caseCount = querySqlite(auditDb, 'select count(*) from contract_test_cases')[0][0];
  const resultCount = querySqlite(auditDb, 'select count(*) from contract_results')[0][0];
  const exceptionCount = querySqlite(auditDb, 'select count(*) from exceptions')[0][0];
  assert.ok(matrixCount >= 7);
  assert.equal(fileCount, 5);
  assert.equal(caseCount, 7);
  assert.equal(resultCount, 5);
  assert.equal(exceptionCount, 1);

  const fileRows = querySqlite(
    auditDb,
    'select file_path, family, module_family, surface, certainty_target, case_count from contract_test_files order by file_path'
  );
  assert.deepEqual(fileRows, [
    ['tests/testthat/test-prot-01c-import-module-contracts.R', 'module_contract', 'import', 'demonolithed_wrappers_and_modules', 'T5_contract_replayed', 3],
    ['tests/testthat/test-prot-03-rollup.R', 'general_unit', 'rollup', 'unspecified_general', 'T4_behavior_replayed', 1],
    ['tests/testthat/test-prot-07b-da-handlers-compat.R', 'compat', 'da_handlers', 'cross_module_and_workflow_segments', 'T5_contract_replayed', 1],
    ['tests/testthat/test-prot-07c-da-results-characterization.R', 'characterization', 'da_results', 'cross_module_and_workflow_segments', 'T5_contract_replayed', 1],
    ['tests/testthat/test-protein-da-golden-master.R', 'golden_master', 'protein_da', 'cross_module_and_workflow_segments', 'T6_end_to_end_verified', 1]
  ]);

  const scenarioRows = querySqlite(
    auditDb,
    "select test_name, primary_scenario from contract_test_cases where contract_test_file_id = (select contract_test_file_id from contract_test_files where file_path = 'tests/testthat/test-prot-01c-import-module-contracts.R') order by test_name"
  );
  assert.deepEqual(scenarioRows, [
    ['mod_prot_import_server exports reports through download seam', 'export_report_side_effects'],
    ['mod_prot_import_server initializes default outputs', 'initialization'],
    ['mod_prot_import_server reports invalid input errors', 'invalid_input']
  ]);

  const resultRows = querySqlite(
    auditDb,
    'select status, test_count, failure_count, skip_count from contract_results order by contract_result_id'
  );
  assert.deepEqual(resultRows, [
    ['passed', 3, 0, 0],
    ['passed', 1, 0, 0],
    ['failed', 1, 1, 0],
    ['passed', 1, 0, 0],
    ['passed', 1, 0, 1]
  ]);

  const exceptionsJson = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json'),
    'utf8'
  ));
  assert.equal(exceptionsJson.exceptions.length, 1);
  assert.equal(exceptionsJson.exceptions[0].exception_type, 'contract_test_failure');

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-summary.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Test files cataloged: `5`/);
  assert.match(summaryMarkdown, /Executed files: `5`/);
  assert.match(summaryMarkdown, /Exceptions emitted: `1`/);
  assert.match(summaryMarkdown, /`module_contract`: `1`/);
});

serialTest('exceptions and curate commands materialize unresolved review state and support auto and manual curation', () => {
  const { baselineRoot, targetRoot } = makeManifestComparisonPair();
  runJson('python3', [
    AUDIT_SCRIPT,
    'manifest',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--manifest-path',
    'tools/refactor/manifest-fidelity.yml'
  ]);

  const report = runJson('python3', [
    AUDIT_SCRIPT,
    'exceptions',
    '--repo-root',
    targetRoot,
    '--auto-curate'
  ]);

  assert.equal(report.summary.total_count, 5);
  assert.equal(report.summary.auto_curated_count, 3);
  assert.equal(report.summary.open_count, 2);
  assert.equal(report.summary.high_open_count, 1);
  assert.equal(report.summary.status_counts.auto_curated, 3);
  assert.equal(report.summary.status_counts.candidate, 2);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const statusRows = querySqlite(
    auditDb,
    'select exception_type, curation_status, disposition from exceptions order by exception_type'
  );
  assert.deepEqual(statusRows, [
    ['comment_only_drift', 'auto_curated', 'equivalent_comment_only'],
    ['manual_merge_expected', 'candidate', null],
    ['missing_target_block', 'candidate', null],
    ['target_resolution_asymmetry', 'auto_curated', 'equivalent_target_resolution_asymmetry'],
    ['whitespace_only_drift', 'auto_curated', 'equivalent_whitespace_only']
  ]);

  const manualMergeKey = querySqlite(
    auditDb,
    "select exception_key from exceptions where exception_type = 'manual_merge_expected'"
  )[0][0];
  const curateResult = runJson('python3', [
    AUDIT_SCRIPT,
    'curate',
    '--repo-root',
    targetRoot,
    '--exception-key',
    manualMergeKey,
    '--curation-status',
    'accepted',
    '--disposition',
    'expected_manual_merge',
    '--owner',
    'audit-bot',
    '--note',
    'Reviewed manual merge expectation'
  ]);

  assert.equal(curateResult.status, 'curated');
  assert.equal(curateResult.updated_count, 1);
  assert.equal(curateResult.open_count, 1);
  assert.equal(curateResult.high_open_count, 1);

  const openReport = runJson('python3', [
    AUDIT_SCRIPT,
    'exceptions',
    '--repo-root',
    targetRoot,
    '--open-only'
  ]);
  assert.equal(openReport.summary.open_count, 1);
  assert.equal(openReport.exceptions.length, 1);
  assert.equal(openReport.exceptions[0].exception_type, 'missing_target_block');

  const latestExceptions = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json'),
    'utf8'
  ));
  assert.equal(latestExceptions.summary.open_count, 1);
  assert.equal(latestExceptions.summary.high_open_count, 1);

  const latestExceptionsMd = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.md'),
    'utf8'
  );
  assert.match(latestExceptionsMd, /High-Severity Open Exceptions/);
  assert.match(latestExceptionsMd, /missing_target_block/);
});

serialTest('closeout command orchestrates integrated audits and emits a readiness summary', () => {
  const { baselineRoot, mainlineRoot, targetRoot } = makeCloseoutComparisonTriplet();
  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'closeout',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--mainline-path',
    mainlineRoot,
    '--target-path',
    targetRoot,
    '--manifest-path',
    'tools/refactor/manifest-closeout.yml',
    '--case-path',
    'tools/refactor/behavior-cases-closeout.json',
    '--test-file',
    'tests/testthat/test-prot-03-rollup.R',
    '--contracts-execute'
  ]);

  assert.equal(summary.mode, 'closeout');
  assert.equal(summary.status, 'ready_for_coverage_campaign');
  assert.equal(summary.proven_parity, true);
  assert.equal(summary.auto_curated_count, 1);
  assert.equal(summary.readiness.coverage_campaign_ready, true);
  assert.equal(summary.readiness.open_exception_count, 0);
  assert.equal(summary.readiness.high_open_exception_count, 0);
  assert.equal(summary.readiness.pinned_baseline_explicit, true);
  assert.equal(summary.readiness.contract_gate, true);
  assert.equal(summary.components.baseline_surface.total_diffs, 0);
  assert.equal(summary.components.mainline_surface.total_diffs, 0);
  assert.equal(summary.components.manifest.exception_count, 1);
  assert.equal(summary.components.behavior.exception_count, 0);
  assert.equal(summary.components.contracts.exception_count, 0);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const runModes = querySqlite(
    auditDb,
    'select mode, count(*) from audit_runs group by mode order by mode'
  );
  assert.deepEqual(runModes, [
    ['behavior', 1],
    ['closeout', 1],
    ['contracts', 1],
    ['manifest', 1],
    ['surface', 2]
  ]);

  const latestExceptions = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json'),
    'utf8'
  ));
  assert.equal(latestExceptions.summary.open_count, 0);
  assert.equal(latestExceptions.exceptions.length, 0);

  const latestCloseout = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-closeout.json'),
    'utf8'
  ));
  assert.equal(latestCloseout.status, 'ready_for_coverage_campaign');

  const closeoutMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-closeout.md'),
    'utf8'
  );
  assert.match(closeoutMarkdown, /Coverage-campaign ready: `true`/);
  assert.match(closeoutMarkdown, /Proven parity: `true`/);
});

serialTest('shared manifest selector library drives verify_refactor and extract_blocks across selector kinds', () => {
  const repoRoot = makeManifestFixtureRepo();
  const outputRoot = makeSandbox('fidelity-extract-');

  const verifyOutput = runText('Rscript', [
    path.resolve(process.cwd(), 'tools/refactor/verify_refactor.R'),
    '--manifest',
    'tools/refactor/manifest-test.yml',
    '--repo-root',
    repoRoot
  ]);
  assert.match(verifyOutput, /Verification passed\./);

  const extractOutput = runText('Rscript', [
    path.resolve(process.cwd(), 'tools/refactor/extract_blocks.R'),
    '--manifest',
    'tools/refactor/manifest-test.yml',
    '--repo-root',
    repoRoot,
    '--output-root',
    outputRoot,
    '--write-targets',
    '--emit-collate',
    'collate-test.txt'
  ]);
  assert.match(extractOutput, /Resolved entries:/);

  const extractedSymbols = fs.readFileSync(path.join(outputRoot, 'R', 'extracted_symbols.R'), 'utf8');
  const extractedS4 = fs.readFileSync(path.join(outputRoot, 'R', 'extracted_s4.R'), 'utf8');
  const extractedAnchor = fs.readFileSync(path.join(outputRoot, 'R', 'extracted_anchor.R'), 'utf8');
  const emittedCollate = fs.readFileSync(path.join(outputRoot, 'collate-test.txt'), 'utf8');

  assert.match(extractedSymbols, /foo <- function/);
  assert.match(extractedSymbols, /helper = function/);
  assert.match(extractedS4, /setMethod\("processThing", "DemoData"/);
  assert.match(extractedAnchor, /# BEGIN BLOCK/);
  assert.match(extractedAnchor, /anchor_value <- 1/);
  assert.match(emittedCollate, /R\/extracted_symbols\.R/);
  assert.match(emittedCollate, /R\/extracted_anchor\.R/);
});

for (const { fn } of queuedTests) {
  fn();
}
