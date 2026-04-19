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
  writeFile(path.join(root, 'R', 'b.R'), `DemoData <- setClass("DemoData", representation(value = "numeric"))

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

function makeSurfaceParityPair() {
  const baselineRoot = makeSandbox('fidelity-surface-parity-baseline-');
  const targetRoot = makeSandbox('fidelity-surface-parity-target-');

  const description = `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Surface parity fixture.
License: MIT
Encoding: UTF-8
Collate:
    'a.R'
    'b.R'
`;
  const namespace = `export(foo)
`;
  const sharedA = `foo <- function(x) {
  x + 1
}
`;
  const sharedB = `dup <- function(y) {
  y
}

dup <- function(y) {
  y
}
`;

  for (const root of [baselineRoot, targetRoot]) {
    writeFile(path.join(root, 'DESCRIPTION'), description);
    writeFile(path.join(root, 'NAMESPACE'), namespace);
    writeFile(path.join(root, 'R', 'a.R'), sharedA);
    writeFile(path.join(root, 'R', 'b.R'), sharedB);
  }

  return { baselineRoot, targetRoot };
}

function makeSurfaceDedupSafePair() {
  const baselineRoot = makeSandbox('fidelity-surface-dedup-baseline-');
  const targetRoot = makeSandbox('fidelity-surface-dedup-target-');

  const description = `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Surface duplicate consolidation fixture.
License: MIT
Encoding: UTF-8
Collate:
    'a.R'
    'b.R'
    'c.R'
`;
  const namespace = `export(foo)
`;

  for (const root of [baselineRoot, targetRoot]) {
    writeFile(path.join(root, 'DESCRIPTION'), description);
    writeFile(path.join(root, 'NAMESPACE'), namespace);
    writeFile(path.join(root, 'R', 'a.R'), `foo <- function(x) {\n  x + 1\n}\n`);
  }

  writeFile(path.join(baselineRoot, 'R', 'b.R'), `dup <- function(y) {\n  y\n}\n`);
  writeFile(path.join(baselineRoot, 'R', 'c.R'), `dup <- function(y) {\n  y\n}\n\nkeeper <- function(z) {\n  z\n}\n`);

  writeFile(path.join(targetRoot, 'R', 'b.R'), `dup <- function(y) {\n  y\n}\n`);
  writeFile(path.join(targetRoot, 'R', 'c.R'), `keeper <- function(z) {\n  z\n}\n`);

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

function makeManifestRepoSelectorTargetPair() {
  const baselineRoot = makeSandbox('fidelity-manifest-repo-target-baseline-');
  const targetRoot = makeSandbox('fidelity-manifest-repo-target-target-');

  writeFile(path.join(baselineRoot, 'R', 'source_symbol.R'), `foo <- function(x) {
  x + 1
}
`);
  writeFile(path.join(targetRoot, 'R', 'target_symbol.R'), `# breadcrumb placeholder for later split
`);
  writeFile(path.join(targetRoot, 'R', 'target_symbol_final.R'), `foo <- function(x) {
  x + 1
}
`);
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-repo-selector.yml'), `version: 1
entries:
  - id: symbol_repo_selector
    source: R/source_symbol.R
    selector: { kind: symbol, value: foo }
    target: R/target_symbol.R
    group: fidelity.repo_target
`);

  return { baselineRoot, targetRoot };
}

function makeManifestTargetSelectorPair() {
  const baselineRoot = makeSandbox('fidelity-manifest-target-selector-baseline-');
  const targetRoot = makeSandbox('fidelity-manifest-target-selector-target-');

  writeFile(path.join(baselineRoot, 'R', 'source_anchor.R'), `# LEGACY START
legacyDispatch <- function(x) {
  x + 1
}

legacyMethod <- function(y) {
  y + 2
}
# LEGACY END
`);
  writeFile(path.join(targetRoot, 'R', 'target_anchor.R'), `# LEGACY START
legacyDispatch <- function(x) {
  x + 1
}

legacyMethod <- function(y) {
  y + 2
}
# NEXT SECTION
`);
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-target-selector.yml'), `version: 1
entries:
  - id: anchor_target_selector
    source: R/source_anchor.R
    selector:
      kind: anchor_range
      start: "^# LEGACY START$"
      end: "^# LEGACY END$"
      end_inclusive: false
    target_selector:
      kind: anchor_range
      start: "^# LEGACY START$"
      end: "^# NEXT SECTION$"
      end_inclusive: false
    target: R/target_anchor.R
    group: fidelity.target_selector
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

function makeCloseoutOverlapTriplet() {
  const baselineRoot = makeSandbox('fidelity-closeout-overlap-baseline-');
  const mainlineRoot = makeSandbox('fidelity-closeout-overlap-mainline-');
  const targetRoot = makeSandbox('fidelity-closeout-overlap-target-');

  const description = `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Closeout overlap fixture package.
License: MIT
Encoding: UTF-8
Collate:
    'core.R'
`;
  const namespace = `export(overlapHelper)
export(exactValue)
`;
  const baselineCore = `overlapHelper <- function(x) {
  x + 1L
}

exactValue <- function(x) {
  x + 1L
}
`;
  const targetCore = `overlapHelper <- function(x) {
  x + 2L
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

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-closeout-overlap.yml'), `version: 1
entries:
  - id: overlap_helper
    source: R/core.R
    selector: { kind: symbol, value: overlapHelper }
    target: R/core.R
    group: fidelity.overlap
`);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'behavior-cases-closeout-overlap.json'), JSON.stringify({
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

function makeCloseoutCatalogCuratedTriplet() {
  const baselineRoot = makeSandbox('fidelity-closeout-curation-baseline-');
  const mainlineRoot = makeSandbox('fidelity-closeout-curation-mainline-');
  const targetRoot = makeSandbox('fidelity-closeout-curation-target-');

  const description = `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Closeout catalog-curation fixture package.
License: MIT
Encoding: UTF-8
Collate:
    'a.R'
    'b.R'
    'c.R'
`;
  const namespace = `export(exactValue)
`;
  const baselineA = `exactValue <- function(x) {
  x + 1L
}
`;
  const baselineB = `dupHelper <- function(x) {
  x
}
`;
  const baselineC = `dupHelper <- function(x) {
  x + 1L
}

keeper <- function(y) {
  y
}
`;
  const targetA = baselineA;
  const targetB = `dupHelper <- function(x) {
  x
}
`;
  const targetC = `keeper <- function(y) {
  y
}
`;

  for (const root of [baselineRoot, mainlineRoot, targetRoot]) {
    writeFile(path.join(root, 'DESCRIPTION'), description);
    writeFile(path.join(root, 'NAMESPACE'), namespace);
  }

  writeFile(path.join(baselineRoot, 'R', 'a.R'), baselineA);
  writeFile(path.join(baselineRoot, 'R', 'b.R'), baselineB);
  writeFile(path.join(baselineRoot, 'R', 'c.R'), baselineC);

  for (const root of [mainlineRoot, targetRoot]) {
    writeFile(path.join(root, 'R', 'a.R'), targetA);
    writeFile(path.join(root, 'R', 'b.R'), targetB);
    writeFile(path.join(root, 'R', 'c.R'), targetC);
  }

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-closeout-curation.yml'), `version: 1
entries:
  - id: exact_value
    source: R/a.R
    selector: { kind: symbol, value: exactValue }
    target: R/a.R
    group: fidelity.exact
`);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'behavior-cases-closeout-curation.json'), JSON.stringify({
    cases: [
      {
        id: 'exact_value',
        family: 'pure_helper',
        entity_key: 'symbol::exactValue',
        fixture_kind: 'inline_setup',
        normalizer_key: 'scalar_integer',
        risk_level: 'low',
        files: ['R/a.R'],
        setup_code: ['x <- 2L'],
        call_expr: 'exactValue(x)',
        normalize_expr: 'as.integer(result)'
      }
    ]
  }, null, 2));

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'fidelity-curations.json'), JSON.stringify({
    curations: [
      {
        id: 'dup_helper_canonicalized',
        exception_key: 'surface_inventory::surface_duplicate_definition::symbol::dupHelper',
        curation_status: 'accepted',
        disposition: 'canonical_duplicate_retired',
        owner: 'fixture',
        review_note: 'Fixture curation accepts the canonical duplicate reduction.'
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

function makeCloseoutCatalogOverlapTriplet() {
  const baselineRoot = makeSandbox('fidelity-closeout-catalog-overlap-baseline-');
  const mainlineRoot = makeSandbox('fidelity-closeout-catalog-overlap-mainline-');
  const targetRoot = makeSandbox('fidelity-closeout-catalog-overlap-target-');

  const description = `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Closeout catalog-overlap fixture package.
License: MIT
Encoding: UTF-8
Collate:
    'core.R'
`;
  const namespace = `export(overlapHelper)
export(exactValue)
`;
  const baselineCore = `overlapHelper <- function(x) {
  x + 1L
}

exactValue <- function(x) {
  x + 1L
}
`;
  const targetCore = `overlapHelper <- function(x) {
  x + 2L
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

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-closeout-catalog-overlap.yml'), `version: 1
entries:
  - id: overlap_helper
    source: R/core.R
    selector: { kind: symbol, value: overlapHelper }
    target: R/core.R
    group: fidelity.overlap
`);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'behavior-cases-closeout-catalog-overlap.json'), JSON.stringify({
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

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'fidelity-curations.json'), JSON.stringify({
    curations: [
      {
        id: 'overlap_helper_manifest_curated',
        exception_key: 'manifest_fidelity::semantic_mismatch::manifest::tools/refactor/manifest-closeout-catalog-overlap.yml::overlap_helper',
        curation_status: 'accepted',
        disposition: 'entrypoint_shell_equivalent_under_tests',
        owner: 'fixture',
        review_note: 'Fixture curation accepts the manifest semantic mismatch so the paired surface drift should collapse as redundant.'
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

function makeCoverageRepo() {
  const root = makeSandbox('fidelity-coverage-');
  writeFile(path.join(root, 'DESCRIPTION'), `Package: CoverageDemo
Title: Coverage Demo
Version: 0.0.1
Description: Demo package for fidelity coverage tests.
License: MIT
Encoding: UTF-8
Suggests:
    testthat
`);
  writeFile(path.join(root, 'NAMESPACE'), `export(foo)
`);
  writeFile(path.join(root, 'R', 'a.R'), `foo <- function(x) {
  if (x > 0) x + 1 else x - 1
}
`);
  writeFile(path.join(root, 'tests', 'testthat', 'test-coverage-demo.R'), `library(testthat)

test_that("foo adjusts positive values", {
  expect_equal(foo(1), 2)
})
`);
  writeFile(path.join(root, 'tools', 'refactor', 'coverage-bundles-fixture.json'), JSON.stringify({
    bundles: [
      {
        bundle_id: 'foo_bundle',
        bundle_type: 'helper',
        source_file: 'R/a.R',
        shared_test_files: ['tests/testthat/test-coverage-demo.R'],
        target_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/a.R',
            line_start: 1,
            line_end: 3
          }
        ]
      }
    ]
  }, null, 2));
  return root;
}

function makeCoverageCompareTriplet() {
  const baselineRoot = makeSandbox('fidelity-coverage-compare-baseline-');
  const targetRoot = makeSandbox('fidelity-coverage-compare-target-');
  const description = `Package: CoverageCompareDemo
Title: Coverage Compare Demo
Version: 0.0.1
Description: Demo package for comparative fidelity coverage tests.
License: MIT
Encoding: UTF-8
Suggests:
    testthat
`;
  const namespace = `export(foo)
`;

  for (const root of [baselineRoot, targetRoot]) {
    writeFile(path.join(root, 'DESCRIPTION'), description);
    writeFile(path.join(root, 'NAMESPACE'), namespace);
  }

  writeFile(path.join(baselineRoot, 'R', 'a.R'), `foo <- function(x) {
  if (x > 0) x + 1 else x - 1
}
`);
  writeFile(path.join(targetRoot, 'R', 'b.R'), `foo <- function(x) {
  if (x > 0) x + 1 else x - 1
}
`);
  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-coverage-demo.R'), `library(testthat)

test_that("foo adjusts positive values", {
  expect_equal(foo(1), 2)
})
`);
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'coverage-compare-bundles.json'), JSON.stringify({
    bundles: [
      {
        bundle_id: 'foo_compare_bundle',
        bundle_type: 'helper_surface',
        shared_test_files: ['tests/testthat/test-coverage-demo.R'],
        baseline_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/a.R',
            line_start: 1,
            line_end: 3
          }
        ],
        target_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/b.R',
            line_start: 1,
            line_end: 3
          }
        ]
      }
    ]
  }, null, 2));
  return { baselineRoot, targetRoot };
}

function makeBundleMapTriplet() {
  const baselineRoot = makeSandbox('fidelity-bundle-map-baseline-');
  const mainlineRoot = makeSandbox('fidelity-bundle-map-mainline-');
  const targetRoot = makeSandbox('fidelity-bundle-map-target-');

  const description = `Package: DemoPkg
Title: Demo
Version: 0.0.1
Description: Bundle-map fixture package.
License: MIT
Encoding: UTF-8
Collate:
    'core.R'
    'source_gap.R'
    'source_gap_target.R'
`;
  const namespace = `export(overlapHelper)
export(exactValue)
`;
  for (const root of [baselineRoot, mainlineRoot, targetRoot]) {
    writeFile(path.join(root, 'DESCRIPTION'), description);
    writeFile(path.join(root, 'NAMESPACE'), namespace);
  }

  writeFile(path.join(baselineRoot, 'R', 'core.R'), `overlapHelper <- function(x) {
  x + 1L
}

exactValue <- function(x) {
  x + 1L
}
`);
  writeFile(path.join(mainlineRoot, 'R', 'core.R'), `overlapHelper <- function(x) {
  x + 2L
}

exactValue <- function(x) {
  x + 1L
}
`);
  writeFile(path.join(targetRoot, 'R', 'core.R'), `overlapHelper <- function(x) {
  x + 2L
}

exactValue <- function(x) {
  x + 1L
}
`);

  writeFile(path.join(baselineRoot, 'R', 'source_gap.R'), `otherHelper <- function(x) {
  x
}
`);
  writeFile(path.join(mainlineRoot, 'R', 'source_gap.R'), `otherHelper <- function(x) {
  x
}
`);
  writeFile(path.join(targetRoot, 'R', 'source_gap.R'), `otherHelper <- function(x) {
  x
}
`);
  writeFile(path.join(baselineRoot, 'R', 'source_gap_target.R'), `# target placeholder before extraction
`);
  writeFile(path.join(mainlineRoot, 'R', 'source_gap_target.R'), `# target placeholder before extraction
`);
  writeFile(path.join(targetRoot, 'R', 'source_gap_target.R'), `gapHelper <- function(x) {
  x
}

gapObserverRegister <- function(x) {
  x + 1L
}
`);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'manifest-closeout-bundle-map.yml'), `version: 1
entries:
  - id: overlap_helper
    source: R/core.R
    selector: { kind: symbol, value: overlapHelper }
    target: R/core.R
    group: fidelity.overlap

  - id: source_gap
    source: R/source_gap.R
    selector: { kind: symbol, value: gapHelper }
    target: R/source_gap_target.R
    group: fidelity.lineage

  - id: source_gap_observer_register
    source: R/source_gap.R
    selector: { kind: symbol, value: gapObserverRegister }
    target: R/source_gap_target.R
    group: fidelity.lineage
`);

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'behavior-cases-closeout-bundle-map.json'), JSON.stringify({
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

  writeFile(path.join(targetRoot, 'tools', 'refactor', 'fidelity-curations.json'), JSON.stringify({
    curations: [
      {
        id: 'bundle_map_overlap_helper_curated',
        exception_key: 'manifest_fidelity::semantic_mismatch::manifest::tools/refactor/manifest-closeout-bundle-map.yml::overlap_helper',
        curation_status: 'accepted',
        disposition: 'helper_equivalent_under_tests',
        owner: 'fixture',
        review_note: 'Fixture curation accepts overlapHelper as helper-equivalent under tests.'
      }
    ]
  }, null, 2));

  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-overlap-helper-characterization.R'), `library(testthat)

test_that("overlapHelper retains curated semantics", {
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
  assert.equal(fs.existsSync(path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage.json')), true);
  assert.equal(fs.existsSync(path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-bundles.json')), true);

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
    'coverage_bundle_members',
    'coverage_bundle_results',
    'coverage_bundles',
    'coverage_exception_links',
    'coverage_files',
    'coverage_runs',
    'coverage_test_links',
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
  assert.equal(entityCount, 6);

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

serialTest('coverage command records a deterministic tool-unavailable result', () => {
  const repoRoot = makeCoverageRepo();
  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage',
    '--repo-root',
    repoRoot,
    '--side',
    'target',
    '--source-path',
    repoRoot,
    '--test-file',
    'tests/testthat/test-coverage-demo.R'
  ], {
    env: {
      ...process.env,
      FIDELITY_COVERAGE_FORCE_UNAVAILABLE: '1'
    }
  });

  assert.equal(summary.mode, 'coverage');
  assert.equal(summary.status, 'tool_unavailable');
  assert.equal(summary.tool_status, 'forced_unavailable');
  assert.equal(summary.test_file_count, 1);
  assert.equal(summary.bundle_count, 0);
  assert.equal(summary.file_count, 0);
  assert.equal(summary.line_percent, null);

  const auditDb = path.join(repoRoot, '.refactor-fidelity-audit', 'audit.db');
  const coverageRunRows = querySqlite(
    auditDb,
    'select side, status, tool_status, test_file_count, bundle_count, file_count from coverage_runs'
  );
  assert.deepEqual(coverageRunRows, [
    ['target', 'tool_unavailable', 'forced_unavailable', 1, 0, 0]
  ]);

  const coverageFileCount = querySqlite(auditDb, 'select count(*) from coverage_files')[0][0];
  const bundleCount = querySqlite(auditDb, 'select count(*) from coverage_bundles')[0][0];
  assert.equal(coverageFileCount, 0);
  assert.equal(bundleCount, 0);

  const summaryMarkdown = fs.readFileSync(
    path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Tool status: `forced_unavailable`/);
  assert.match(summaryMarkdown, /Line coverage: `n\/a`/);
});

serialTest('coverage command persists bundle and file coverage from fixture payload', () => {
  const repoRoot = makeCoverageRepo();
  const fixturePath = path.join(repoRoot, 'coverage-fixture.json');
  writeFile(fixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: repoRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/a.R',
        line_percent: 66.6667,
        lines_total: 3,
        lines_covered: 2,
        covered_lines: [1, 2],
        uncovered_lines: [3]
      }
    ],
    line_percent: 66.6667,
    lines_total: 3,
    lines_covered: 2,
    notes: ['fixture payload']
  }, null, 2));

  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage',
    '--repo-root',
    repoRoot,
    '--side',
    'target',
    '--source-path',
    repoRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-bundles-fixture.json'
  ], {
    env: {
      ...process.env,
      FIDELITY_COVERAGE_FIXTURE_JSON: fixturePath
    }
  });

  assert.equal(summary.mode, 'coverage');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.tool_status, 'ok');
  assert.equal(summary.test_file_count, 1);
  assert.equal(summary.bundle_count, 1);
  assert.equal(summary.file_count, 1);
  assert.equal(summary.lines_total, 3);
  assert.equal(summary.lines_covered, 2);
  assert.equal(summary.bundle_status_counts.completed, 1);
  assert.equal(summary.sample_bundle_results[0].bundle_id, 'foo_bundle');
  assert.equal(summary.sample_bundle_results[0].lines_total, 3);
  assert.equal(summary.sample_bundle_results[0].lines_covered, 2);

  const auditDb = path.join(repoRoot, '.refactor-fidelity-audit', 'audit.db');
  const runModes = querySqlite(auditDb, 'select mode, count(*) from audit_runs group by mode order by mode');
  assert.deepEqual(runModes, [['coverage', 1]]);

  const coverageRunRows = querySqlite(
    auditDb,
    'select status, tool_status, test_file_count, bundle_count, file_count, round(line_percent, 1), lines_total, lines_covered from coverage_runs'
  );
  assert.deepEqual(coverageRunRows, [
    ['completed', 'ok', 1, 1, 1, 66.7, 3, 2]
  ]);

  const coverageFileRows = querySqlite(
    auditDb,
    'select file_path, round(line_percent, 1), lines_total, lines_covered from coverage_files'
  );
  assert.deepEqual(coverageFileRows, [['R/a.R', 66.7, 3, 2]]);

  const bundleRows = querySqlite(
    auditDb,
    'select bundle_id, bundle_type, source_file_path from coverage_bundles'
  );
  assert.deepEqual(bundleRows, [['foo_bundle', 'helper', 'R/a.R']]);

  const memberRows = querySqlite(
    auditDb,
    'select entity_key, file_path, line_start, line_end from coverage_bundle_members'
  );
  assert.deepEqual(memberRows, [['symbol::foo', 'R/a.R', 1, 3]]);

  const bundleResultRows = querySqlite(
    auditDb,
    'select status, round(line_percent, 1), lines_total, lines_covered from coverage_bundle_results'
  );
  assert.deepEqual(bundleResultRows, [['completed', 66.7, 3, 2]]);

  const testLinkRows = querySqlite(
    auditDb,
    'select coverage_bundle_id is not null, test_file_path, source from coverage_test_links order by source, test_file_path'
  );
  assert.deepEqual(testLinkRows, [
    [1, 'tests/testthat/test-coverage-demo.R', 'bundle_manifest']
  ]);

  const latestCoverage = JSON.parse(fs.readFileSync(
    path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage.json'),
    'utf8'
  ));
  assert.equal(latestCoverage.status, 'completed');
  assert.equal(latestCoverage.bundle_count, 1);

  const summaryMarkdown = fs.readFileSync(
    path.join(repoRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Bundle manifest: `tools\/refactor\/coverage-bundles-fixture\.json`/);
  assert.match(summaryMarkdown, /Line coverage: `66\.7%`/);
});

serialTest('coverage-compare command persists baseline and target bundle deltas from shared tests', () => {
  const { baselineRoot, targetRoot } = makeCoverageCompareTriplet();
  const baselineFixturePath = path.join(targetRoot, 'baseline-coverage-fixture.json');
  const targetFixturePath = path.join(targetRoot, 'target-coverage-fixture.json');
  writeFile(baselineFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: baselineRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/a.R',
        line_percent: 33.3333,
        lines_total: 3,
        lines_covered: 1,
        covered_lines: [1],
        uncovered_lines: [2, 3]
      }
    ],
    line_percent: 33.3333,
    lines_total: 3,
    lines_covered: 1,
    notes: ['baseline fixture payload']
  }, null, 2));
  writeFile(targetFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: targetRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/b.R',
        line_percent: 100,
        lines_total: 3,
        lines_covered: 3,
        covered_lines: [1, 2, 3],
        uncovered_lines: []
      }
    ],
    line_percent: 100,
    lines_total: 3,
    lines_covered: 3,
    notes: ['target fixture payload']
  }, null, 2));

  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-compare',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-compare-bundles.json'
  ], {
    env: {
      ...process.env,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_BASELINE: baselineFixturePath,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_TARGET: targetFixturePath
    }
  });

  assert.equal(summary.mode, 'coverage_compare');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.bundle_count, 1);
  assert.equal(summary.test_file_count, 1);
  assert.equal(summary.baseline.line_percent_display, '33.3%');
  assert.equal(summary.target.line_percent_display, '100.0%');
  assert.equal(summary.bundle_status_counts.completed, 1);
  assert.equal(summary.delta_counts.improved, 1);
  assert.equal(summary.bundle_results[0].bundle_id, 'foo_compare_bundle');
  assert.equal(summary.bundle_results[0].baseline_line_percent.toFixed(1), '33.3');
  assert.equal(summary.bundle_results[0].target_line_percent.toFixed(1), '100.0');
  assert.equal(summary.bundle_results[0].line_percent_delta.toFixed(1), '66.7');

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const runModes = querySqlite(auditDb, 'select mode, count(*) from audit_runs group by mode order by mode');
  assert.deepEqual(runModes, [['coverage_compare', 1]]);

  const coverageRunRows = querySqlite(
    auditDb,
    'select side, status, tool_status, bundle_count, file_count from coverage_runs order by side'
  );
  assert.deepEqual(coverageRunRows, [
    ['baseline', 'completed', 'ok', 1, 1],
    ['target', 'completed', 'ok', 1, 1]
  ]);

  const memberRows = querySqlite(
    auditDb,
    `select coverage_bundles.side, coverage_bundle_members.file_path
       from coverage_bundle_members
       join coverage_bundles on coverage_bundles.coverage_bundle_id = coverage_bundle_members.coverage_bundle_id
      where coverage_bundles.run_id = (select run_id from audit_runs where mode = "coverage_compare")
      order by coverage_bundles.side`
  );
  assert.deepEqual(memberRows, [
    ['baseline', 'R/a.R'],
    ['target', 'R/b.R']
  ]);

  const comparisonRow = querySqlite(
    auditDb,
    `select coverage_bundle_results.status, round(coverage_bundle_results.line_percent, 1), coverage_bundle_results.result_json
       from coverage_bundle_results
       join coverage_bundles on coverage_bundles.coverage_bundle_id = coverage_bundle_results.coverage_bundle_id
      where coverage_bundles.run_id = (select run_id from audit_runs where mode = "coverage_compare")
        and coverage_bundles.side = "comparison"`
  )[0];
  assert.equal(comparisonRow[0], 'completed');
  assert.equal(comparisonRow[1], 100.0);
  const comparisonPayload = JSON.parse(comparisonRow[2]);
  assert.equal(comparisonPayload.baseline.line_percent.toFixed(1), '33.3');
  assert.equal(comparisonPayload.target.line_percent.toFixed(1), '100.0');
  assert.equal(comparisonPayload.line_percent_delta.toFixed(1), '66.7');

  const latestCoverageCompare = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage-compare.json'),
    'utf8'
  ));
  assert.equal(latestCoverageCompare.status, 'completed');
  assert.equal(latestCoverageCompare.delta_counts.improved, 1);

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage-compare.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Refactor Fidelity Coverage Compare/);
  assert.match(summaryMarkdown, /Delta Counts/);
});

serialTest('coverage-compare explicit-tests-only narrows the family test subset', () => {
  const { baselineRoot, targetRoot } = makeCoverageCompareTriplet();
  writeFile(path.join(targetRoot, 'tests', 'testthat', 'test-extra-coverage-demo.R'), `library(testthat)

test_that("extra coverage fixture test", {
  expect_true(TRUE)
})
`);
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'coverage-compare-bundles-explicit.json'), JSON.stringify({
    bundles: [
      {
        bundle_id: 'foo_compare_bundle',
        bundle_type: 'helper_surface',
        shared_test_files: [
          'tests/testthat/test-coverage-demo.R',
          'tests/testthat/test-extra-coverage-demo.R'
        ],
        baseline_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/a.R',
            line_start: 1,
            line_end: 3
          }
        ],
        target_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/b.R',
            line_start: 1,
            line_end: 3
          }
        ]
      }
    ]
  }, null, 2));

  const baselineFixturePath = path.join(targetRoot, 'baseline-coverage-explicit-fixture.json');
  const targetFixturePath = path.join(targetRoot, 'target-coverage-explicit-fixture.json');
  writeFile(baselineFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: baselineRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/a.R',
        line_percent: 33.3333,
        lines_total: 3,
        lines_covered: 1,
        covered_lines: [1],
        uncovered_lines: [2, 3]
      }
    ],
    line_percent: 33.3333,
    lines_total: 3,
    lines_covered: 1,
    notes: ['baseline explicit fixture']
  }, null, 2));
  writeFile(targetFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: targetRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/b.R',
        line_percent: 100,
        lines_total: 3,
        lines_covered: 3,
        covered_lines: [1, 2, 3],
        uncovered_lines: []
      }
    ],
    line_percent: 100,
    lines_total: 3,
    lines_covered: 3,
    notes: ['target explicit fixture']
  }, null, 2));

  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-compare',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-compare-bundles-explicit.json',
    '--explicit-tests-only',
    '--test-file',
    'tests/testthat/test-coverage-demo.R'
  ], {
    env: {
      ...process.env,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_BASELINE: baselineFixturePath,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_TARGET: targetFixturePath
    }
  });

  assert.equal(summary.mode, 'coverage_compare');
  assert.equal(summary.test_file_count, 1);
  assert.deepEqual(summary.selected_test_files, ['tests/testthat/test-coverage-demo.R']);
});

serialTest('coverage-evidence treats baseline test failures under target-selected tests as comparison gaps', () => {
  const { baselineRoot, targetRoot } = makeCoverageCompareTriplet();
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'coverage-evidence-bundles-baseline-failures.json'), JSON.stringify({
    bundles: [
      {
        bundle_id: 'foo_compare_bundle',
        bundle_type: 'helper_surface',
        shared_test_files: ['tests/testthat/test-coverage-demo.R'],
        baseline_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/a.R',
            line_start: 1,
            line_end: 20
          }
        ],
        target_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/b.R',
            line_start: 1,
            line_end: 20
          }
        ]
      }
    ]
  }, null, 2));

  const baselineFixturePath = path.join(targetRoot, 'baseline-coverage-evidence-failures-fixture.json');
  const targetFixturePath = path.join(targetRoot, 'target-coverage-evidence-failures-fixture.json');
  writeFile(baselineFixturePath, JSON.stringify({
    status: 'test_failures',
    tool_status: 'testthat_failures',
    project_root: baselineRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    test_count: 1,
    failure_count: 1,
    skip_count: 0,
    failed_tests: ['shared compatibility test'],
    files: [
      {
        file_path: 'R/a.R',
        line_percent: 0,
        lines_total: 20,
        lines_covered: 0,
        covered_lines: [],
        uncovered_lines: Array.from({ length: 20 }, (_, idx) => idx + 1)
      }
    ],
    line_percent: 0,
    lines_total: 20,
    lines_covered: 0,
    notes: ['baseline selected tests fail']
  }, null, 2));
  writeFile(targetFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: targetRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    test_count: 1,
    failure_count: 0,
    skip_count: 0,
    failed_tests: [],
    files: [
      {
        file_path: 'R/b.R',
        line_percent: 95,
        lines_total: 20,
        lines_covered: 19,
        covered_lines: Array.from({ length: 19 }, (_, idx) => idx + 1),
        uncovered_lines: [20]
      }
    ],
    line_percent: 95,
    lines_total: 20,
    lines_covered: 19,
    notes: ['target selected tests pass']
  }, null, 2));

  const compareSummary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-compare',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-evidence-bundles-baseline-failures.json'
  ], {
    env: {
      ...process.env,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_BASELINE: baselineFixturePath,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_TARGET: targetFixturePath
    }
  });

  assert.equal(compareSummary.status, 'completed_with_incomplete_sides');
  assert.equal(compareSummary.baseline.status, 'test_failures');
  assert.equal(compareSummary.baseline.failure_count, 1);
  assert.equal(compareSummary.target.status, 'completed');
  assert.equal(compareSummary.target.failure_count, 0);

  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-evidence',
    '--repo-root',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-evidence-bundles-baseline-failures.json',
    '--coverage-compare-run-id',
    compareSummary.run_id
  ]);

  assert.equal(summary.status, 'completed_with_gaps');
  assert.equal(summary.low_coverage_bundle_count, 0);
  assert.equal(summary.comparison_gap_bundle_count, 1);
  assert.equal(summary.exception_count, 1);
  assert.equal(summary.sample_bundles[0].coverage_gate_status, 'comparison_gap');
  assert.equal(summary.bundles[0].baseline_status, 'test_failures');
  assert.equal(summary.bundles[0].target_status, 'completed');

  const latestCompareMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage-compare.md'),
    'utf8'
  );
  assert.match(latestCompareMarkdown, /Failing tests: `1`/);
  assert.match(latestCompareMarkdown, /shared compatibility test/);
});

serialTest('coverage-evidence command materializes low-coverage bundle evidence and exceptions', () => {
  const { baselineRoot, targetRoot } = makeCoverageCompareTriplet();
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'coverage-evidence-bundles.json'), JSON.stringify({
    bundles: [
      {
        bundle_id: 'foo_compare_bundle',
        bundle_type: 'helper_surface',
        shared_test_files: ['tests/testthat/test-coverage-demo.R'],
        baseline_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/a.R',
            line_start: 1,
            line_end: 20
          }
        ],
        target_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/b.R',
            line_start: 1,
            line_end: 20
          }
        ]
      }
    ]
  }, null, 2));
  const baselineFixturePath = path.join(targetRoot, 'baseline-coverage-evidence-fixture.json');
  const targetFixturePath = path.join(targetRoot, 'target-coverage-evidence-fixture.json');
  writeFile(baselineFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: baselineRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/a.R',
        line_percent: 95,
        lines_total: 20,
        lines_covered: 19,
        covered_lines: Array.from({ length: 19 }, (_, idx) => idx + 1),
        uncovered_lines: [20]
      }
    ],
    line_percent: 95,
    lines_total: 20,
    lines_covered: 19,
    notes: ['baseline evidence fixture']
  }, null, 2));
  writeFile(targetFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: targetRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/b.R',
        line_percent: 60,
        lines_total: 20,
        lines_covered: 12,
        covered_lines: Array.from({ length: 12 }, (_, idx) => idx + 1),
        uncovered_lines: Array.from({ length: 8 }, (_, idx) => idx + 13)
      }
    ],
    line_percent: 60,
    lines_total: 20,
    lines_covered: 12,
    notes: ['target evidence fixture']
  }, null, 2));

  runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-compare',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-evidence-bundles.json'
  ], {
    env: {
      ...process.env,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_BASELINE: baselineFixturePath,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_TARGET: targetFixturePath
    }
  });

  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-evidence',
    '--repo-root',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-evidence-bundles.json'
  ]);

  assert.equal(summary.mode, 'coverage_evidence');
  assert.equal(summary.status, 'completed_with_gaps');
  assert.equal(summary.bundle_count, 1);
  assert.equal(summary.low_coverage_bundle_count, 1);
  assert.equal(summary.regression_bundle_count, 1);
  assert.equal(summary.justified_bundle_count, 1);
  assert.equal(summary.package_coverage.target_line_percent_display, '60.0%');
  assert.equal(summary.package_coverage.gate_status, 'below_target');
  assert.equal(summary.gate_counts.below_target, 1);
  assert.equal(summary.delta_counts.regressed, 1);
  assert.equal(summary.exception_count, 2);
  assert.equal(summary.below_target_bundles[0].bundle_id, 'foo_compare_bundle');
  assert.equal(summary.below_target_bundles[0].threshold_display, '90.0%');

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const runModes = querySqlite(auditDb, 'select mode, count(*) from audit_runs group by mode order by mode');
  assert.deepEqual(runModes, [
    ['coverage_compare', 1],
    ['coverage_evidence', 1]
  ]);

  const exceptionRows = querySqlite(
    auditDb,
    'select audit_layer, exception_type, severity from exceptions order by exception_type'
  );
  assert.deepEqual(exceptionRows, [
    ['coverage_evidence', 'coverage_target_below_threshold', 'medium'],
    ['coverage_evidence', 'package_coverage_below_threshold', 'medium']
  ]);

  const evidenceBundleRows = querySqlite(
    auditDb,
    `select side, coverage_gate_status
       from coverage_bundles
      where run_id = (select run_id from audit_runs where mode = "coverage_evidence")`
  );
  assert.deepEqual(evidenceBundleRows, [['evidence', 'below_target']]);

  const evidenceResultRow = querySqlite(
    auditDb,
    `select status, round(line_percent, 1), result_json
       from coverage_bundle_results
      where run_id = (select run_id from audit_runs where mode = "coverage_evidence")`
  )[0];
  assert.equal(evidenceResultRow[0], 'below_target');
  assert.equal(evidenceResultRow[1], 60.0);
  const evidencePayload = JSON.parse(evidenceResultRow[2]);
  assert.equal(evidencePayload.threshold_pct.toFixed(1), '90.0');
  assert.equal(evidencePayload.regressed_vs_baseline, true);

  const latestCoverageEvidence = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage-evidence.json'),
    'utf8'
  ));
  assert.equal(latestCoverageEvidence.status, 'completed_with_gaps');
  assert.equal(latestCoverageEvidence.low_coverage_bundle_count, 1);

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage-evidence.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Below Target Bundles/);
  assert.match(summaryMarkdown, /Regressed Bundles/);
});

serialTest('coverage-evidence flags baseline-target comparison gaps and suppresses package gate on subset runs', () => {
  const { baselineRoot, targetRoot } = makeCoverageCompareTriplet();
  writeFile(path.join(targetRoot, 'tools', 'refactor', 'coverage-evidence-comparison-bundles.json'), JSON.stringify({
    bundles: [
      {
        bundle_id: 'foo_compare_gap',
        bundle_type: 'wrapper_entrypoint',
        shared_test_files: ['tests/testthat/test-coverage-demo.R'],
        baseline_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/a.R',
            line_start: 1,
            line_end: 20
          }
        ],
        target_members: [
          {
            entity_key: 'symbol::foo',
            file_path: 'R/b.R',
            line_start: 1,
            line_end: 20
          }
        ]
      },
      {
        bundle_id: 'unused_compare_gap',
        bundle_type: 'wrapper_entrypoint',
        shared_test_files: ['tests/testthat/test-coverage-demo-extra.R'],
        baseline_members: [
          {
            entity_key: 'symbol::bar',
            file_path: 'R/a.R',
            line_start: 21,
            line_end: 25
          }
        ],
        target_members: [
          {
            entity_key: 'symbol::bar',
            file_path: 'R/b.R',
            line_start: 21,
            line_end: 25
          }
        ]
      }
    ]
  }, null, 2));
  const baselineFixturePath = path.join(targetRoot, 'baseline-coverage-comparison-gap-fixture.json');
  const targetFixturePath = path.join(targetRoot, 'target-coverage-comparison-gap-fixture.json');
  writeFile(baselineFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: baselineRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/a.R',
        line_percent: 0,
        lines_total: 20,
        lines_covered: 0,
        covered_lines: [],
        uncovered_lines: Array.from({ length: 20 }, (_, idx) => idx + 1)
      }
    ],
    line_percent: 0,
    lines_total: 20,
    lines_covered: 0,
    notes: ['baseline comparison-gap fixture']
  }, null, 2));
  writeFile(targetFixturePath, JSON.stringify({
    status: 'completed',
    tool_status: 'ok',
    project_root: targetRoot,
    test_files: ['tests/testthat/test-coverage-demo.R'],
    files: [
      {
        file_path: 'R/b.R',
        line_percent: 95,
        lines_total: 20,
        lines_covered: 19,
        covered_lines: Array.from({ length: 19 }, (_, idx) => idx + 1),
        uncovered_lines: [20]
      }
    ],
    line_percent: 95,
    lines_total: 20,
    lines_covered: 19,
    notes: ['target comparison-gap fixture']
  }, null, 2));

  const compareSummary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-compare',
    '--repo-root',
    targetRoot,
    '--baseline-path',
    baselineRoot,
    '--target-path',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-evidence-comparison-bundles.json',
    '--bundle-id',
    'foo_compare_gap'
  ], {
    env: {
      ...process.env,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_BASELINE: baselineFixturePath,
      FIDELITY_COVERAGE_COMPARE_FIXTURE_JSON_TARGET: targetFixturePath
    }
  });

  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'coverage-evidence',
    '--repo-root',
    targetRoot,
    '--bundle-manifest',
    'tools/refactor/coverage-evidence-comparison-bundles.json',
    '--coverage-compare-run-id',
    compareSummary.run_id,
    '--bundle-id',
    'foo_compare_gap'
  ]);

  assert.equal(summary.mode, 'coverage_evidence');
  assert.equal(summary.status, 'completed_with_gaps');
  assert.equal(summary.bundle_count, 1);
  assert.equal(summary.low_coverage_bundle_count, 0);
  assert.equal(summary.comparison_gap_bundle_count, 1);
  assert.equal(summary.regression_bundle_count, 0);
  assert.equal(summary.exception_count, 1);
  assert.equal(summary.package_coverage.gate_status, 'not_applicable_subset');
  assert.equal(summary.gate_counts.comparison_gap, 1);
  assert.equal(summary.comparison_gap_bundles[0].bundle_id, 'foo_compare_gap');
  assert.equal(summary.comparison_gap_bundles[0].target_line_percent_display, '95.0%');
  assert.equal(summary.comparison_gap_bundles[0].baseline_line_percent_display, '0.0%');

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const exceptionRows = querySqlite(
    auditDb,
    'select audit_layer, exception_type, severity from exceptions order by exception_type'
  );
  assert.deepEqual(exceptionRows, [
    ['coverage_evidence', 'coverage_comparison_gap', 'medium']
  ]);

  const evidenceResultRow = querySqlite(
    auditDb,
    `select status, round(line_percent, 1), result_json
       from coverage_bundle_results
      where run_id = (select run_id from audit_runs where mode = "coverage_evidence")`
  )[0];
  assert.equal(evidenceResultRow[0], 'comparison_gap');
  assert.equal(evidenceResultRow[1], 95.0);
  const evidencePayload = JSON.parse(evidenceResultRow[2]);
  assert.equal(evidencePayload.threshold_pct.toFixed(1), '80.0');
  assert.equal(evidencePayload.baseline_line_percent.toFixed(1), '0.0');
  assert.equal(evidencePayload.target_line_percent.toFixed(1), '95.0');

  const latestCoverageEvidence = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage-evidence.json'),
    'utf8'
  ));
  assert.equal(latestCoverageEvidence.package_coverage.gate_status, 'not_applicable_subset');
  assert.equal(latestCoverageEvidence.comparison_gap_bundle_count, 1);

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-coverage-evidence.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Comparison Gap Bundles/);
});

serialTest('bundle-map command collapses curated parity evidence into ranked coverage bundles', () => {
  const { baselineRoot, mainlineRoot, targetRoot } = makeBundleMapTriplet();
  const closeout = runJson('python3', [
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
    'tools/refactor/manifest-closeout-bundle-map.yml',
    '--case-path',
    'tools/refactor/behavior-cases-closeout-bundle-map.json',
    '--test-file',
    'tests/testthat/test-overlap-helper-characterization.R',
    '--contracts-execute'
  ]);

  assert.equal(closeout.status, 'ready_for_coverage_campaign');

  const summary = runJson('python3', [
    AUDIT_SCRIPT,
    'bundle-map',
    '--repo-root',
    targetRoot
  ]);

  assert.equal(summary.mode, 'bundle_map');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.bundle_count, 3);
  assert.equal(summary.linked_exception_count, 4);
  assert.equal(summary.bundle_type_counts.helper_surface, 1);
  assert.equal(summary.bundle_type_counts.lineage_family, 2);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const bundleRows = querySqlite(
    auditDb,
    'select bundle_type, linked_exception_count, priority_rank from coverage_bundles where run_id = (select run_id from audit_runs where mode = \"bundle_map\") order by priority_rank'
  );
  assert.deepEqual(bundleRows, [
    ['helper_surface', 2, 1],
    ['lineage_family', 1, 2],
    ['lineage_family', 1, 3]
  ]);

  const exceptionLinkRows = querySqlite(
    auditDb,
    'select count(*) from coverage_exception_links where run_id = (select run_id from audit_runs where mode = \"bundle_map\")'
  );
  assert.equal(exceptionLinkRows[0][0], 4);

  const testLinkRows = querySqlite(
    auditDb,
    'select test_file_path from coverage_test_links where run_id = (select run_id from audit_runs where mode = \"bundle_map\") order by test_file_path'
  );
  assert.deepEqual(testLinkRows, [['tests/testthat/test-overlap-helper-characterization.R']]);

  const bundleManifest = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-bundles.json'),
    'utf8'
  ));
  assert.equal(bundleManifest.bundle_count, 3);
  assert.equal(bundleManifest.bundles[0].linked_exception_count, 2);
  assert.equal(bundleManifest.bundles[0].shared_test_files[0], 'tests/testthat/test-overlap-helper-characterization.R');
  const lineageBundles = bundleManifest.bundles.filter((bundle) => bundle.bundle_type === 'lineage_family');
  lineageBundles.sort((left, right) => {
    const leftStart = left.target_members[0]?.line_start ?? 0;
    const rightStart = right.target_members[0]?.line_start ?? 0;
    return leftStart - rightStart;
  });
  assert.equal(lineageBundles.length, 2);
  assert.deepEqual(lineageBundles.map((bundle) => bundle.baseline_paths), [['R/source_gap.R'], ['R/source_gap.R']]);
  assert.deepEqual(lineageBundles.map((bundle) => bundle.target_paths), [['R/source_gap_target.R'], ['R/source_gap_target.R']]);
  assert.deepEqual(
    lineageBundles.map((bundle) => bundle.target_members),
    [
      [
        {
          entity_key: 'symbol::gapHelper',
          file_path: 'R/source_gap_target.R',
          line_start: 1,
          line_end: 3,
          metadata: {
            side: 'target',
            selector_kind: 'symbol',
            comparison_status: 'source_missing_target_resolved'
          }
        }
      ],
      [
        {
          entity_key: 'symbol::gapObserverRegister',
          file_path: 'R/source_gap_target.R',
          line_start: 5,
          line_end: 7,
          metadata: {
            side: 'target',
            selector_kind: 'symbol',
            comparison_status: 'source_missing_target_resolved'
          }
        }
      ]
    ]
  );
  assert.deepEqual(lineageBundles.map((bundle) => bundle.baseline_members), [[], []]);

  const summaryMarkdown = fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-bundles.md'),
    'utf8'
  );
  assert.match(summaryMarkdown, /Priority Bundles/);
  assert.match(summaryMarkdown, /helper_surface/);
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
  assert.equal(summary.diff_counts.collate_drift || 0, 0);
  assert.equal(summary.exception_count, 5);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const fileSurfaceCount = querySqlite(auditDb, 'select count(*) from inventory_files')[0][0];
  const exportSurfaceCount = querySqlite(auditDb, 'select count(*) from inventory_exports')[0][0];
  const diffCount = querySqlite(auditDb, 'select count(*) from inventory_diffs')[0][0];
  const exceptionCount = querySqlite(auditDb, 'select count(*) from exceptions')[0][0];
  assert.equal(fileSurfaceCount, 5);
  assert.equal(exportSurfaceCount, 6);
  assert.equal(diffCount, 11);
  assert.equal(exceptionCount, 5);

  const diffClasses = querySqlite(
    auditDb,
    'select diff_class, entity_key from inventory_diffs order by diff_class, entity_key'
  );
  assert.deepEqual(diffClasses, [
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
  assert.match(summaryMarkdown, /Total diffs: `11`/);
  assert.match(summaryMarkdown, /Exception candidates: `5`/);
  assert.match(summaryMarkdown, /`missing_definition`: `1`/);
  assert.match(summaryMarkdown, /`moved_definition`: `3`/);
});

serialTest('surface command ignores symmetric duplicate definitions when both sides share the same duplicate shape', () => {
  const { baselineRoot, targetRoot } = makeSurfaceParityPair();
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

  assert.equal(summary.total_diffs, 0);
  assert.equal(summary.open_diff_count, 0);
  assert.equal(summary.exception_count, 0);
});

serialTest('surface command ignores duplicate-shape reductions when semantic fingerprints are preserved', () => {
  const { baselineRoot, targetRoot } = makeSurfaceDedupSafePair();
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

  assert.equal(summary.total_diffs, 0);
  assert.equal(summary.open_diff_count, 0);
  assert.equal(summary.exception_count, 0);
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
  assert.equal(summary.resolver_counts.selector, 4);
  assert.equal(summary.resolver_counts.raw_text_search, 1);
  assert.equal(summary.resolver_counts.file_normalized_text_search || 0, 0);
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
    ['anchor_normalized', 'normalized_only', 'selector'],
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
  assert.equal(summary.status_counts.source_missing_target_resolved, 1);
  assert.equal(summary.resolver_counts.selector, 1);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const comparisonRows = querySqlite(
    auditDb,
    'select entry_id, status, target_resolver from manifest_entries join manifest_comparisons using (manifest_entry_id)'
  );
  assert.deepEqual(comparisonRows, [[
    'source_gap',
    'source_missing_target_resolved',
    'selector'
  ]]);

  const exceptionRows = querySqlite(
    auditDb,
    'select exception_type, severity, reason from exceptions order by exception_type'
  );
  assert.deepEqual(exceptionRows, [[
    'source_lineage_gap_target_resolved',
    'low',
    'Baseline/source block could not be resolved via `selector:symbol`, but the target block resolved via `selector`: baseline/source selector no longer resolves directly; target block resolved via selector'
  ]]);

  const latestExceptions = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-exceptions.json'),
    'utf8'
  ));
  assert.equal(latestExceptions.summary.open_count, 1);
  assert.equal(latestExceptions.exceptions[0].exception_type, 'source_lineage_gap_target_resolved');
});

serialTest('manifest command resolves stale target paths through repo-wide selector fallback', () => {
  const { baselineRoot, targetRoot } = makeManifestRepoSelectorTargetPair();
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
    'tools/refactor/manifest-repo-selector.yml'
  ]);

  assert.equal(summary.mode, 'manifest');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.entry_count, 1);
  assert.equal(summary.exception_count, 1);
  assert.equal(summary.status_counts.raw_exact, 1);
  assert.equal(summary.resolver_counts.repo_selector_search, 1);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const comparisonRows = querySqlite(
    auditDb,
    'select entry_id, status, target_resolver from manifest_entries join manifest_comparisons using (manifest_entry_id)'
  );
  assert.deepEqual(comparisonRows, [[
    'symbol_repo_selector',
    'raw_exact',
    'repo_selector_search'
  ]]);

  const exceptionRows = querySqlite(
    auditDb,
    'select exception_type, severity from exceptions order by exception_type'
  );
  assert.deepEqual(exceptionRows, [[
    'target_resolution_asymmetry',
    'low'
  ]]);
});

serialTest('manifest command supports asymmetric target selectors for moved anchor ranges', () => {
  const { baselineRoot, targetRoot } = makeManifestTargetSelectorPair();
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
    'tools/refactor/manifest-target-selector.yml'
  ]);

  assert.equal(summary.mode, 'manifest');
  assert.equal(summary.status, 'completed');
  assert.equal(summary.entry_count, 1);
  assert.equal(summary.exception_count, 0);
  assert.equal(summary.status_counts.raw_exact, 1);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const comparisonRows = querySqlite(
    auditDb,
    'select entry_id, status, target_resolver from manifest_entries join manifest_comparisons using (manifest_entry_id)'
  );
  assert.deepEqual(comparisonRows, [[
    'anchor_target_selector',
    'raw_exact',
    'selector'
  ]]);
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

serialTest('closeout auto-curates redundant surface drifts when manifest semantic review already covers the same entity', () => {
  const { baselineRoot, mainlineRoot, targetRoot } = makeCloseoutOverlapTriplet();
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
    'tools/refactor/manifest-closeout-overlap.yml',
    '--case-path',
    'tools/refactor/behavior-cases-closeout-overlap.json',
    '--test-file',
    'tests/testthat/test-prot-03-rollup.R',
    '--contracts-execute'
  ]);

  assert.equal(summary.status, 'blocked');
  assert.equal(summary.auto_curated_count, 1);
  assert.equal(summary.readiness.open_exception_count, 1);
  assert.equal(summary.readiness.high_open_exception_count, 1);
  assert.equal(summary.readiness.surface_gate, true);
  assert.equal(summary.readiness.manifest_gate, false);
  assert.equal(summary.readiness.open_exception_layer_counts.manifest_fidelity, 1);
  assert.equal(summary.readiness.open_exception_layer_counts.surface_inventory || 0, 0);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const exceptionRows = querySqlite(
    auditDb,
    "select exception_type, curation_status, disposition, curation_rule from exceptions order by exception_type"
  );
  assert.deepEqual(exceptionRows, [
    ['semantic_mismatch', 'candidate', null, null],
    ['surface_definition_drift', 'auto_curated', 'covered_by_manifest_semantic_review', 'auto:surface_definition_drift_manifest_overlap']
  ]);
});

serialTest('closeout applies persistent curation catalog entries by exception key', () => {
  const { baselineRoot, mainlineRoot, targetRoot } = makeCloseoutCatalogCuratedTriplet();
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
    'tools/refactor/manifest-closeout-curation.yml',
    '--case-path',
    'tools/refactor/behavior-cases-closeout-curation.json',
    '--test-file',
    'tests/testthat/test-prot-03-rollup.R',
    '--contracts-execute'
  ]);

  assert.equal(summary.status, 'ready_for_coverage_campaign');
  assert.equal(summary.proven_parity, true);
  assert.equal(summary.catalog_curated_count, 1);
  assert.equal(summary.readiness.open_exception_count, 0);
  assert.equal(summary.readiness.high_open_exception_count, 0);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const exceptionRows = querySqlite(
    auditDb,
    "select exception_type, curation_status, disposition, owner, curation_rule from exceptions order by exception_type"
  );
  assert.deepEqual(exceptionRows, [
    ['surface_duplicate_definition', 'accepted', 'canonical_duplicate_retired', 'fixture', 'catalog:tools/refactor/fidelity-curations.json:dup_helper_canonicalized']
  ]);

  const latestCloseout = JSON.parse(fs.readFileSync(
    path.join(targetRoot, '.refactor-fidelity-audit', 'reports', 'latest-closeout.json'),
    'utf8'
  ));
  assert.equal(latestCloseout.catalog_curated_count, 1);
});

serialTest('closeout clears redundant surface drift after catalog-curating the covering manifest mismatch', () => {
  const { baselineRoot, mainlineRoot, targetRoot } = makeCloseoutCatalogOverlapTriplet();
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
    'tools/refactor/manifest-closeout-catalog-overlap.yml',
    '--case-path',
    'tools/refactor/behavior-cases-closeout-catalog-overlap.json',
    '--test-file',
    'tests/testthat/test-prot-03-rollup.R',
    '--contracts-execute'
  ]);

  assert.equal(summary.status, 'ready_for_coverage_campaign');
  assert.equal(summary.proven_parity, true);
  assert.equal(summary.catalog_curated_count, 1);
  assert.equal(summary.auto_curated_count, 1);
  assert.equal(summary.readiness.open_exception_count, 0);
  assert.equal(summary.readiness.high_open_exception_count, 0);

  const auditDb = path.join(targetRoot, '.refactor-fidelity-audit', 'audit.db');
  const exceptionRows = querySqlite(
    auditDb,
    "select exception_type, curation_status, disposition, owner, curation_rule from exceptions order by exception_type"
  );
  assert.deepEqual(exceptionRows, [
    ['semantic_mismatch', 'accepted', 'entrypoint_shell_equivalent_under_tests', 'fixture', 'catalog:tools/refactor/fidelity-curations.json:overlap_helper_manifest_curated'],
    ['surface_definition_drift', 'auto_curated', 'covered_by_manifest_semantic_review', null, 'auto:surface_definition_drift_manifest_overlap']
  ]);
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
