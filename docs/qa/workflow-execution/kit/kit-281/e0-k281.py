# E0 (KIT-281): the verbatim parts of sim-k281.R, score-k281.R and l1-k281.R against KIT-211's scripts on staging and KIT-197's at
# ef63c918 (text only; git and the standard library). Usage, from the package root:
#   python3 -I docs/qa/workflow-execution/kit/kit-281/e0-k281.py
# Blocks are located by their MARKERS, not by line numbers (the headers differ in length). It asserts only identities of code the
# ticket declares unchanged (section 6); everything else is declared new or adapted code and is listed in script-diff.md.
#   sim:   the region from the first "# ==== verbatim:" marker to the last "# ==== end verbatim ====" marker (KIT-211 :24-279,
#          which is KIT-197's lines plus KIT-211's two guarded lines) is byte-identical; and these functions, by name, are identical:
#          bFilterKept, packageBKept, lfcDeviation, aggregateSubset, writeTsv, num, readPart, partPaths, writeTask, assemble, sha256.
#   score: every "verbatim: score-k197.R:a-b" block equals KIT-197's lines a to b and KIT-211's block; these functions, by name, equal
#          KIT-211's: limitsOf, q1, fmt, mdTable, resampleStats, summariseStats, bootstrapWeights, jackSe, quant, num, writeTsv.
#   l1:    every "verbatim: l1-k197.R:a-b" block equals KIT-197's lines a to b.
import hashlib, re, subprocess, sys

COMMIT = 'ef63c9180'
KIT197 = 'docs/qa/workflow-execution/kit/kit-197/'
KIT211 = 'docs/qa/workflow-execution/kit/kit-211/'
KIT281 = 'docs/qa/workflow-execution/kit/kit-281/'
SHA211 = {'sim-k211.R': '51c2d4672fc7e758b565c4e74e886eac5f805751871e0a47feaa382c518011b4',
          'score-k211.R': 'aefc66baa7dc42a86f5f69247911e51b8701200b1e7d95851bec676ae1110d54',
          'l1-k211.R': 'a387ac1e74d44a054106b235a25848af71f941b3b2da903bed24e1123baba374'}
ok = True


def read(path):
    return open(path, encoding='utf-8').read()


def old197(name):
    return subprocess.check_output(['git', 'show', '%s:%s%s' % (COMMIT, KIT197, name)]).decode().split('\n')


def sha_text(text):
    return hashlib.sha256(text.encode()).hexdigest()


def function_block(lines, name):
    """The lines of `name <- function ...`: from its first line to the line before the next top-level statement or comment,
    trailing blank lines stripped."""
    start = next((i for i, l in enumerate(lines) if l.startswith(name + ' <- function')), None)
    if start is None:
        return None
    end = start + 1
    while end < len(lines) and not re.match(r'^[^\s})]', lines[end]):
        end += 1
    block = lines[start:end]
    while block and block[-1] == '':
        block.pop()
    return block


def verbatim_region(lines):
    first = next(i for i, l in enumerate(lines) if l.startswith('# ==== verbatim:'))
    last = max(i for i, l in enumerate(lines) if l == '# ==== end verbatim ====')
    return lines[first:last + 1]


def named_blocks(lines, prefix):
    """(marker, a, b, body) of every '# ==== verbatim: <prefix>:a-b at ef63c918 ====' block."""
    out, i = [], 0
    while i < len(lines):
        m = re.match(r'# ==== verbatim: %s:(\d+)-(\d+) at ef63c918 ====' % re.escape(prefix), lines[i])
        if m:
            j = i + 1
            while lines[j] != '# ==== end verbatim ====':
                j += 1
            out.append((int(m.group(1)), int(m.group(2)), lines[i + 1:j]))
            i = j
        i += 1
    return out


# KIT-211's files are the declared sources: their SHA-256 must be the ticket's
for name, want in SHA211.items():
    got = sha_text(read(KIT211 + name))
    good = got == want
    ok &= good
    print('%s: SHA-256 %s... %s' % (name, got[:16], 'as declared' if good else 'NOT as declared (%s...)' % want[:16]))

# ---- sim-k281.R ----
s211 = read(KIT211 + 'sim-k211.R').split('\n')
s281 = read(KIT281 + 'sim-k281.R').split('\n')
r211, r281 = verbatim_region(s211), verbatim_region(s281)
good = r211 == r281
ok &= good
print('sim-k281.R: verbatim region %d lines (KIT-211 :24-279) %s' % (len(r281), 'byte-identical' if good else 'DIFFERS'))
print('  SHA-256 of the region (newline-terminated):', sha_text('\n'.join(r281) + '\n'))
for name in ['bFilterKept', 'packageBKept', 'lfcDeviation', 'aggregateSubset', 'writeTsv', 'num', 'readPart', 'partPaths',
             'writeTask', 'assemble', 'sha256']:
    a, b = function_block(s211, name), function_block(s281, name)
    good = a is not None and a == b
    ok &= good
    print('  sim function %-16s %s' % (name, 'identical' if good else 'DIFFERS or absent'))

# ---- score-k281.R ----
c197 = old197('score-k197.R')
c211 = read(KIT211 + 'score-k211.R').split('\n')
c281 = read(KIT281 + 'score-k281.R').split('\n')
blocks = named_blocks(c281, 'score-k197.R')
ok &= len(blocks) == 3  # num+writeTsv (one block), jackSe, quant: the four functions of the ticket
print('score-k281.R: %d verbatim blocks (3 expected: num and writeTsv share one)' % len(blocks))
for a, b, body in blocks:
    good = body == c197[a - 1:b] and body in [x[2] for x in named_blocks(c211, 'score-k197.R')]
    ok &= good
    print('  block score-k197.R:%d-%d (%d lines): %s' % (a, b, len(body), 'identical to KIT-197 and KIT-211' if good else 'DIFFERS'))
for name in ['limitsOf', 'q1', 'fmt', 'mdTable', 'resampleStats', 'summariseStats', 'bootstrapWeights', 'jackSe', 'quant', 'num', 'writeTsv']:
    a, b = function_block(c211, name), function_block(c281, name)
    good = a is not None and a == b
    ok &= good
    print('  score function %-16s %s' % (name, 'identical' if good else 'DIFFERS or absent'))

# ---- l1-k281.R ----
l197 = old197('l1-k197.R')
l281 = read(KIT281 + 'l1-k281.R').split('\n')
blocks = named_blocks(l281, 'l1-k197.R')
ok &= len(blocks) >= 1
print('l1-k281.R: %d verbatim blocks' % len(blocks))
for a, b, body in blocks:
    good = body == l197[a - 1:b]
    ok &= good
    print('  block l1-k197.R:%d-%d (%d lines): %s' % (a, b, len(body), 'identical' if good else 'DIFFERS'))
print('E0:', 'pass' if ok else 'FAIL')
sys.exit(0 if ok else 1)
