# E0: the verbatim blocks of sim-k211.R, score-k211.R and l1-k211.R against KIT-197's scripts at ef63c918 (text only; git
# and the standard library). Usage, from the package root: python3 -I docs/qa/workflow-execution/kit/kit-211/e0-k211.py
# sim-k211.R holds KIT-197's sim-k197.R lines 26 to 279 unchanged except two inserted lines (Av's guarded variance line
# after the variance draw; Ag's guarded planting line directly after `beta[where] <- signs`); the blocks of score-k211.R
# and l1-k211.R (markers "verbatim: <file>:<a>-<b>") equal the lines a to b of KIT-197's file exactly.
import difflib, hashlib, re, subprocess, sys

COMMIT = 'ef63c9180'
KIT197 = 'docs/qa/workflow-execution/kit/kit-197/'
KIT211 = 'docs/qa/workflow-execution/kit/kit-211/'
MARKED = ['    if (identical(cellname, "Av")) s2 <- s2 * 2^((20 - mu) / 2)',
          '    if (identical(cellname, "Ag")) beta <- 3 * beta']


def old(name):
    return subprocess.check_output(['git', 'show', '%s:%s%s' % (COMMIT, KIT197, name)]).decode().split('\n')


def sha(lines):
    return hashlib.sha256(('\n'.join(lines) + '\n').encode()).hexdigest()


ok = True
# sim-k211.R
new = open(KIT211 + 'sim-k211.R').read().split('\n')
first = next(i for i, l in enumerate(new) if l.startswith('# ==== verbatim: filter-effect.R:23-234'))
last = max(i for i, l in enumerate(new) if l == '# ==== end verbatim ====')
region = new[first:last + 1]
ref = old('sim-k197.R')[25:279]
ops = [o for o in difflib.SequenceMatcher(None, ref, region, autojunk=False).get_opcodes() if o[0] != 'equal']
inserted = [l for o in ops if o[0] == 'insert' for l in region[o[3]:o[4]]]
good = all(o[0] == 'insert' for o in ops) and inserted == MARKED
ok &= good
print('sim-k211.R: lines %d to %d of the file = sim-k197.R lines 26 to 279 (%d lines) plus %d inserted: %s'
      % (first + 1, last + 1, len(ref), len(inserted), 'as declared' if good else 'NOT as declared'))
for l in inserted:
    print('  inserted:', l)
print('  SHA-256 of the copy (newline-terminated):', sha([l for l in region if l not in MARKED]))
print('  SHA-256 of the block with its two marked lines:', sha(region))
# score-k211.R and l1-k211.R
for fname, orig_name in (('score-k211.R', 'score-k197.R'), ('l1-k211.R', 'l1-k197.R')):
    lines = new_lines = open(KIT211 + fname).read().split('\n')
    ref = old(orig_name)
    i = 0
    n = 0
    while i < len(lines):
        m = re.match(r'# ==== verbatim: %s:(\d+)-(\d+) at ef63c918 ====' % re.escape(orig_name), lines[i])
        if m:
            a, b = int(m.group(1)), int(m.group(2))
            j = i + 1
            while lines[j] != '# ==== end verbatim ====':
                j += 1
            body = lines[i + 1:j]
            same = body == ref[a - 1:b]
            ok &= same
            n += 1
            print('%s: block %s-%s (%d lines): %s; SHA-256 %s' % (fname, a, b, len(body), 'identical' if same else 'DIFFERS', sha(body)))
            i = j
        i += 1
    ok &= n > 0
    print('%s: %d verbatim blocks' % (fname, n))
print('E0:', 'pass' if ok else 'FAIL')
sys.exit(0 if ok else 1)
