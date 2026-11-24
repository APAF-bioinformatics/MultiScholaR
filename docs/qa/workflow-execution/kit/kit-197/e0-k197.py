# E0: the verbatim blocks of sim-k197.R against filter-effect.R at 6416dbd6 (text only).
import re, sys
sim = open(sys.argv[1]).read().split('\n')
fe = open(sys.argv[2]).read().split('\n')
i = 0; blocks = []; total_diff = []
while i < len(sim):
    m = re.match(r'# ==== verbatim: filter-effect\.R:(\d+)-(\d+) at 6416dbd6', sim[i])
    if m:
        a, b = int(m.group(1)), int(m.group(2)); j = i + 1
        while sim[j] != '# ==== end verbatim ====': j += 1
        body = sim[i+1:j]; ref = fe[a-1:b]
        diffs = [(a+k, ref[k], body[k]) for k in range(max(len(ref), len(body))) if k >= len(ref) or k >= len(body) or ref[k] != body[k]]
        blocks.append((a, b, len(body), len(ref), [d[0] for d in diffs])); total_diff += diffs; i = j
    i += 1
for bl in blocks: print('block %d-%d: %d lines (source %d), differing source lines: %s' % bl)
print('differing lines in total:', len(total_diff))
for d in total_diff: print('line', d[0]); print('  source:', d[1]); print('  copy:  ', d[2])
