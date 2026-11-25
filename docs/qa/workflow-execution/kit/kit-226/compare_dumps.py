"""KIT-226 control C: compare the W1 dumps (R code text) and the prot_dia GUI-entry records leaf by leaf.

Reads text only (no R). Usage:
  python3 -I compare_dumps.py <stored test file> <stored prot_dia.json> <base dir> <head dir 1> <head dir 2> <out tsv>
A leaf is a non-list value (each element of a c() vector is its own leaf); its path is the list names from `omicEnrichW1Frozen[["<case>"]]` down (W1), or the JSON
keys from `steps` down (GUI entry). Writes one row per leaf that differs anywhere: source, path, stored, base, head1,
head2, base_eq_stored, heads_equal, changed.
"""
import json
import re
import sys


def tokens(text):
    pat = re.compile(r'\s+|#[^\n]*|"(?:\\.|[^"\\])*"|`[^`]*`|[A-Za-z_.][A-Za-z0-9_.:]*|-?[0-9][0-9.eE+-]*L?|<-|[()\[\],=]|.')
    pos = 0
    out = []
    line = 1
    while pos < len(text):
        m = pat.match(text, pos)
        tok = m.group(0)
        if not tok.isspace() and not tok.startswith('#'):
            out.append((tok, line))
        line += tok.count('\n')
        pos = m.end()
    return out


class Parser:
    def __init__(self, toks):
        self.t = toks
        self.i = 0

    def peek(self):
        return self.t[self.i][0] if self.i < len(self.t) else None

    def take(self):
        tok = self.t[self.i]
        self.i += 1
        return tok

    def expr(self):
        """Returns ('list', [(name, value)], line) or ('leaf', text, line)."""
        tok, line = self.take()
        # A list, or a c() vector whose elements are leaves of their own (named or by position).
        if tok in ('list', 'c') and self.peek() == '(':
            self.take()
            items = []
            while self.peek() != ')':
                name = ''
                if self.i + 1 < len(self.t) and self.t[self.i + 1][0] == '=':
                    name = self.take()[0].strip('`')
                    self.take()
                items.append((name, self.expr()))
                if self.peek() == ',':
                    self.take()
            self.take()
            return ('list', items, line)
        parts = [tok]
        if self.peek() == '(':
            depth = 0
            while True:
                t, _ = self.take()
                parts.append(t)
                if t == '(':
                    depth += 1
                elif t == ')':
                    depth -= 1
                    if depth == 0:
                        break
        return ('leaf', ' '.join(parts), line)


def flatten(node, path, out):
    if node[0] == 'list':
        if not node[1]:
            out[path] = ('list()', node[2])
        for k, (name, value) in enumerate(node[1]):
            flatten(value, path + ('%s' % (name or '[%d]' % (k + 1)),), out)
    else:
        out[path] = (node[1], node[2])


def w1_cases(text):
    """{case: {path: (leaf, line)}} for each `omicEnrichW1Frozen[["case"]] <- list(...)` in the text."""
    cases = {}
    for m in re.finditer(r'^omicEnrichW1Frozen\[\["([^"]+)"\]\] <- ', text, re.M):
        start = m.end()
        offset_line = text.count('\n', 0, start)
        p = Parser([(t, l + offset_line) for t, l in tokens(text[start:])])
        node = p.expr()
        flat = {}
        flatten(node, (m.group(1),), flat)
        cases[m.group(1)] = flat
    return cases


def json_flat(x, path, out):
    if isinstance(x, dict):
        if not x:
            out[path] = '{}'
        for k, v in x.items():
            json_flat(v, path + (k,), out)
    elif isinstance(x, list):
        if not x:
            out[path] = '[]'
        for k, v in enumerate(x):
            json_flat(v, path + ('[%d]' % (k + 1),), out)
    else:
        out[path] = json.dumps(x)


def main():
    stored_r, stored_json, base, head1, head2, out = sys.argv[1:7]
    rows = []
    stored = {}
    for case in w1_cases(open(stored_r).read()).values():
        stored.update(case)
    dumps = []
    for d in (base, head1, head2):
        flat = {}
        for case in w1_cases(open(d + '/test-omic-enrich-module-kit107.R').read()).values():
            flat.update(case)
        dumps.append(flat)
    keys = sorted(set(stored) | set().union(*[set(x) for x in dumps]))
    for k in keys:
        vals = [stored.get(k, ('<absent>', 0))[0]] + [x.get(k, ('<absent>', 0))[0] for x in dumps]
        if len(set(vals)) > 1:
            rows.append(['W1', '/'.join(k), str(stored.get(k, ('', ''))[1])] + vals)
    oracle = {}
    json_flat(json.load(open(stored_json))['steps'], (), oracle)
    recs = []
    for d in (base, head1, head2):
        flat = {}
        json_flat(json.load(open(d + '/gui-entry/prot_dia.json'))['steps'], (), flat)
        recs.append(flat)
    keys = sorted(set(oracle) | set().union(*[set(x) for x in recs]))
    for k in keys:
        vals = [oracle.get(k, '<absent>')] + [x.get(k, '<absent>') for x in recs]
        if len(set(vals)) > 1:
            rows.append(['prot_dia', '/'.join(k), ''] + vals)
    with open(out, 'w') as f:
        f.write('\t'.join(['source', 'path', 'stored_line', 'stored', 'base', 'head1', 'head2', 'base_eq_stored',
                           'heads_equal', 'changed']) + '\n')
        for r in rows:
            stored_v, b, h1, h2 = r[3:7]
            f.write('\t'.join(r + [str(b == stored_v), str(h1 == h2), str(b != h1)]) + '\n')
    print(len(rows), 'differing leaves written to', out)


if __name__ == '__main__':
    main()
