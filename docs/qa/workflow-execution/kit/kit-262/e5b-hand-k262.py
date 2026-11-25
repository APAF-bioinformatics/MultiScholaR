"""KIT-262 E5b hand values (local, standard library only): the constructed curves of e5b-k262.R, with the elbow, the
separation and the composite of each rule worked independently of the package's code.
Run: python3 -I docs/qa/workflow-execution/kit/kit-262/e5b-hand-k262.py"""
import math


def elbow(delta, eps=0.05, min_effect=0.05):
    finite = [(k + 1, d) for k, d in enumerate(delta) if math.isfinite(d)]
    mx = max(d for _, d in finite)
    if mx < min_effect:
        return 1
    return min(k for k, d in finite if d >= mx * (1 - eps))


def composite(sep, k, w, cap):
    if sep <= 0:
        return 0.0
    if cap == 1:
        pen = 0 if k <= 1 else 1
    elif k <= cap:
        pen = w * (k - 1) / (cap - 1)
    else:
        pen = w + (1 - w) * (1 - math.exp(-(k - cap)))
    pen = max(0.0, min(1.0, pen))
    return sep * (1 - pen)


cases = {
    # name: (cap, [(percentage, delta curve, c)])
    "A_readmitted_outranks": (2, [(5, [.1, .2, .3, .4, .5, .6, .6, .6], 8), (10, [.1, .2, .3, .3, .3, .3, .3, .3], 2)]),
    "B_readmitted_does_not": (2, [(5, [.3, .6, .9, .9, .9, .9], 8), (10, [.0, .05, .1, .15, .2, .2], 2)]),
    "C_truncated_elbow_below": (3, [(5, [.5, .5, .5, 1.0, 1.0, 1.0], 6)]),
    "D_truncated_separation_nonpositive": (2, [(5, [-.1, -.05, .4, .4, .4, .4], 6)]),
}
for name, (cap, rows) in cases.items():
    print(name, "cap", cap)
    for pct, delta, c in rows:
        k = elbow(delta)
        sep = max(delta)
        kb = elbow(delta[:cap])
        sepb = max(delta[:cap])
        print("  pct %d: k %d sep %.4f c %d | P usable %s composite %.6f | C1 usable %s composite %.6f k_used %d | C2 k_b %d sep_b %.4f composite_b %.6f usable %s"
              % (pct, k, sep, c, k <= c, composite(sep, k, .5, cap), min(k, cap) <= c, composite(sep, k, .5, cap), min(k, cap),
                 kb, sepb, composite(sepb, kb, .5, cap), kb <= c))
