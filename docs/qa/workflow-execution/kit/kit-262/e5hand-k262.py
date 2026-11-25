"""KIT-262 E5 hand values (local, standard library only): the Wilson limits and the gain threshold that score-k262.R's
self-test uses, computed independently. Run: python3 -I docs/qa/workflow-execution/kit/kit-262/e5hand-k262.py"""
import math
import sys
from statistics import NormalDist

nd = NormalDist()
z95 = nd.inv_cdf(0.975)
z_gain = nd.inv_cdf(1 - 0.05 / 180)


def wilson(x, n, z):
    p = x / n
    d = 1 + z * z / n
    centre = (p + z * z / (2 * n)) / d
    half = z * math.sqrt(p * (1 - p) / n + z * z / (4 * n * n)) / d
    return centre - half, centre + half


ok = True
print("z95 = %.6f; z_gain = %.6f" % (z95, z_gain))
ok &= abs(z_gain - 3.4524) < 5e-4
for x in (5, 6, 7, 8):
    lo, hi = wilson(x, 40, z_gain)
    print("one-sided gain limit at %d of 40: lower %.6f" % (x, lo))
ok &= wilson(7, 40, z_gain)[0] >= 0.05 > wilson(6, 40, z_gain)[0]
lo, hi = wilson(2, 40, z95)
print("two-sided 95%% at 2 of 40: %.4f to %.4f" % (lo, hi))
ok &= abs(lo - 0.0138) < 5e-4 and abs(hi - 0.1650) < 5e-4
# smallest count whose lower limit reaches 0.05, at R = 40: the ticket says 7
smallest = next(x for x in range(41) if wilson(x, 40, z_gain)[0] >= 0.05)
print("smallest count of 40 with lower limit >= 0.05: %d" % smallest)
ok &= smallest == 7
sys.exit(0 if ok else 1)
