# KIT-150: remove one run's folder after its stages (disk ceiling), guarded: the folder must be a direct child of the
# runs folder and named as the plan names runs (<size>_<build>_r<rep>[_failures]); anything else is refused.
# Usage: python3 -I k150-remove.py <runs folder> <run folder>
import os
import re
import shutil
import sys

root = os.path.realpath(sys.argv[1])
target = os.path.realpath(sys.argv[2])
name = os.path.basename(target)
if os.path.dirname(target) != root or not re.fullmatch(r"(S|M|L|Lp|P1|P2|P3)_(base|A|B1|B3|C|RC)_r[0-9]+(_failures)?", name):
    sys.exit(f"k150-remove: refused {target} (not a run folder under {root})")
if os.path.isdir(target):
    shutil.rmtree(target)
print(f"k150-remove: removed {name}")
