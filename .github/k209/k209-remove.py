# KIT-209: remove one journey's folder after the journey (disk), guarded: the folder must be a direct child of the runs
# folder and named as the plan names journeys (<unit>_<build>_r<rep>[_<inject>]); anything else is refused.
# KIT-150's k150-remove.py, adapted. Usage: python3 -I k209-remove.py <runs folder> <journey folder>
import os
import re
import shutil
import sys

root = os.path.realpath(sys.argv[1])
target = os.path.realpath(sys.argv[2])
name = os.path.basename(target)
pattern = r"(U[0-9]{1,2}|S[12]|CAL|D[1-7]|X10)_(base|C)_r[0-9]+(_(fail_memory|fail_time|package_error))?"
if os.path.dirname(target) != root or not re.fullmatch(pattern, name):
    sys.exit(f"k209-remove: refused {target} (not a journey folder under {root})")
if os.path.isdir(target):
    shutil.rmtree(target)
print(f"k209-remove: removed {name}")
