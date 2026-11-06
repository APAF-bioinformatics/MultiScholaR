#!/bin/bash
# all.sh: KIT-100's three screen journeys in sequence (one attempt each, protocol v1.9), each committed and pushed by run.sh.
S=/home/doktersmol/Projects/APAF/MultiScholaR-dia-proof-20260911-VhC1n5/kit/scripts-100
for m in kit-journey-100a kit-journey-100b kit-journey-100c; do bash $S/run.sh $m > $S/run-${m#kit-journey-}.out 2>&1; tail -4 $S/run-${m#kit-journey-}.out; done
echo ALL DONE
