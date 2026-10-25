#!/bin/zsh
# Bench: clean the previous LO's IQ, capture + analyse the next LO, print the flush log, dither map and grouped signals.
#   bench/live_lo/step.sh <prev_MHz> <MHz> [lo_step args...]
cd "$(dirname "$0")/../.."
prev=$1; f=$2; shift 2
for x in $prev $((prev + 1)); do .venv/bin/python bench/live_lo/lo_step.py --clean $x >/dev/null 2>&1; done
log=/private/tmp/claude-501/-Users-srikanth-Documents-modules-spectrum-intelligence/c94d801a-125c-43e1-8ae4-81b504e1452f/scratchpad/lo.log
nice -n 19 .venv/bin/python bench/live_lo/lo_step.py --freq $f --ppm 10 "$@" > $log 2>&1
grep '"gain_db"' $log | .venv/bin/python -c "import json,sys; d=json.loads(sys.stdin.read()); print('gain', d['gain_db'], 'rms', [round(b['rms']) for b in d['buffers']], 'clip', [b.get('clip') for b in d['buffers'] if 'clip' in b])"
grep -E "^\s+[0-9]{3,4}\.[0-9]-" $log | grep -v "REAL$"
.venv/bin/python bench/live_lo/lo_summary.py $f
