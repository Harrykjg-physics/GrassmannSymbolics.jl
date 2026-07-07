#!/usr/bin/env bash
set -u
cd /home/jgkong/project/2026/0707/hubbard
rm -f success.sentinel
./run_with_memory_guard.sh --limit-gb 20 --log run.log -- bash -lc '. /home/jgkong/bin/run.sh /home/jgkong/project/2026/0707/hubbard/tmp/run_hubbard_validation.jl; wait'
status=$?
echo RUN_STATUS:$status
test -f success.sentinel && echo SUCCESS_SENTINEL_PRESENT || echo SUCCESS_SENTINEL_MISSING
tail -n 220 run.log
exit $status
