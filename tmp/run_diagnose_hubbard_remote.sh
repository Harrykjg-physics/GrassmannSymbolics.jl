#!/usr/bin/env bash
set -u
cd /home/jgkong/project/2026/0707/hubbard
./run_with_memory_guard.sh --limit-gb 20 --log diagnose.log -- bash -lc '. /home/jgkong/bin/run.sh /home/jgkong/project/2026/0707/hubbard/tmp/diagnose_hubbard.jl; wait'
tail -n 180 diagnose.log
