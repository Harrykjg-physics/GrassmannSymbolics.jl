#!/usr/bin/env bash
set -u
cd /home/jgkong/project/2026/0707/hubbard
./run_with_memory_guard.sh --limit-gb 20 --log measure.log -- bash -lc '. /home/jgkong/bin/run.sh /home/jgkong/project/2026/0707/hubbard/tmp/diagnose_hubbard_measures.jl; wait'
tail -n 260 measure.log
