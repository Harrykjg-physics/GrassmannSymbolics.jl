#!/usr/bin/env bash
set -u
cd /home/jgkong/project/2026/0707/njl
./run_with_memory_guard.sh --limit-gb 20 --log diagnose.log -- bash -lc '. /home/jgkong/bin/run.sh /home/jgkong/project/2026/0707/njl/tmp/diagnose_njl.jl; wait'
tail -n 160 diagnose.log
