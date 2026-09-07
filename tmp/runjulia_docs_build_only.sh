#!/bin/sh
set -eu
julia --project=. tmp/docs_build_only_validate.jl > docs_build.out 2>&1
echo "JULIA_EXIT:$?" > docs_build.exit
