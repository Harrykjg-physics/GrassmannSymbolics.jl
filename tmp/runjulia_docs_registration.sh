#!/bin/sh
set -eu
julia --project=. tmp/docs_registration_validate.jl > validation.out 2>&1
echo "JULIA_EXIT:$?" > validation.exit
