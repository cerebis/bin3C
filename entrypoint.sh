#!/bin/bash --login
# The --login ensures the bash configuration is loaded,

# Temporarily disable strict mode and activate conda:
set +euo pipefail
source /venv/bin/activate
# enable strict mode:
set -euo pipefail

# exec the final command:
cmd="$1"; shift
exec $cmd "$@"

