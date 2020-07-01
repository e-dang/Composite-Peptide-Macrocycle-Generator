#!/bin/sh

export SRC_DIR=$(cd $(dirname "${BASH_SOURCE}"); pwd)

exec python $(dirname "${SRC_DIR}")/scripts/cpmg_entry_point.py "$@"