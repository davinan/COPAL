#!/usr/bin/env bash
# Source this before copal predict.
REPO=/mnt/home/dnan/projects/luxr/src/publication_repo
export AF3_CONTAINER=/software/containers/af3.sif
export AF3_DB_PATH=/databases/lab/af3_DB
export AF3_ROOT=$REPO/external/alphafold3
export AF3_MODEL_DIR=/databases/alphafold
export BOLTZ2_BIN=$REPO/external/boltz/.venv/bin/boltz
export RF3_INSTALL=$REPO/external/foundry
export RF3_CKPT=$RF3_INSTALL/rf3_latest.pt
export COPAL_AF3_SCRIPT=$REPO/external/wrappers/af3_native.sh
export COPAL_BOLTZ2_SCRIPT=$REPO/external/wrappers/boltz2.sh
export PATH=/software/hmmer/latest/bin:$PATH
