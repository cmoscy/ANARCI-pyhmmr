#!/usr/bin/env bash
# Copy pipeline outputs into the packaged source tree (lib/python/anarci).
# Run from the repository root after build_pipeline/RUN_pipeline.sh succeeds.

set -euo pipefail

ROOT="$(cd "$(dirname "$0")/.." && pwd)"
PIPELINE="${ROOT}/build_pipeline"
PKG="${ROOT}/lib/python/anarci"
HMM_DEST="${PKG}/dat/HMMs"

GERMLINES_SRC="${PIPELINE}/curated_alignments/germlines.py"
HMM_SRC="${PIPELINE}/HMMs/ALL.hmm"

if [[ ! -f "${GERMLINES_SRC}" ]]; then
  echo "Error: ${GERMLINES_SRC} not found. Run build_pipeline/RUN_pipeline.sh first." >&2
  exit 1
fi

if [[ ! -f "${HMM_SRC}" ]]; then
  echo "Error: ${HMM_SRC} not found. Run build_pipeline/RUN_pipeline.sh first." >&2
  exit 1
fi

mkdir -p "${HMM_DEST}"
cp "${GERMLINES_SRC}" "${PKG}/germlines.py"
cp -f "${PIPELINE}"/HMMs/ALL.hmm* "${HMM_DEST}/"

echo "Synced packaged data:"
echo "  ${PKG}/germlines.py"
echo "  ${HMM_DEST}/"
echo ""
echo "Next: commit the updated files and smoke-test, e.g.:"
echo "  python -c \"from anarci import anarci; print(anarci([('t','EVQLQQSGAEVVRSGASVKLSCTASGFNIKDYYIHWVKQRPEKGLEWIGWIDPEIGDTEYVPKFQGKATMTADTSSNTAYLQLSSLTSEDTAVYYCNAGHDYDRGRFPYWGQGTLVTVSA')], scheme='imgt', output=False))\""
