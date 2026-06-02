# Maintainer pipeline: regenerate packaged HMMs and germlines

This directory rebuilds `germlines.py` and the HMM database from IMGT germline data. **End-user installs do not run this pipeline** — they use the prebuilt files under `lib/python/anarci/`.

## When to run

- IMGT germline data has been updated and you want to refresh the fork
- You are bootstrapping a checkout that has no packaged data under `lib/python/anarci/dat/HMMs/`
- You need to verify or reproduce how the bundled assets were produced

Do **not** run this as part of `pip install` or `uv sync`.

## Requirements

- Network access (downloads from IMGT via `RipIMGT.py`)
- [HMMER3](http://hmmer.org/) (`hmmbuild`, `hmmpress`) on `PATH`
- MUSCLE (used by `FormatAlignments.py`; see upstream `bin/` or system install)
- Python 3 with dependencies needed by the pipeline scripts

## Steps

From the repository root:

```bash
cd build_pipeline
bash RUN_pipeline.sh
cd ..
bash build_pipeline/sync_packaged_data.sh
```

`sync_packaged_data.sh` copies:

- `build_pipeline/curated_alignments/germlines.py` → `lib/python/anarci/germlines.py`
- `build_pipeline/HMMs/ALL.hmm*` → `lib/python/anarci/dat/HMMs/`

Then commit the updated files under `lib/python/anarci/` and tag a release if appropriate.

Intermediate directories (`IMGT_sequence_files`, `muscle_alignments`, `curated_alignments`, `HMMs` under `build_pipeline/`) are gitignored.

## Runtime note (pyhmmer fork)

Installed ANARCI uses **pyhmmer** for alignments. System `hmmscan` is not required for the Python API or the `ANARCI` CLI in this fork. The pipeline still uses HMMER3 command-line tools to **build** the packaged `ALL.hmm` database.
