# Contributing

Thanks for your interest in this pipeline.

## Before you open a pull request

Please open an issue first describing what you'd like to change. These scripts
correspond to a specific published analysis, so changes that alter results need
to be discussed before they're merged.

## Ground rules for this repository

**Never commit absolute paths.** Every site-specific location belongs in
`config.sh`, which is git-ignored. If your change needs a new path, add it to
`config.sh.example` with a comment explaining what the directory must contain,
and reference it as `${YOUR_VAR}` in the script.

**No data files.** This repository holds code only. Sequence data, VCFs,
genotypes, and phenotypes are hosted externally and linked from the README.
This matters especially for any dataset under a data use agreement.

**Keep sampling reproducible.** Anything that draws a random subset must accept
a seed and default to the one used in the published run.

## Reporting a problem

Issues are welcome, particularly for missing steps or commands that fail on a
scheduler other than Slurm. Please include the script name, the command you ran,
and the error output.
