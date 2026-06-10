# explore-drug-targets

A tool to explore OpenTargets data to rank compounds by selectivity and adverse effect profile for a given target.

1. filters `molecules.tsv.gz` to clinical/approved compounds,
2. counts the number of targets per compound (selectivity),
3. aggregates `adverseEffects.tsv.gz` LLR statistics per compound, and
4. ranks compounds by increasing number of targets and increasing maximum LLR.

## Setup

Clone the repository and restore the R environment:

```bash
git clone <repo-url>
cd explore-drug-targets
R -e "install.packages('renv'); renv::restore()"
```

## Usage

```bash
Rscript R/drug_target_cli.R --target <gene> [options]
```

| Option | Type | Default | Description |
| --- | --- | --- | --- |
| `--target` | string | required | Ensembl gene ID (e.g. `ENSG00000183873`) or HUGO symbol (e.g. `SCN5A`) |
| `--minPhase` | integer | 3 | Minimum clinical trial phase; approved compounds always included |
| `--top` | integer | 10 | Number of top-ranked compounds to print |
| `--output` | string | auto | Output filename (no extension); defaults to a timestamped name |
| `--plot` | flag | FALSE | Save a scatterplot of all ranked compounds |

## Examples

```bash
# example run
Rscript R/drug_target_cli.R \
  --target ENSG00000146648 \
  --minPhase 4 \
  --top 10 \
  --output results_egfr \
  --plot

# impersonated user
scripts/user_query_test.sh
```

## Project strutcure

- `data/` # input TSV files
- `R/` # CLI script
- `scripts/` # example user workflow
- `results/` # output TSVs and plots (not tracked)
