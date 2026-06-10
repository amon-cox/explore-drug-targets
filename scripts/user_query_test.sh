#!/usr/bin/env bash

# Impersonated user: exploring drug candidates for SCN5A (a cardiac sodium channel)
# Queries by HUGO symbol and Ensembl ID, with optiosn for output.

# default settings, target required
echo "=== Query 1: Top 10 compounds for SCN5A (default settings) ==="
Rscript R/drug_target_cli.R --target SCN5A
echo ""

# same target with Ensembl ID, adjusted filters, custom file name, and plot
echo "=== Query 2: Top 5 compounds for ENSG00000183873 (Phase 4) ==="
Rscript R/drug_target_cli.R --target ENSG00000183873 --top 5 --minPhase 4 --output ENSG00000183873_results --plot
echo ""
