#!/bin/sh
# =======================================================================================
#
# This script extracts outputs from ELM-FATES simulations.

# Rutuja Chitra-Tarak (Tue Aug 3, 2021)
# =======================================================================================

# 6. To find which cases are complete (output/Filter.txt) and which are not (output/Missing.txt), run:

python src/create.filter.py

 
# 7. To extract outputs from ELM ensembles (output/exrtact/elm_daily_outputs.txt), run:

python src/extract.output.py

