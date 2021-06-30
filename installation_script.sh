#!/bin/bash
conda install -c conda-forge mamba

mamba create -f install/condaenv.yaml

conda activate 4dh

Rscript install/install_R.R

