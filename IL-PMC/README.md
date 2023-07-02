# Potential of monolayer charge

## Publication

These scripts generate inputs, run calculations, and analyze results described in the "Potential of Monolayer Charge" preprint: https://arxiv.org/abs/2301.13681.

## Cloning from this repository

git clone -n --depth=1 --filter=tree:0 https://github.com/vilab-tartu/NaRIBaS

cd NaRIBaS

git sparse-checkout set --no-cone IL-PMC

git checkout

## Results and Analysis

Actual results are available on demand from vladislav.ivanistsev@ut.ee.

Publication Figures are stored in the Figures folder.

In the Analysis folder, the ASE database (pmc8.db) and Jupyter notebook (analysis8.jpynb) contain all data and scripts needed to reproduce two Figures and Table from https://arxiv.org/abs/2301.13681.

## Software installation

conda update -n base -c conda-forge conda

conda create --name gpaw2023

conda activate gpaw2023

conda install -c conda-forge openmpi

conda install -c anaconda jupyter

conda install pandas

conda install -c conda-forge gpaw

conda install -c conda-forge dftd4-python

conda install -c conda-forge jupyterlab

conda install -c conda-forge chargemol

pip install git+https://github.com/funkymunkycool/Cube-Toolz.git

If needed, tweak GPAW according to https://doublelayer.eu/vilab/2022/03/30/bayesian-error-estimation-for-rpbe/

## Running the calculations

Change the scripts to allow running at a cluster with a queueing system.

cd NaRIBaS

chmod -R 755 *

STEP0: Generate models

./naribas STEP0-model.sh

STEP1: Optimize models

./naribas STEP1-optimization.sh

STEP2: Single point calc

./naribas STEP2-singlepoint.sh

STEP2: cDFT calculation

./naribas STEP2-cdft.sh

STEP3: Calculate properties

./naribas STEP3-properties.sh

STEP4: Add data to database

./naribas STEP4-database.sh
