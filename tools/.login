#!/bin/tcsh

alias ls 'ls --color'

## 1. set environmental variables
setenv USERDIR /usr/projects/higrad/$user
setenv SCRATCHDIR /usr/projects/higrad/$user

setenv CASE_ROOT    $USERDIR/E3SM_cases

setenv E3SM_ROOT    $SCRATCHDIR/E3SM.DRM/E3SM
setenv ELM_ROOT     $SCRATCHDIR/E3SM.DRM/E3SM/components/elm/src
setenv FATES_ROOT   $SCRATCHDIR/E3SM.DRM/E3SM/components/elm/src/external_models/fates
setenv RUN_ROOT     $SCRATCHDIR/E3SM.DRM/runs
setenv ARCHIVE_ROOT $SCRATCHDIR/E3SM.DRM/archive
setenv DIN_LOC_ROOTF /lustre/scratch4/turquoise/$user/Data/Data.E3SM/inputdata # root directory of all CIME and component input data

## 2. Load modules.
setenv CXX "CC"
setenv CC "cc"
module load gcc/11.2.0
module swap PrgEnv-gnu/8.4.0 PrgEnv-gnu/8.3.3
module load cmake/3.20.3 friendly-testing mkl cray-mpich/8.1.26 cray-hdf5-parallel/1.12.2.1 cray-netcdf-hdf5parallel/4.9.0.1 python/3.9-anaconda-2021.11

