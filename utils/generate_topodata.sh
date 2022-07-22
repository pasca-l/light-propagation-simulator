#!/bin/bash

# $1 takes in the name of data folder as string

python3 topography_ssp.py $1
python3 topography_dod.py $1
python3 gaussian_fitting.py $1
