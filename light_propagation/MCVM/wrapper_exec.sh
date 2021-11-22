#!/bin/bash

python3 mkdir_data.py
gcc montecarlo.c -lm -std=c99
./a.out
python3 push_data.py $1
