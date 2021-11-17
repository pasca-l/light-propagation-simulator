#!/bin/bash

python3 mkdir_data.py
gcc -lm -std=c99 montecarlo.c
./a.out
python3 push_data.py $1
