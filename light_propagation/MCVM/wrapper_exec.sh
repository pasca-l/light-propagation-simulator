#!/bin/bash

python mkdir_data.py
gcc -lm montecarlo.c
./a.out
python push_data.py $1
