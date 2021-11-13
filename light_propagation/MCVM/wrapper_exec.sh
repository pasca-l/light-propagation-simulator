#!/bin/bash

python mkdir_results.py
gcc -lm montecarlo.c
./a.out
