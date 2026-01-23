#!/usr/bin/python3

import argparse
import os
import subprocess
import sys

def benchmark_1():
    model = "SLIQR"
    print("Generating benchmarks for different software")
    try:
        subprocess.run([f'julia benchmark-1/populate.jl {model}'], shell=True, text=True)
    except:
        print("Beda")
    
    print("Running a simple benchmark for every software")
    try:
        subprocess.run([f'python src/run.py -p "benchmark-1 & {model}" -t 3600 -m 10'], shell=True, text=True)
    except:
        print("Beda")
    
    for method in ['f4-direct', 'f4-flat', 'ffmodstd', 'paramgb', 'slimgb']:
        os.path.isfile(f"benchmark-1/{method}/{model}_output.txt")
        with open(f"benchmark-1/{method}/{model}_output.txt") as io:
            content = io.readlines()
        print(f"{method}: OK")

def benchmark_2():
    model = "SLIQR"
    try:
        subprocess.run([f'julia benchmark-2/main.jl {model}'], shell=True, text=True)
    except:
        print("Beda")

def benchmark_3():
    model = "SLIQR"
    try:
        subprocess.run([f'julia benchmark-3/populate.jl {model}'], shell=True, text=True)
    except:
        print("Beda")
    try:
        subprocess.run([f'python src/run.py -p "benchmark-3 & simplify & {model}.jl" -t 3600 -m 10'], shell=True, text=True)
    except:
        print("Beda")
    try:
        subprocess.run([f'python src/run.py -p "benchmark-3 & maple_simplify & {model} & generate_gens" -t 3600 -m 10'], shell=True, text=True)
        subprocess.run([f'python src/run.py -p "benchmark-3 & maple_simplify & run_simplify & {model} & .mpl" -t 3600 -m 10'], shell=True, text=True)
    except:
        print("Beda")
    try:
        subprocess.run([f'python src/run.py -p "benchmark-3 & input_stats & {model}" -t 3600 -m 10'], shell=True, text=True)
    except:
        print("Beda")
    try:
        subprocess.run([f'python src/run.py -p "benchmark-3 & independence & {model}" -t 3600 -m 10'], shell=True, text=True)
    except:
        print("Beda")
    try:
        subprocess.run([f'python benchmark-3/collect.jl'], shell=True, text=True)
    except:
        print("Beda")
        
# benchmark_1()
# benchmark_2()
benchmark_3()
