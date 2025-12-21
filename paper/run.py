#!/usr/bin/python3

import argparse
import os
import subprocess
import sys

skipped = ['NFkB', 'Covid2', 'Akt', 'TumorPillis', 'TumorHu', 'LeukaemiaLeon2021', 'MAPK_5out_bis', 'cLV2', 'QWWC']

DIRS    = {
    '1' : [
        'bench1/f4-direct', 'bench1/f4-direct', 'bench1/f4-direct', 'bench1/f4-direct', 'bench1/f4-direct'
    ],
    '3' : [
    'bench3/results'
    ]
}
ext_to_software = {
    '.jl': 'julia', 
    '.sing': 'Singular --cpus=1 --ticks-per-sec=1000',
    '.mpl': '/opt/maple2025/bin/maple'
}

def pattern_occurs(pattern, string):
    conjunctions = pattern.split("&")
    for c in conjunctions:
        disjunctions = c.split("|")
        holds = False
        for d in disjunctions:
            d = d.strip()
            if d in string:
                holds = True
        if not holds:
            return False
    return True

def main(args):
    all_files = []
    for directory_path in DIRS[args.bench]:
        if not os.path.exists(directory_path):
            continue
        for root, dirs, files in os.walk(directory_path):
            if pattern_occurs(".ipynb", root):
                continue
            for file in files:
                if not os.path.splitext(file)[1] in ext_to_software.keys():
                    continue
                full_path = os.path.join(root, file)
                if not pattern_occurs(args.pattern, full_path):
                    continue
                skip = False
                for name in skipped:
                    if pattern_occurs(name, full_path):
                        print(f'Skipping {full_path}')
                        skip = True
                if skip:
                    continue
                all_files.append(full_path)
    
    print("Will be running these files: ", all_files)
    
    for i, file_path in enumerate(all_files):
        file, ext = os.path.splitext(os.path.basename(file_path))
        log_name = file + "_output.txt"
        log_path = os.path.join(os.path.dirname(file_path), log_name)
        software = ext_to_software[ext]
        print(f"({i+1:2}/{len(all_files)}) {file_path}: Running")
        with open(log_path, 'w') as log_io:
            try:
                result = subprocess.run(
                    [f"ulimit -v {args.memory * 2**20} && timeout -s KILL {args.timeout}s {software} {file_path}"],
                    timeout=args.timeout,
                    stderr=log_io,
                    stdout=log_io,
                    shell=True,
                    text=True,
                )
            except subprocess.TimeoutExpired:
                print(f"({i+1:2}/{len(all_files)}) {file_path}: Process timed out after {args.timeout} seconds.")
            except Exception as e:
                print(f"({i+1:2}/{len(all_files)}) {file_path}: An error occurred: {e}")
            finally:
                log_io.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-b", "--bench", default='3', 
                    help='Benchmark to run.')
    parser.add_argument("-p", "--pattern", default="", 
                    help='Optional pattern to filter.')
    parser.add_argument("-t", "--timeout", default="60", 
                    help='Timeout, in seconds.')
    parser.add_argument("-m", "--memory", default="20", 
                    help='Memory limit, in GB.')
    args = parser.parse_args()
    
    args.memory = int(args.memory)
    args.timeout = int(args.timeout)
    
    main(args)