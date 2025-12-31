#!/usr/bin/python3

import argparse
import os
import subprocess
import sys

skipped = ['NFkB', 'Covid2', 'Akt', 'TumorPillis', 'TumorHu', 'LeukaemiaLeon2021', 'MAPK_5out_bis', 'cLV2', 'QWWC']

ext_to_software = {
    '.jl': 'julia', 
    '.sing': 'Singular --cpus=1 --ticks-per-sec=1000',
    '.mpl': '/opt/maple2025/bin/maple'
}

def pattern_occurs(pattern, string):
    conj = pattern.split("&")
    for c in conj:
        disj = c.split("|")
        holds = False
        for d in disj:
            d = d.strip()
            if d in string:
                holds = True
        if not holds:
            return False
    return True

def walk(args, path):
    if not os.path.exists(path):
        return []
    if pattern_occurs(".ipynb", path):
        return []
    if os.path.isfile(path):
        filename = os.path.basename(path)
        if not os.path.splitext(filename)[1] in ext_to_software.keys():
            return []
        if not pattern_occurs(args.pattern, path):
            return []
        for name in skipped:
            if pattern_occurs(name, path):
                print(f'Skipping {path}')
                return []
        return [path]
    else:
        experience = []
        for fork in os.listdir(path):
            experience += walk(args, os.path.join(path, fork))
        return experience
    
def main(args):
    all_files = walk(args, os.getcwd())
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
    parser.add_argument("-p", "--pattern", default="", 
                    help='Pattern to filter. The format is CNF, the syntax is & and |.')
    parser.add_argument("-t", "--timeout", default="60", 
                    help='Timeout, in seconds.')
    parser.add_argument("-m", "--memory", default="20", 
                    help='Memory limit, in GB.')
    args = parser.parse_args()
    
    args.memory = int(args.memory)
    args.timeout = int(args.timeout)
    
    main(args)