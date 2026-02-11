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
    conj = pattern.split(",")
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
        print(f"({i+1:2}/{len(all_files)}) Running: {file_path}")
        print(f"({i+1:2}/{len(all_files)}) Results will be written to: {log_path}")
        ulimit = "" if args.memory == 0 else f"ulimit -v {args.memory * 2**20} && "
        with open(log_path, 'w') as log_io:
            try:
                result = subprocess.run(
                    # [f"{ulimit}timeout -s KILL {args.timeout}s {software} {file_path}"],
                    [f"{ulimit} {software} {file_path}"],
                    timeout=args.timeout,
                    stderr=log_io,
                    stdout=log_io,
                    shell=True,
                    text=True,
                )
            except subprocess.TimeoutExpired:
                print(f"({i+1:2}/{len(all_files)}) Process timed out after {args.timeout} seconds: {file_path}")
            except Exception as e:
                print(f"({i+1:2}/{len(all_files)}) An error occurred: {e} : {file_path}")
            finally:
                log_io.close()

if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("-p", "--pattern", default="", 
                    help='Pattern to filter. Can use & and |.')
    parser.add_argument("-t", "--timeout", default="3600", 
                    help='Timeout, in seconds.')
    parser.add_argument("-m", "--memory", default="0", 
                    help='Memory limit, in GB.')
    parser.add_argument("--maple", default=ext_to_software[".mpl"], help='Command to run Maple (e.g., /opt/maple2025/bin/maple).')
    parser.add_argument("--singular", default="Singular --cpus=1 --ticks-per-sec=1000", help=ext_to_software[".sing"])
    args = parser.parse_args()
    
    args.memory = int(args.memory)
    args.timeout = int(args.timeout)

    ext_to_software['.sing'] = args.singular
    ext_to_software['.mpl'] = args.maple
    
    main(args)