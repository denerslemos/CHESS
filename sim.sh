#!/bin/sh
#!/usr/bin/env python
#SBATCH --time=20:00:00
#SBATCH --ntasks=1
#SBATCH --error=error.dat

export LARGE_FILES="true"

infile=$1

./run.py -o "$infile"

