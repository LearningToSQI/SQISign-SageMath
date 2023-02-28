"""
Run with:
    sage -python benchmarks/benchmark_keygen.py
If you want to skip debugging asserts, run
    sage -python -O benchmarks/benchmark_keygen.py

If you get an error about modules not being found
you may have to edit the $PYTHONPATH variable.

I fixed it by adding

export PYTHONPATH=${PWD}

To the file ~/.sage/sagerc
"""

# Python imports
import cProfile
import pstats

# Local imports
from SQISign import SQISign

sqisign = SQISign()

cProfile.run("sqisign.keygen()", "sqisign_keygen.cProfile")
p = pstats.Stats("sqisign_keygen.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(int(40))
