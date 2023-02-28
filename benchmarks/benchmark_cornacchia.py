"""
Run with:
    sage -python benchmarks/benchmark_cornacchia.py
If you want to skip debugging asserts, run
    sage -python -O benchmarks/benchmark_cornacchia.py

If you get an error about modules not being found
you may have to edit the $PYTHONPATH variable.

I fixed it by adding

export PYTHONPATH=${PWD}

To the file ~/.sage/sagerc
"""

# Python imports
import cProfile
import pstats
import time

# Sagemath imports
from sage.all import randint

# Local imports
from utilities import Cornacchia

total_time = 0
iter_count = 10000
for _ in range(iter_count):
    x = randint(0, 2**256)
    t0 = time.time()
    sol = Cornacchia(x, 1)
    total_time += (time.time() - t0)

print(f"Average time: {total_time / iter_count:5f}")
