"""
Run with:
    sage -python benchmarks/benchmark_BiDLP.py
If you want to skip debugging asserts, run
    sage -python -O benchmarks/benchmark_BiDLP.py

If you get an error about modules not being found
you may have to edit the $PYTHONPATH variable.

I fixed it by adding

export PYTHONPATH=${PWD}

To the file ~/.sage/sagerc
"""

# Python imports
import cProfile
import pstats

# Sagemath imports
from sage.all import randint

# Local imports
from setup import *
from isogenies import torsion_basis, BiDLP

P, Q = torsion_basis(E0, T_prime)

cProfile.run("torsion_basis(E0, T_prime)", "torsion_basis.cProfile")
p = pstats.Stats("torsion_basis.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(int(40))

x = randint(0, T_prime)
R = x * P + Q

cProfile.run("BiDLP(R, P, Q, T_prime)", "BiDLP.cProfile")
p = pstats.Stats("BiDLP.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(int(40))
