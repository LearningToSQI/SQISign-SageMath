"""
Run with:
    sage -python benchmarks/benchmark_torsionbasis.py
If you want to skip debugging asserts, run
    sage -python -O benchmarks/benchmark_torsionbasis.py

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

# Local imports
from isogenies import generate_point_order_D, torsion_basis
from setup import *


t_new = 0
iteration = 50
for _ in range(iteration):
    tmp = time.time()
    generate_point_order_D(E0, T_prime)
    t_new += time.time() - tmp

t_new = t_new / iteration

print(f"New time generate_point_order_D: {t_new:.5f}")

t_new = 0
iteration = 50
for _ in range(iteration):
    tmp = time.time()
    torsion_basis(E0, T_prime)
    t_new += time.time() - tmp

t_new = t_new / iteration

print(f"New time torsion_basis: {t_new:.5f}")


# cProfile.run("generate_point_order_D_old(E0, T_prime)", 'generate_point_order_D_old.cProfile')
# p = pstats.Stats('generate_point_order_D_old.cProfile')
# p.strip_dirs().sort_stats("cumtime").print_stats(int(10))

# cProfile.run("generate_point_order_D(E0, T_prime)", 'generate_point_order_D.cProfile')
# p = pstats.Stats('generate_point_order_D.cProfile')
# p.strip_dirs().sort_stats("cumtime").print_stats(int(10))
