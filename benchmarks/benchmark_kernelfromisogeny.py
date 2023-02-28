"""
Run with:
    sage -python benchmark_kernelfromisogeny.py
If you want to skip debugging asserts, run 
    sage -python -O benchmark_kernelfromisogeny.py
"""

# Python imports
import cProfile
import pstats
import time

# Sage imports
from sage.all import randint, factor

# local imports
from setup import *
from isogenies import (
    torsion_basis,
    EllipticCurveIsogenyFactored,
    kernel_from_isogeny_prime_power
)

# Benchmark new isogeny computation
for D in [Dc, T_prime]:
    print(f"Degree = {factor(D)}")

    torsion_time = time.time()
    P, Q = torsion_basis(E0, D)
    print(f"Basis took: {time.time() - torsion_time:.5f}")

    x = randint(0, D)
    K = P + Q
    K._order = D

    isogeny_time = time.time()
    ϕ = E0.isogeny(K, algorithm="factored")
    print(f"Old Isogeny took: {time.time() - isogeny_time:.5f}")

    isogeny_time = time.time()
    ϕ_new = EllipticCurveIsogenyFactored(E0, K, order=D)
    print(f"New Isogeny took: {time.time() - isogeny_time:.5f}")

    print()

    assert ϕ.codomain().is_isomorphic(ϕ_new.codomain())

cProfile.run("kernel_from_isogeny_prime_power(ϕ)", "kernel_from_isogeny_prime_power.cProfile")
p = pstats.Stats("kernel_from_isogeny_prime_power.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(int(40))
