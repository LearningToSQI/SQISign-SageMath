"""
Run with:
    sage -python benchmarks/benchmark_idealtoisogeny.py
If you want to skip debugging asserts, run
    sage -python -O benchmarks/benchmark_idealtoisogeny.py

If you get an error about modules not being found
you may have to edit the $PYTHONPATH variable.

I fixed it by adding

export PYTHONPATH=${PWD}

To the file ~/.sage/sagerc
"""

# Python imports
import cProfile
import pstats

# SageMath imports
from sage.all import ZZ

# Local imports
from deuring import IdealToIsogenyFromKLPT
from KLPT import RepresentIntegerHeuristic
from utilities import inert_prime
from setup import E0, O0, Bτ, eτ, p, l, ω


# =============================#
# This is just SQISign keygen #
# =============================#

Nl = l**eτ
while True:
    Nτ = inert_prime(Bτ, -ZZ(ω**2))
    # We need the product to be large enough for
    # RepresentIntegerHeuristic.
    if Nτ * Nl > 2 * p:
        break

# Compute an endomorphism γ of norm Nτ l^eτ
# Nτ < Bτ =
γ = None
g = 0
while γ is None:
    γ = RepresentIntegerHeuristic(Nτ * Nl, parity=True)
Iτ = O0 * γ + O0 * Nτ
Jτ = O0 * γ.conjugate() + O0 * Nl

# Compute the secret isogeny τ
I_trivial = O0.unit_ideal()
ϕ_trivial = E0.isogeny(E0(0))

cProfile.run(
    "IdealToIsogenyFromKLPT(Jτ, I_trivial, ϕ_trivial, I_prime=Iτ)",
    "IdealToIsogenyFromKLPT.cProfile",
)
p = pstats.Stats("IdealToIsogenyFromKLPT.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(int(40))
