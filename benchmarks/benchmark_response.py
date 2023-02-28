"""
Run with:
    sage -python benchmarks/benchmark_response.py
If you want to skip debugging asserts, run
    sage -python -O benchmarks/benchmark_response.py

If you get an error about modules not being found
you may have to edit the $PYTHONPATH variable.

I fixed it by adding

export PYTHONPATH=${PWD}

To the file ~/.sage/sagerc
"""

import cProfile
import pstats

# Local imports
from SQISign import SQISign
from setup import *

sqisign = SQISign()

# Keygen
keypair, keypair_ideals = sqisign.keygen()
Iτ, Jτ = keypair_ideals
EA, τ_prime = keypair

# Commitment
Iψ, ψ_ker = sqisign.commitment()
ϕ_ker = sqisign.challenge(ψ_ker)

# Response
cProfile.run(
    "sqisign.response(keypair, Iτ, Jτ, Iψ, ψ_ker, ϕ_ker)", "sqisign_response.cProfile"
)
p = pstats.Stats("sqisign_response.cProfile")
p.strip_dirs().sort_stats("cumtime").print_stats(int(100))
