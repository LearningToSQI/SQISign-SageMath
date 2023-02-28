"""
Run with:
    sage -python benchmarks/benchmark_sqisign.py
If you want to skip debugging asserts, run
    sage -python -O benchmarks/benchmark_sqisign.py

If you get an error about modules not being found
you may have to edit the $PYTHONPATH variable.

I fixed it by adding

export PYTHONPATH=${PWD}

To the file ~/.sage/sagerc
"""

import cProfile
import pstats
import time


from SQISign import SQISign


def print_info(str):
    print("="*80)
    print(f"{str}".center(80))
    print("="*80)


# Start the profiler
pr = cProfile.Profile()
pr.enable()

#
#
#

prover = SQISign()
verifier = SQISign()

print_info("Starting SQISign")
sqisign_time = time.time()

# Keygen
print_info("Starting Keygen")
keygen_time = time.time()
prover.keygen()
EA = prover.export_public_key()
print_info(f"Keygen took {time.time() - keygen_time:5f}")

# Commitment
print_info("Starting Commitment")
commitment_time = time.time()
E1 = prover.commitment()
print_info(f"Commitment took {time.time() - commitment_time:5f}")

# Challenge
print_info("Starting Challenge")
challenge_time = time.time()
ϕ_ker = verifier.challenge(E1)
print_info(f"Challenge took {time.time() - challenge_time:5f}")

# Response
print_info("Starting Response")
response_time = time.time()
S = prover.response(ϕ_ker)
print_info(f"Response took {time.time() - response_time:5f}")

# Verification
print_info("Starting Verification")
verify_time = time.time()
response_valid = verifier.verify_response(EA, E1, S, ϕ_ker)

# Check verification
print_info(f"Verification took {time.time() - verify_time:5f}")
assert response_valid, "SQISign response was not valid"
print(f"INFO [SQISign]: SQISign was successful!")

# All finished!
print_info(f"SQISign took {time.time() - sqisign_time:5f}")

pr.disable()
pr.dump_stats("sqisign.cProfile")
p = pstats.Stats('sqisign.cProfile')
p.strip_dirs().sort_stats("cumtime").print_stats(250)
