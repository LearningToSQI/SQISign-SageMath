"""
Similar to how `example_SQISign.sage is written, but with more
assertions and debugging prints. 
"""

# Python imports
import time

# Local imports
from SQISign import SQISign
from utilities import print_info
from setup import *

prover = SQISign()
verifier = SQISign()

print_info("Starting SQISign")
sqisign_time = time.time()

# Keygen
print_info("Starting Keygen")
keygen_time = time.time()
prover.keygen()
print_info(f"Keygen took {time.time() - keygen_time:5f}")

# Unpack and check keygen
EA = prover.pk
τ_prime, Iτ, Jτ = prover.sk
assert τ_prime.degree() == Jτ.norm()
assert gcd(Dc, T_prime) == 1

# Commitment
print_info("Starting Commitment")
commitment_time = time.time()
E1 = prover.commitment()
print_info(f"Commitment took {time.time() - commitment_time:5f}")

# Check commitment secret values
ψ_ker, ψ, Iψ = prover.commitment_secrets
assert Iψ.norm() == T_prime
assert ψ_ker.order() == T_prime
assert ψ.degree() == T_prime

# Challenge
print_info("Starting Challenge")
challenge_time = time.time()
ϕ_ker = verifier.challenge(E1)
print_info(f"Challenge took {time.time() - challenge_time:5f}")

# Check Challenge
assert ϕ_ker.order() == Dc

# Response
print_info("Starting Response")
response_time = time.time()
S = prover.response(ϕ_ker)
print_info(f"Response took {time.time() - response_time:5f}")

# Verification
print_info("Starting Verification")
verify_time = time.time()
EA = prover.export_public_key()
response_valid = verifier.verify_response(EA, E1, S, ϕ_ker)

# Check verification
print_info(f"Verification took {time.time() - verify_time:5f}")
assert response_valid, "SQISign response was not valid"
print(f"INFO [SQISign]: SQISign was successful!")

# All finished!
print_info(f"SQISign took {time.time() - sqisign_time:5f}")
