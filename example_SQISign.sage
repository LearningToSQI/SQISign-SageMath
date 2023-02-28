"""
Example of SQISign as a one-round interactive identification
protocol.

We imagine two parties, `prover` and `verifier`. The `prover`
demonstrates knowledge of the endomorphism ring End(EA) in the 
following way:

- The prover's public key is an elliptic curve EA, and their secret
is the isogeny E0 → EA, where End(E0) is known to everyone.

- The prover then makes a commitment by computing a second secret 
isogeny ψ : E0 → E1 and sends the codomain to the verifier.

- The verifier makes a challenge ϕ: E1 → E2 and sends ϕ to the prover

- The prover responds with σ : EA → E2, which is done via knowledge 
of End(EA) (through knowing End(E0) and τ : E0 → EA). The prover
sends σ to the verifier. 

- If ϕ_dual ∘ σ : EA → E1 is a cyclic isogeny, the verifier returns
true and false otherwise
"""

# Python imports
import time

# Local imports
from SQISign import SQISign
from utilities import print_info

# SQISign is a protocol between a prover and verifier
prover = SQISign()
verifier = SQISign()

print_info("Starting SQISign")
sqisign_time = time.time()

# The prover generates their keypair and makes a commitment
# which is a secret isogeny ψ : E0 → E1 and sends the codomain
# of ψ to the verifier
print_info("Computing Keypair")
prover.keygen()
EA = prover.export_public_key()

print_info("Computing Commitment")
E1 = prover.commitment()

# The verifier receives the commitment and makes a random
# challenge. Rather than sending the isogeny, the verifier
# simply sends the generator of the isogeny ϕ : E1 → E2
print_info("Computing Challenge")
phi_ker = verifier.challenge(E1)

# The verifier makes a response to the challenge, which is
# an isogeny σ : EA → E2 of degree l^e. The prover compresses
# σ to a bitstring S and sends this to the verifier
print_info("Computing Response")
S = prover.response(phi_ker)

# The verifier uses the prover's public key and their response
# S and checks if σ is an isogeny EA → E2 of degree l^e and
# whether ϕ_dual ∘ σ : EA → E1 is cyclic
print_info("Validating response")
valid = verifier.verify_response(EA, E1, S, phi_ker)

print_info(f"SQISign example worked: {valid}")
print_info(f"SQISign took {time.time() - sqisign_time:5f}")
