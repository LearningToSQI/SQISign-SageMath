"""
Given parameter sets well designed for SQISign computes all the
global variables needed for various sub-algorithms.

TODO: making so many things global is ugly, we should be initialising 
      classes and passing things around. Will take a bit of refactoring
      so could be done when we try and include the optimisations for
      E/Fp2 and its twist
"""

# Sage imports
from sage.all import (
    proof,
    GF,
    EllipticCurve,
    QuaternionAlgebra,
    ceil,
    log,
    round,
    var,
    gcd,
)

# Local imports
from parameters import params

# Speedy and still (mostly) correct
proof.all(False)


def construct_fields_and_curves(p):
    """
    Given the SQISign prime, compute the finite
    fields and elliptic curves we need. Additionally
    we make `sqrt_minus_one` ∈ Fp4 to allow us to
    compute the twisting endomorphism.
    """
    # Construct finite fields
    x = var("x")
    Fp2 = GF(p**2, name="z2", modulus=x**2 + 1)
    Fp4 = Fp2.extension(2, name="z4")

    # √-1 in Fp4 is needed for automorphisms
    sqrt_minus_one = min(Fp4(-1).sqrt(all=True))  # Deterministic sqrt

    # Supersingular curve with End(E0) known
    E0 = EllipticCurve(Fp4, [1, 0])
    E0.set_order((p**2 - 1) ** 2, num_checks=0)

    return Fp2, Fp4, E0, sqrt_minus_one


def construct_quaternion_alg_and_order(p, q):
    """
    Compute the quaternion algebra and maximal
    order O0. We (currently) will always have 
    p ≡ 3 mod 4 and q = 1.
    """
    # Quaternion algebra things
    B = QuaternionAlgebra(-q, -p)
    i, j, k = B.gens()

    # TODO: derive ω from function depending on p
    # R[ω] is the distinguished quadratic subring 
    # of the algebra B_{p, ∞}.
    # When p ≡ 3 mod 4, we have ω = i
    assert p % 4 == 3, "Currently only p ≡ 3 mod 4 is supported"
    ω = i

    # Maximal order O, and a corresponding left ideal I
    O0 = B.quaternion_order([1, i, (i + j) / 2, (1 + k) / 2])

    return B, i, j, k, ω, O0


def construct_security_parameters(p):
    """
    Given p, derive bounds for other
    security parameters

    p  = 2^(2*λ)
    Bτ ≃ 2^(λ/2)
    eτ ≃ 2*λ
    """
    λ_security = ceil(log(p, 2) / 2)
    Bτ = 2 ** (λ_security // 2)  # 2^(λ/2)
    eτ = 2 * λ_security

    return λ_security, Bτ, eτ


def construct_heuristics(p, l):
    """
    Various KLPT sub-functions require elements
    to be larger or smaller than bounds. Derive
    them here and pass them through as globals.
    """
    sqrt_p = round(p.sqrt())
    logp = log(p, l)
    loglogp = log(logp, l)

    prime_norm_heuristic = l ** ceil(logp / 2 + loglogp)
    represent_heuristic = ceil(logp / 4) - 4

    return sqrt_p, logp, loglogp, prime_norm_heuristic, represent_heuristic


"""
Given security parameters generate everything we need
"""
# Extract values from dictionary
p, q, l, available_torsion, T, e, Δ, T_prime, Dc, f_step_max = params.values()

# Construct fields and curves
Fp2, Fp4, E0, sqrt_minus_one = construct_fields_and_curves(p)

# Construct Quaternion algebra and maximal order
B, i, j, k, ω, O0 = construct_quaternion_alg_and_order(p, q)

# Construct security parameters
λ_security, Bτ, eτ = construct_security_parameters(p)

# Construct some fixed constants for heuristics
# (Mainly for KLPT)
sqrt_p, logp, loglogp, prime_norm_heuristic, represent_heuristic = construct_heuristics(
    p, l
)

# Check algebra invariants
assert q == 1, "Only q = 1 is supported"
assert p % 4 == 3, "p % 4 ≡ 1 is not currently supported"
assert p % 12 == 7, "p % 12 ̸≡ 7 is not currently supported"

# Check soundness of torsion subgroups
assert available_torsion % l == 0, "l does not divide the available torsion"
assert available_torsion % T == 0, "T does not divide the available torsion"
assert available_torsion % T_prime == 0, "T_prime does not divide the available torsion"
assert available_torsion % Dc == 0, "Dc does not divide the available torsion"
assert (
    available_torsion % 2**f_step_max == 0
), "2^f_step_max does not divide the available torsion"
assert gcd(Dc, T_prime) == 1, "Challenge and commitment torsions are not coprime"

# Check field and curves
assert E0.a_invariants() == (0, 0, 0, 1, 0), "Only E : y^2 = x^3 + x is supported"
assert E0.is_supersingular(), "E0 is not supersingular"
assert sqrt_minus_one**2 == -1, "sqrt_minus_one is not the square root of minus one!"

# Check Quaternion Algebra and Order
assert O0.discriminant() == B.discriminant(), "O0 is not a maximal order"

# Check security parameters
assert Dc > l**λ_security, "Dc is not large enough"
