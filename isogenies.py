"""
Helper functions for isogeny computations. Mainly, functions to
compute torsion basis, isogenies, kernels generators of kernels
and dual isogenies.

Mainly `EllipticCurveIsogenyFactored` is a wrapper around the usual
factored isogeny in SageMath but also allows a call to the velusqrt
algorithm when the prime isogeny is large. In practice, we seem to 
need l > 2000 to see a benefit (probably because we work in Fp4?) 

We also include the sparse strat for isogenies of prime powers, this
should be in Sage by default soon (TM).

The fix in ticket #34732 should help with Velu performance, so need
to retest once working with 9.8.

TODO:

Throughout this code we assume that we are working on E / Fp^4. 
When we move back to Fp2, we will need to modify this code
"""

# Sage Imports
from sage.all import (
    ZZ,
    cached_function,
    EllipticCurveIsogeny,
    factor,
    set_random_seed,
)
from sage.schemes.elliptic_curves.hom_velusqrt import EllipticCurveHom_velusqrt
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

# Local imports
from utilities import has_order_D
from pari_interface import discrete_log_pari

# =========================================== #
# Compute points of order D and Torsion Bases #
# =========================================== #

def generate_random_point(E, seed=None):
    """
    E0.random_point() returns either P
    or -P with equal probability.

    We always select the element with
    smaller y-coordinate to make this
    deterministic.
    """
    # Allow a setting of the seed to
    # ensure the same point is always returned
    if seed is not None:
        set_random_seed(seed)

    P = E.random_element()
    return min(P, -P)


def generate_point_order_D(E, D):
    """
    Input:  An elliptic curve E / Fp4
            An integer D dividing (p^2 - 1)
    Output: A point P of order D.
    """
    p = E.base().characteristic()
    n = (p**2 - 1) // D
    for _ in range(1000):
        P = n * generate_random_point(E)

        # Case when we randomly picked
        # a point in the n-torsion
        if P.is_zero():
            continue

        # Check that P has order exactly D
        if has_order_D(P, D):
            P._order = ZZ(D)
            return P

    raise ValueError(f"Never found a point of order D...")


def generate_linearly_independent_point(E, P, D, canonical=False):
    """
    Input:  An elliptic curve E / Fp4
            A point P ∈ E[D]
            An integer D dividing (p^2 - 1)
    Output: A point Q such that E[D] = <P, Q>
    """
    # This ensures all random points are deterministically
    # generated
    if canonical:
        set_random_seed(0)

    for _ in range(2000):
        # Generate a random point of order D
        Q = generate_point_order_D(E, D)

        # Make sure the point is linearly independent
        pair = P.weil_pairing(Q, D, algorithm="pari")
        if has_order_D(pair, D, multiplicative=True):
            Q._order = ZZ(D)
            return Q

    raise ValueError(f"Never found a linearly independent point...")


@cached_function
def torsion_basis(E, D, canonical=False):
    """
    Generate basis of E(Fp^4)[D] of supersingular curve

    Optional   canonical: bool
               Makes torsion generation deterministic.
    """
    p = E.base().characteristic()

    # Ensure D divides the curve's order
    if (p**2 - 1) % D != 0:
        print(f"{factor(D) = }")
        print(f"{factor(p**2 - 1) = }")
        raise ValueError(f"D must divide the point's order")

    # This ensures all random points are deterministically
    # generated
    if canonical:
        set_random_seed(0)

    P = generate_point_order_D(E, D)
    Q = generate_linearly_independent_point(E, P, D)

    return P, Q

# ================================== #
# Compute composite degree isogenies #
# ================================== #

def EllipticCurveIsogenyFactored(E, P, order=None, velu_bound=400):
    """
    Works similarly to EllipticCurveHom_composite
    but with two main additions:

    Introduces a sparse strategy for prime power
    isogenies, taken from
    https://trac.sagemath.org/ticket/34239
    This should be default soon (9.8 maybe)

    For primes l > 400, we use velusqrt as
    the algorithm. This bound was found by testing
    in tests/test_isogenies.sage

    Additionally, we allow `order` as an optional parameter
    and `velu_bound` controls when sqrtvelu kicks in
    """

    def EllipticCurveHom_velusqrt_setorder(P):
        """
        To speed things up, we manually set the order
        assuming all curves have order (p^2 - 1)^2

        I think this is fixed for 9.8, but not everyone
        will be running the latest SageMath version.
        """
        E = P.curve()
        p = E.base().characteristic()
        E._order = ZZ((p**2 - 1) ** 2)
        return EllipticCurveHom_velusqrt(E, P)

    def evaluate_factored_isogeny(phi_list, P):
        """
        Given a list of isogenies, evaluates the
        point for each isogeny in the list
        """
        for phi in phi_list:
            P = phi(P)
        return P

    def sparse_isogeny_prime_power(P, l, e, split=0.8, velu_bound=2000):
        """
        Compute chain of isogenies quotienting
        out a point P of order l**e
        https://trac.sagemath.org/ticket/34239
        """
        if l > velu_bound:
            isogeny_algorithm = lambda Q, l: EllipticCurveHom_velusqrt_setorder(Q)
        else:
            isogeny_algorithm = lambda Q, l: EllipticCurveIsogeny(
                Q.curve(), Q, degree=l, check=False
            )

        def recursive_sparse_isogeny(Q, k):
            assert k
            if k == 1:  # base case
                return [isogeny_algorithm(Q, l)]

            k1 = int(k * split + 0.5)
            k1 = max(1, min(k - 1, k1))  # clamp to [1, k-1]

            Q1 = l**k1 * Q
            L = recursive_sparse_isogeny(Q1, k - k1)

            Q2 = evaluate_factored_isogeny(L, Q)
            R = recursive_sparse_isogeny(Q2, k1)

            return L + R

        return recursive_sparse_isogeny(P, e)

    # Ensure P is a point on E
    if P.curve() != E:
        raise ValueError(f"The supplied kernel must be a point on the curve E")

    if order:
        P._order = ZZ(order)
    cofactor = P.order()

    # Deal with isomorphisms
    if cofactor == 1:
        return EllipticCurveIsogeny(P.curve(), P)

    ϕ_list = []
    for l, e in cofactor.factor():
        # Compute point Q of order l^e
        D = ZZ(l**e)
        cofactor //= D
        Q = cofactor * P

        # Manually setting the order means 
        # Sage won't try and do it for each
        # l-isogeny in the iteration
        Q._order = D

        # Use Q as kernel of degree l^e isogeny
        ψ_list = sparse_isogeny_prime_power(Q, l, e, velu_bound=velu_bound)

        # Map P through chain length e of l-isogenies
        P = evaluate_factored_isogeny(ψ_list, P)
        ϕ_list += ψ_list

    return EllipticCurveHom_composite.from_factors(ϕ_list)


# ===================================== #
# Fast computation of dual given kernel #
# ===================================== #

def dual_isogeny_and_kernel(ϕ, R, order=None):
    """
    Compute the dual isogeny given an
    isogeny and its kernel.
    Inspired by ia.cr/2019/499

    Input: An isogeny ϕ : E -> E / <R> of degree D
           The generator R of the ker(ϕ)

    Output: The dual of ϕ and its kernel
    """
    E1 = ϕ.domain()
    E2 = ϕ.codomain()

    if not order:
        D = ZZ(ϕ.degree())
    else:
        D = ZZ(order)

    S = generate_linearly_independent_point(E1, R, D)

    ker_phi_dual = ϕ(S)
    ϕ_dual = EllipticCurveIsogenyFactored(E2, ker_phi_dual, order=D)
    Eϕ = ϕ_dual.codomain()

    # Correcting isomorphism
    if Eϕ != E1:
        iso = Eϕ.isomorphism_to(E1)
        ϕ_dual = iso * ϕ_dual

    return ϕ_dual, ker_phi_dual


def dual_isogeny(ϕ, R, order=None):
    """
    Wrapper function for when we only want the
    dual isogeny but not the dual isogeny's
    kernel
    """
    ϕ_dual, _ = dual_isogeny_and_kernel(ϕ, R, order=order)
    return ϕ_dual


# ===================================== #
# Compute kernel generators of degree D #
# ===================================== #

def generate_kernels_prime_power(E, l, e):
    """
    Compute a list of points of order
    exactly l**e which should generate
    all isogenies of degree l**e
    """
    D = l**e
    P, Q = torsion_basis(E, D)

    # Send points of the form
    # [x]P + Q
    K = Q
    yield K
    for _ in range(D - 1):
        K += P
        K._order = ZZ(D)
        yield K

    # Send points of the form
    # P + [lx]Q
    K = P
    yield K
    lQ = l * Q
    for _ in range(D // l - 1):
        K += lQ
        K._order = ZZ(D)
        yield K


def generate_two_torsion_kernels(E):
    """
    Generates the kernels of order 2
    by generating a point of order
    dividing 2^k, and performing at most k
    doublings until we reach a point of order 2.

    Additionally, the only points of order
    2 are P, Q and P + Q

    So if we find P, and another point K
    of order 2 then K is either Q or P + Q
    such that K + P is either Q + P or Q

    This skips needing to compute pairings
    for linear independence
    """
    p = E.base().characteristic()

    # Compute odd cofactor
    n = p**2 - 1
    f = n.valuation(2)
    n = n // 2**f

    for _ in range(1000):
        P = n * E.random_point()
        if not P.is_zero():
            Pdbl = 2 * P
            while not Pdbl.is_zero():
                P = Pdbl
                Pdbl = 2 * P
            yield P
            break

    for _ in range(1000):
        Q = n * E.random_point()
        if not Q.is_zero():
            Qdbl = 2 * Q
            while not Qdbl.is_zero():
                Q = Qdbl
                Qdbl = 2 * Q
            if Q != P:
                yield Q
                break
    yield P + Q


def generate_kernels_division_polynomial(E, l):
    """
    Generate all kernels which generate cyclic isogenies
    of degree l from the curve E.

    Kernel generators are found by computing the roots
    of the l-th division polynomial and lifting these values
    to points on the elliptic curve.
    """
    f = E.division_polynomial(l)
    xs = [x for x, _ in f.roots()]

    for x in xs:
        K = E.lift_x(x)
        K._order = ZZ(l)
        yield K

# ===================================== #
#  Fast DLP solving using Weil pairing  #
# ===================================== #

def DLP(P, G, D):
    """
    Given two elliptic curve points G, P
    such that P = [x]G find x by using
    Weil pairings to allow the dlog to
    be performed over Fp4

    This can be between 30-100% faster than
    calling P.discrete_log(Q).

    TODO:
    Only supported for prime powers.
    """
    # Find a linearly independent point for
    # Weil pairing
    Q = generate_linearly_independent_point(G.curve(), G, D)

    # Map from elliptic curves to Fp4
    g = G.weil_pairing(Q, D, algorithm="pari")
    p = P.weil_pairing(Q, D, algorithm="pari")

    return discrete_log_pari(p, g, D)


def BiDLP(R, P, Q, D):
    """
    Given a basis P,Q of E[D] finds
    a,b such that R = [a]P + [b]Q.

    Uses the fact that 
        e([a]P + [b]Q, [c]P + [d]Q) = e(P,Q)^(ad-bc)
    """
    # e(P,Q)
    pair_PQ = P.weil_pairing(Q, D, algorithm="pari")

    # Write R = aP + bQ for unknown a,b
    # e(R, Q) = e(P, Q)^a
    pair_a = R.weil_pairing(Q, D, algorithm="pari")

    # e(R,-P) = e(P, Q)^b
    pair_b = R.weil_pairing(-P, D, algorithm="pari")

    # Now solve the dlog on Fq
    a = discrete_log_pari(pair_a, pair_PQ, D)
    b = discrete_log_pari(pair_b, pair_PQ, D)

    return a, b

# ===================================== #
#   Compute a kernel from an isogeny    #
# ===================================== #

def kernel_from_isogeny_prime_power(ϕ):
    """
    Given a prime-power degree isogeny ϕ
    computes a generator of its kernel
    """
    E = ϕ.domain()
    D = ϕ.degree()

    # Deal with isomorphisms
    if D == 1:
        return E(0)

    # Generate a torsion basis of E[D]
    P, Q = torsion_basis(E, D)

    # Compute the image of P,Q
    imP, imQ = ϕ(P), ϕ(Q)

    # Ensure we can use Q as a base
    # for the discrete log
    # TODO:
    # Here we assume D is a prime power,
    # this could fail for the general case
    if not has_order_D(imQ, D):
        P, Q = Q, P
        imP, imQ = imQ, imP

    # Solve the discrete log such that
    # imP = -[x]imQ. `DLP` uses Weil 
    # pairing to shift the dlog to Fp4.
    x = DLP(-imP, imQ, D)

    return P + x * Q
