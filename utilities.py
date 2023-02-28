"""
Helper functions for SQISign which seemed general enough to 
be utilities.

Includes:

- Cornacchia's algorithm for expressing integers as x^2 + d*y^2, 
- Algorithms for finding inert primes
- An efficient function for determining whether an element has order D, 
  suitable for additive and multiplicative groups.
"""

# Sage imports
from sage.all import (
    random_prime,
    factor,
    gp,
    ZZ,
    prod,
    is_pseudoprime,
    kronecker,
    cached_function,
    gcd,
    randint,
    choice,
    two_squares
)

# ======================================= #
#  Cornacchia's alg. and helper functions #
# ======================================= #

def sum_of_squares_friendly(n):
    """
    We can write any n = x^2 + y^2 providing that there
    are no prime power factors p^k | n such that 
    p ≡ 3 mod 4 and k odd.

    We can reject bad n by checking these conditions.
    """
    # We don't want any prime factors in n such that p ≡ 3 mod 4
    # Technically, we could allow these if they appear with even
    # exp. but this seems good enough.
    # We take the produce p ≡ 3 mod 4 for p < 500
    bad_prime_prod = 9758835872005600424541432018720696484164386911251445952809873289484557900064682994157437552186699082331048989
    if gcd(bad_prime_prod, n) != 1:
        return False

    # Now we consider the odd part of n and try and determine
    # if there are bad factors
    n_val = n.valuation(2)
    n_odd = n // (2**n_val)

    # First we remove all good prime factors p ≡ 1 mod 4
    # We take the produce p ≡ 1 mod 4 for p < 500
    good_prime_prod = 6396589037802337253083696670901044588941174775898621074430463772770231888063102454871836215129505
    g = gcd(good_prime_prod, n_odd)
    while g != 1:
        n_odd = n_odd // g
        g = gcd(n_odd , g)

    # Whatever is left is either composite, or a large ish-prime
    # If n_odd % 4 == 3, then there is certainly a prime factor 
    # p^k | n_odd such that p = 3 mod 4 and k is odd, so this is
    # no good
    if n_odd % 4 == 3:
        return False

    # We now have two cases, either n_odd is a prime p ≡ 1 mod 4
    # and we have a solution, or n_odd is a composite which is 
    # good, or has an even number of primes p ≡ 3 mod 4. The second
    # case is bad and we will only know by fully factoring it, so
    # to avoid the expensive factoring, we will reject these cases
    return is_pseudoprime(n_odd)

def cornacchia_friendly(n, B=2**20, d=None):
    """
    Checks whether the input for Cornacchia 
    is relatively easy to factor, to ensure it's
    fast to solve / reject.
    """
    if n < 0:
        return False

    # When we have to find m = x^2 + y^2, which
    # is the case for SQISign, we can perform a
    # few additional checks which seems to help
    # with performance
    if d == 1:
        return sum_of_squares_friendly(n)

    if n < 2**160:
        return True
    
    # Factor n by trial division up to B.
    # p is the part that did not get factored
    p = factor(n, limit=B)[-1][0]
    return p < 2**160 or is_pseudoprime(p)


def Cornacchia(M, d):
    """
    Given integers M,d computes
    x,y such that

    M = x^2 + dy^2
    """
    # First test if number is cornacchia friendly (cornacchia factors the number)
    if not cornacchia_friendly(M, d=d):
        # Empty list is what Q.qfbsolve
        # returns if it fails, so lets be
        # consistent I guess :)
        return []

    # When d=1 we can use SageMath's two_squares
    # which seems faster than calling Q.qfbsolve(ZZ(M))
    if d == 1:
        # two_squares raises value errors on bad
        # M, so we need to catch these
        try:
            sol = two_squares(M)
        except:
            sol = []
        return sol

    # General case:
    # Use Pari to define quadratic binary form
    # f(x,y) = x^2 + d*y^2
    Q = gp.Qfb(ZZ(1), ZZ(0), ZZ(d))
    sol = Q.qfbsolve(ZZ(M))
    # Convert from Pari interface to Sage Integer
    return list(map(ZZ, sol))

# ======================================= #
#   Compute random inert primes in R[ω]   #
# ======================================= #

def is_inert(ω, p):
    """
    Given an element ω ∈ B, and a prime p, returns if
    p is inert in ℤ[ω]

    For ℤ[ω], a prime is inert when the Kronecker symbol:
        (-d / p) = -1
    """

    def discriminant(ω):
        """
        Computes the discriminant of an element of a
        quaternion algebra, denoted Δ(ω)
        """
        Trd = ω.reduced_trace()
        Nrd = ω.reduced_norm()
        if Trd not in ZZ or Nrd not in ZZ:
            raise ValueError(f"The element ω must be integral")

        return ZZ(Trd**2 - 4 * Nrd)

    def recover_d(ω):
        """
        Given an integral element of a Quaternion Algebra,
        and its discriminant Δ(ω), computes the integer d
        such that:

        Δ(ω) = -d    if d ≡ 3 mod 4
             = -4d   otherwise
        """
        Δ = discriminant(ω)
        if Δ % 4 == 0:
            return -(Δ // 4)
        return -Δ

    d = recover_d(ω)
    return kronecker(-d, p) == -1


def inert_prime(bound, d):
    """
    Input: An upper bound for the output prime
           d, the integer such that ω^2 = -d for
           ℤ[ω]

    Output: A prime < bound which is inert in
            the ring ℤ[ω]
    """
    while True:
        p = random_prime(bound)
        if kronecker(-d, p) == -1:
            return p


# =============================================== #
#  Random selection, adapted from the Magma impl. #
#  https://github.com/SQISign/sqisign-magma       #
# =============================================== #

def generate_bounded_divisor(bound, T, facT):
    """
    Finds L2 > bound powersmooth
    We want the result to divide T
    """
    L2, facL2 = 1, dict()

    # Compute a random factor L2 | T
    # and create dict to keep track 
    # of exp.
    for pi, ei in facT:
        exp = ei - randint(0, ei)
        L2 *= pi**exp
        facL2[pi] = exp

    # Mul random factors from T until L2 
    # is large enough
    while L2 <= bound:
        pi, ei = choice(facT)
        if facL2[pi] < ei:
            L2 *= pi
            facL2[pi] += 1

    return ZZ(L2)

# ================================================== #
#  Code to check whether a group element has order D #
# ================================================== #

def batch_cofactor_mul_generic(G_list, pis, group_action, lower, upper):
    """
    Input:  A list of elements `G_list`, such that
                G is the first entry and the rest is empty
                in the sublist G_list[lower:upper]
            A list `pis` of primes p such that
                their product is D
            The `group_action` of the group
            Indices lower and upper
    Output: None
`
    NOTE: G_list is created in place
    """

    # check that indices are valid
    if lower > upper:
        raise ValueError(f"Wrong input to cofactor_multiples()")

    # last recursion step does not need any further splitting
    if upper - lower == 1:
        return

    # Split list in two parts,
    # multiply to get new start points for the two sublists,
    # and call the function recursively for both sublists.
    mid = lower + (upper - lower + 1) // 2
    cl, cu = 1, 1
    for i in range(lower, mid):
        cu = cu * pis[i]
    for i in range(mid, upper):
        cl = cl * pis[i]
    # cl = prod(pis[lower:mid])
    # cu = prod(pis[mid:upper])

    G_list[mid] = group_action(G_list[lower], cu)
    G_list[lower] = group_action(G_list[lower], cl)

    batch_cofactor_mul_generic(G_list, pis, group_action, lower, mid)
    batch_cofactor_mul_generic(G_list, pis, group_action, mid, upper)


@cached_function
def has_order_constants(D):
    """
    Helper function, finds constants to
    help with has_order_D
    """
    D = ZZ(D)
    pis = [p for p, _ in factor(D)]
    D_radical = prod(pis)
    Dtop = D // D_radical
    return Dtop, pis


def has_order_D(G, D, multiplicative=False):
    """
    Given an element G in a group, checks if the
    element has order exactly D. This is much faster
    than determining its order, and is enough for 
    many checks we need when computing the torsion
    bonus.

    We allow both additive and multiplicative groups
    """
    # For the case when we work with elements of Fp^k
    if multiplicative:
        group_action = lambda a, k: a**k
        is_identity = lambda a: a.is_one()
        identity = G.parent()(1)
    # For the case when we work with elements of E / Fp^k
    else:
        group_action = lambda a, k: k * a
        is_identity = lambda a: a.is_zero()
        identity = G.curve()(0)

    if is_identity(G):
        return False

    D_top, pis = has_order_constants(D)

    # If G is the identity after clearing the top
    # factors, we can abort early
    Gtop = group_action(G, D_top)
    if is_identity(Gtop):
        return False

    G_list = [identity for _ in range(len(pis))]
    G_list[0] = Gtop

    # Lastly we have to determine whether removing any prime 
    # factors of the order gives the identity of the group
    if len(pis) > 1:
        batch_cofactor_mul_generic(G_list, pis, group_action, 0, len(pis))
        if not all([not is_identity(G) for G in G_list]):
            return False

    return True


# ================ #
#  Print Debugging #
# ================ #

def print_info(str, banner="="):
    """
    Print information with a banner to help
    with visibility during debug printing
    """
    print(banner * 80)
    print(f"{str}".center(80))
    print(banner * 80)