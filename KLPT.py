"""
A KLPT implementation specialised for use for SQISign.

Many algorithms are taken from the from the SQISign paper:

    SQISign: compact post-quantum signatures from quaternions and isogenies
    Luca De Feo, David Kohel, Antonin Leroux, Christophe Petit, and Benjamin Wesolowski,
    https://ia.cr/2020/1240

With additional understanding and sub-algorithms from the original KLPT
paper:

    On the quaternion l-isogeny path problem
    David Kohel, Kristin Lauter, Christophe Petit, Jean-Pierre Tignol
    https://ia.cr/2014/505

Generally, we expect only to import EquivalentSmoothIdealHeuristic 
and SigningKLPT which are different forms of the KLPT algorithm used 
throughout SQISign.

As an implementation note, we try and append "Heuristic" to function 
names when we only know solutions can be derived heuristically. If no 
solution is found, `None` is returned.
"""

# Python imports
from random import choice

# Sage imports
from sage.all import (
    gcd,
    ZZ,
    Zmod,
    log,
    ceil,
    floor,
    sqrt,
    flatten,
    factor,
    Matrix,
    vector,
    randint,
    inverse_mod,
    prod,
    choice,
    CRT,
)

# Local imports
from ideals import (
    chi,
    is_integral,
    is_cyclic,
    reduced_basis,
    pullback_ideal,
    pushforward_ideal,
    equivalent_left_ideals,
    quaternion_basis_gcd,
    eichler_order_from_ideal,
    ideal_generator,
    make_cyclic,
    small_equivalent_ideal,
    quadratic_norm,
)
from utilities import Cornacchia, generate_bounded_divisor, is_inert
from lattices import generate_close_vectors, generate_small_norm_quat
from setup import (
    O0,
    p,
    q,
    j,
    ω,
    l,
    e,
    logp,
    loglogp,
    prime_norm_heuristic,
    represent_heuristic,
)

# ========================================== #
# Functions for solving finding equivalent   #
# prime norm ideals                          #
# ========================================== #


def generate_small_norm_quat_random(Ibasis, coeff_bound, search_bound):
    """
    Pick a random linear combination from Ibasis to compute an element
    α ∈ B_{0, ∞}
    """
    for _ in range(search_bound):
        xs = [randint(-coeff_bound, coeff_bound) for _ in range(len(Ibasis))]
        if gcd(xs) != 1:
            continue
        α = sum(αi * x for αi, x in zip(Ibasis, xs))

        yield α


def prime_norm_algebra_element(
    nI,
    Ibasis,
    coeff_bound,
    search_bound,
    previous=set(),
    allowed_factors=None,
    random_elements=False,
):
    """
    Find an element α ∈ B_{0, ∞} with small,
    prime scaled norm.

    Optional: `allowed_factors` allows the norm to
               be composite, where it is expected that
               the result is a large prime multiplied by
               small factors dividing allowed_factors.
    """

    if random_elements:
        small_elements = generate_small_norm_quat_random(
            Ibasis, coeff_bound, search_bound
        )
    else:
        max_norm_bound = prime_norm_heuristic
        small_elements = generate_small_norm_quat(
            Ibasis, max_norm_bound, count=search_bound
        )

    for α in small_elements:
        α_norm = ZZ(α.reduced_norm()) // nI

        # Even norms can be rejected early
        # as allowed_factors is either odd
        # or None.
        if α_norm % 2 == 0:
            continue

        # We can allow α to have composite norm
        # if small factors are within the KLPT
        # target norm T.
        α_norm_reduced = α_norm
        if allowed_factors:
            g = gcd(α_norm, allowed_factors)
            if g != 1:
                α_norm_reduced //= g

        # If we've failed with this norm before
        # continue
        if α_norm_reduced in previous:
            continue

        # Check if the element has prime norm
        if α_norm_reduced.is_pseudoprime():
            # Check if the prime small enough
            if α_norm_reduced < prime_norm_heuristic:
                return α, α_norm
    return None, None


def EquivalentPrimeIdealHeuristic(
    I, previous=set(), allowed_factors=None, random_elements=False
):
    """
    Given an ideal I with norm nI, attempts
    to find an equivalent ideal J with prime norm.

    If unsuccessful, returns None
    """
    # TODO: what's a good initial small search?
    coeff_bound = max((floor(logp / 10)), 7)

    # TODO: what's a good search bound?
    search_bound = max(coeff_bound**4, 4096)

    # Norm of Ideal
    nI = ZZ(I.norm())

    # Compute the Minkowski reduced basis
    Ibasis = reduced_basis(I)

    # Find an element with small prime norm
    α, N = prime_norm_algebra_element(
        nI,
        Ibasis,
        coeff_bound,
        search_bound,
        previous=previous,
        allowed_factors=allowed_factors,
        random_elements=random_elements,
    )

    if α is None:
        print(f"DEBUG [EquivalentPrimeIdealHeuristic] No equivalent prime found")
        return None, None, None

    assert ZZ(α.reduced_norm()) // nI == N
    assert α in I

    # Compute the ideal given α
    J = chi(α, I)

    return J, N, α


def RepresentIntegerHeuristic(M, parity=False):
    """
    Algorithm 1 (Page 8)

    Given an integer M, with M > p, attempts to
    find a random element γ with norm M.

    If no element is found after `bound` tries,
    returns none
    """

    def RepresentInteger(M, z, t, parity=False):
        M_prime = M - p * quadratic_norm(z, t)
        two_squares = Cornacchia(M_prime, -ZZ(ω**2))
        if two_squares:
            x, y = two_squares
            if parity and (x + t) % 2 == 0 and (y + z) % 2 == 0:
                return None
            return x + ω * y + j * (z + ω * t)
        # No solution for the given M
        return None

    if M <= p:
        raise ValueError(f"Can only represent integers M > p.")
    m = max(floor(sqrt(M / (p * (1 + q)))), 5)

    # TODO: how many times should we try?
    for _ in range(m**2):
        z = randint(-m, m)
        t = randint(-m, m)
        γ = RepresentInteger(M, z, t, parity=parity)

        if γ is not None:
            # Found a valid solution, return
            assert γ.reduced_norm() == M, "The norm is incorrect"
            assert γ in O0, "The element is not contained in O0"
            return γ

    # No solution found, return None
    print(f"DEBUG [RepresentIntegerHeuristic]: No solution found")
    return None


def EquivalentRandomEichlerIdeal(I, Nτ):
    """
    Algorithm 6 (SQISign paper)

    Input:  I a left O-ideal
    Output: K ∼ I of norm coprime with Nτ
    """
    nI = I.norm()
    # TODO: what should the size of `bound` be
    bound = 10

    # Step 1: find an element ωS such that Nτ is inert in ℤ[ωS]
    O = I.left_order()
    while True:
        ωS = sum([randint(-bound, bound) * b for b in O.basis()])
        if is_inert(ωS, Nτ):
            break

    # Step 2: find a random element γ in I such that n(γ)/n(I) is coprime with Nτ
    while True:
        γ = sum([randint(-bound, bound) * b for b in I.basis()])
        if gcd(γ.reduced_norm() // nI, Nτ) == 1:
            break

    # Step 3: select a random class (C : D) ∈ P1(Z/Nτ Z).
    x = randint(0, Nτ)
    if x == p:
        C, D = 1, 0
    else:
        C, D = x, 1

    # Step 4: set β = (C + ωSD)γ.
    β = (C + ωS * D) * γ

    # Step 5: return K = χI(β)
    return chi(β, I)


# ========================================== #
# Functions for solving the (Ideal/Eichler)  #
# Mod Constraint                             #
# ========================================== #


def solve_mod_constraint_kernel(quat_list):
    """
    Helper function which given a list of 8 quaternions
    constructs a matrix and computes its kernel.

    Used in `IdealModConstraint` and `EichlerModConstraint`
    to compute C, D.
    """
    # Split each element of B into its coefficients
    # a ∈ B = t + x i + y j + z k: t,x,y,x ∈ Q
    matrix_coefficients = flatten([x.coefficient_tuple() for x in quat_list])

    # We now have 32 elements, four coefficients from each
    # element
    assert len(matrix_coefficients) == 8 * 4

    # Create an 8x4 matrix over QQ from the coefficients
    mat = Matrix(8, 4, matrix_coefficients)

    # Clear the denominator and work over ZZ
    mat, _ = mat._clear_denom()

    # Compute its kernel
    kernel = mat.kernel()

    # Following the SQISign MAGMA code, we pick
    # the third row.
    row = kernel.basis()[2]
    C, D = row[2:4]  # row = A,B,C,D,_,_,_,_

    return C, D


def IdealModConstraint(I, γ):
    """
    Section 4.3

    Given an ideal I of norm N and an element
    γ of a maximal order O0 of norm NM finds
    (C0 : D0) in P^1(Z/NZ) such that

    μ0 = j(C0 + ωD0)

    Then verifies that γ*μ0 in I
    """

    N = I.norm()
    # First four elements constructed from γ
    matrix_quaternions = [γ * N, γ * N * ω, γ * j, γ * j * ω]

    # Next four elements are the basis of the ideal I
    matrix_quaternions += list(I.basis())

    # We should now have a list of 8 elements of B
    assert len(matrix_quaternions) == 8

    # Generate a 8x4 matrix matrix from these quaternions
    C0, D0 = solve_mod_constraint_kernel(matrix_quaternions)
    μ0 = j * (C0 + ω * D0)

    # Check we made it correctly
    assert γ * μ0 in I

    if C0 == 0 or D0 == 0:
        print(f"DEBUG [IdealModConstraint]: Something is maybe going to break...")
        print(f"{C0 = }")
        print(f"{D0 = }")
    return C0, D0


def EichlerModConstraint(L, EichlerL, γ1, γ2):
    """
    Input:  An ideal L, its corresponding Eichler order
            ℤ + L and two elements of a quaternion algebra
            γ1, γ2
    Output: (C : D) ∈ P1(Z/n(L)Z) such that the element
            μ1 = j(C + ωD)
            γ1*μ1*γ2 ∈ Z + L

    Taken from Section 6.2

    Given an ideal L of norm N and two elements
    γ1, γ2 we find μ1 such that γ1*μ1*γ2 is in
    the Eichler Order.
    """
    nL = ZZ(L.norm())
    assert gcd(γ1.reduced_norm(), nL) == 1
    assert gcd(γ2.reduced_norm(), nL) == 1

    # Construct the matrix, our C,D are extracted from the kernel.
    # First four elements constructed from γ
    matrix_quaternions = [γ1 * nL, γ1 * nL * ω, γ1 * j, γ1 * j * ω]

    # Next four elements are the basis of the ideal I
    matrix_quaternions += [b * γ2.conjugate() for b in EichlerL.basis()]

    # We should now have a list of 8 elements of B
    assert len(matrix_quaternions) == 8

    C1, D1 = solve_mod_constraint_kernel(matrix_quaternions)
    μ1 = j * (C1 + ω * D1)

    assert γ1 * μ1 * γ2 in EichlerL
    return C1, D1


def check_ideal_mod_constraint(L, γ, C, D):
    """
    For the case when we work with prime
    power norm (we currently only support 2^k)
    we need to make sure either C or D can be
    invertible.

    This function tries to correct for even factors
    while maintaining γ*μ0 is in L.
    """
    # When both C and D are even we can
    # try and scale to make at most one
    # odd.
    while C % 2 == 0 and D % 2 == 0:
        C = C // 2
        D = D // 2
        μ0 = j * (C + ω * D)
        if γ * μ0 not in L:
            # Could not halve C, D and
            # maintain γ*μ in L
            return False, None, None

    # Nrd(μ0) must be coprime with N for
    # inversion
    μ0 = j * (C + ω * D)
    if gcd(μ0.reduced_norm(), L.norm()) != 1:
        return False, None, None

    return True, C, D


# ============================================== #
# Functions for solving the Strong Approximation #
# ============================================== #


def strong_approximation_construct_lattice(N, C, D, λ, L, small_power_of_two=False):
    """
    Constructing the lattice basis and target vector following PS18

    When N is a power of two, we have `small_power_of_two`
    set to True. When this is the case we double the
    modulus and halve coeff_z, coeff_t
    """

    if small_power_of_two:
        coeff_z = p * λ * C
        coeff_t = p * λ * D
    else:
        coeff_z = 2 * p * λ * C
        coeff_t = 2 * p * λ * D

    cst_term = L - p * λ**2 * quadratic_norm(C, D)

    if small_power_of_two:
        cst_term, check = cst_term.quo_rem(2 * N)
    else:
        cst_term, check = cst_term.quo_rem(N)
    assert check == 0

    coeff_t_inv = inverse_mod(coeff_t, N)

    zp0 = 0
    tp0 = ZZ(cst_term * coeff_t_inv % N)

    # Construct the lattice basis
    lattice_basis = N * Matrix(ZZ, [[1, ZZ((-coeff_z * coeff_t_inv) % N)], [0, N]])

    # Construct the target vector
    target = λ * vector([ZZ(C), ZZ(D)]) + N * vector([zp0, tp0])
    return lattice_basis, target, zp0, tp0


def strong_approximation_lattice_heuristic(N, C, D, λ, L, small_power_of_two=False):
    """
    Constructs a lattice basis and then looks for
    close vectors to the target.

    Allows for optimising output from pN^4 to pN^3,
    which helps keep the norm small and hence the
    degree of the isogenies small
    """

    # We really only expect this for the case when N = 2^k
    swap = False
    if D == 0 or gcd(D, N) != 1:
        C, D = D, C
        swap = True

    # Construct the lattice
    lattice_basis, target, zp0, tp0 = strong_approximation_construct_lattice(
        N, C, D, λ, L, small_power_of_two=small_power_of_two
    )

    # Generate vectors close to target
    close_vectors = generate_close_vectors(lattice_basis, -target, p, L)

    xp, yp = None, None
    for close_v in close_vectors:
        zp, tp = close_v
        assert zp % N == 0, "Can't divide zp by N"
        assert tp % N == 0, "Can't divide tp by N"

        zp = ZZ(zp / N) + zp0
        tp = ZZ(tp / N) + tp0
        M = L - p * quadratic_norm(λ * C + zp * N, λ * D + tp * N)
        M, check = M.quo_rem(N**2)
        assert check == 0, "Cant divide by N^2"

        if M < 0:
            continue

        # Try and find a solution to
        # M = x^2 + y^2
        two_squares = Cornacchia(ZZ(M), -ZZ(ω**2))
        if two_squares:
            xp, yp = two_squares
            break

    if xp is None:
        # Never found vector which had a valid solution
        return None

    # Use solution to construct element μ
    # μ = λ*j*(C + D*ω) + N*(xp + ω*yp + j*(zp + ω*tp))

    # If we swapped earlier, swap again!
    if swap:
        C, D = D, C
        tp, zp = zp, tp

    μ = N * xp + N * yp * ω + (λ * C + N * zp) * j + (λ * D + N * tp) * j * ω

    # Check that Nrd(μ) == L
    # and that μ is in O0
    assert μ.reduced_norm() == L
    assert μ in O0
    return μ


def StrongApproximationHeuristic(
    N, C, D, facT, composite_factors=None, small_power_of_two=False
):
    """
    Algorithm 2 (Page 9)

    Given an N such and two integers C, D such that
    μ0 = j(C + ωD), we find a
      μ = λμ0 + N μ1
    such that the norm of μ is a divisor of facT

    This function works when:
    - Case 1: When N is prime
    - Case 2: When N is a power of two
    - Case 3: When N is composite.

    For Case 2, the optional bool small_power_of_two must be True
    For Case 3, the factors of N must be included in the
    optional argument: composite_factors
    """

    # For Case 2
    if small_power_of_two:
        # When we work with N = 2^k, we
        # work mod 2N and halve the lattice
        K = Zmod(2 * N)

    # For Case 1 and Case 3
    else:
        K = Zmod(N)

    # Check whether T is large enough to find
    # solutions
    T = prod([l**e for l, e in facT])

    # Case 3: N = N1*N2
    # For the case in SigningKLPT when N is the product of two
    # primes, we have a fixed L2. So, we try this and return
    # None if it fails
    if composite_factors:
        NL, Nτ = composite_factors
        assert NL * Nτ == N, "Supplied factors are incorrect"

        # Bound is fixed
        L2 = T

        # Make sure we can compute λ^2
        try:
            λλ = K(L2) / K(p * quadratic_norm(C, D))
        except Exception as e:
            # p*(C^2 + D^2) is not invertible mod N
            print(f"ERROR [StrongApproximationHeuristic]: {e}")
            return None

        # Recover the square root with CRT
        # We supply C, D such that we can take the sqrt
        # no need to check this is possible
        λL = ZZ(Zmod(NL)(λλ).sqrt())
        λτ = ZZ(Zmod(Nτ)(λλ).sqrt())
        λ = CRT([λL, λτ], [NL, Nτ])

        # Not sure why, but this sometimes happens...
        if λ == 0:
            return None

        # If the lattice computation fails this returns None
        return strong_approximation_lattice_heuristic(
            N, C, D, λ, L2, small_power_of_two=small_power_of_two
        )

    # Case 1: N is prime
    # Case 2: N = l^2
    # Otherwise, we can pick a bunch of L2 and check if
    # any of them work.
    tested = set()

    # Debugging make sure the bound is ok for given T
    bound = ceil((p * N**3))
    if bound > T:
        print(
            f"DEBUG [StrongApproximationHeuristic] The prime norm is too large, "
            f"no valid divisors. N ~ 2^{N.nbits()}, T ~ 2^{T.nbits()}"
        )
        return None

    # Try and find μ by testing a bunch of different L2 | T
    for _ in range(100):
        L2 = generate_bounded_divisor(bound, T, facT)
        if L2 in tested:
            continue
        tested.add(L2)

        # Check given L2 whether we can compute λ
        try:
            λλ = K(L2) / K(p * quadratic_norm(C, D))
        except:
            # p*(C^2 + D^2) is not invertible mod N
            # Not sure why this is happening...
            return None

        # Ensure we can compute λ from λ^2
        if not λλ.is_square():
            continue

        λ = ZZ(λλ.sqrt())

        # Not sure why, but this sometimes happens...
        if λ == 0:
            continue

        μ = strong_approximation_lattice_heuristic(
            N, C, D, λ, L2, small_power_of_two=small_power_of_two
        )
        # Sometimes we have no valid solutions from the lattice,
        # so we try again with a new L2.
        if μ:
            return μ
    # All attempts for L2 failed. Pick a new N by regenerating γ
    return None


# ====================================== #
# KLPT algorithm, specialised for SQISign #
# ====================================== #

# ================ #
# Helper functions #
# ================ #


def compute_L1_KLPT(N, facT):
    """
    Helper function for KLPT which is
    used to find necessary output size of
    `RepresentIntegerHeuristic`.
    """
    L1 = 1
    while N * L1 < represent_heuristic * p:
        l, e = choice(facT)
        L1 *= l
        facT.remove((l, e))
        if e > 1:
            facT.append((l, e - 1))
    return L1, facT


def equivalent_prime_ideal_wrapper_KLPT(
    I,
    T,
    previously_seen,
    small_power_of_two=False,
    equivalent_prime_ideal=None,
    allowed_factors=None,
):
    """
    Handles various cases for KLPT algorithm.

    Case 1: when the input has a small power of two norm,
    we can skip computing a prime ideal altogether

    Case 2: we already know a prime norm ideal, so we
    just return this

    Case 3: we need to compute a prime norm ideal,
    run EquivalentPrimeIdealHeuristic()

    Optional: when we supply allowed_factors as input,
    we allow the output to be a composite. We expect
    this to be some large prime with some other smaller
    factors which must divide allowed_factors.
    This makes it easier to find L, and we can correct
    for these small factors later.
    """
    L0 = 1
    # When we are working with nI = 2^k, with 2^k << p
    # we need to proceed differently.
    # Instead of computing an equivalent prime norm ideal L
    # we instead work with `I` directly.
    if small_power_of_two:
        # We make sure the ideal is a fixed point for the
        # action of (R/2R)^*
        i = I.quaternion_algebra().gens()[0]
        if I + O0 * 2 != O0 * (i + 1) + O0 * 2:
            I = O0.left_ideal([b * (i + 1) for b in I.basis()])

        # Use the input I as L
        L = I
        N = ZZ(L.norm())
        return L, N, L0, previously_seen

    # During keygen, we already know a good small norm ideal
    # which we can pass in to skip the below block.
    if equivalent_prime_ideal:
        L = equivalent_prime_ideal
        N = ZZ(L.norm())
        return L, N, L0, previously_seen

    # TODO: how many times should we try?
    # Could we just run this once by setting good inner bounds?
    for _ in range(10):
        L, N, _ = EquivalentPrimeIdealHeuristic(
            I, previous=previously_seen, allowed_factors=allowed_factors
        )

        # Found a suitable equivalent ideal
        if L is not None:
            break

    # If we never find a prime norm ideal, throw an error
    if L is None:
        raise ValueError(
            "Could not find a prime norm ideal, need bigger T or smaller norm output..."
        )

    previously_seen.add(N)

    if allowed_factors:
        # Absorb the small primes from D(T^2) into L0
        L0 = gcd(N, allowed_factors)
        # Restore N to be prime
        N //= L0
        assert N.is_pseudoprime(), "N is not prime"

    return L, N, L0, previously_seen


def strong_approximation_wrapper_KLPT(L, N, L0, facT, small_power_of_two=False):
    """
    Helper function:

    Given a prime norm ideal L, tries to solve
    the StrongApproximationPrime while obeying
    various bounds and edge cases for SQISign
    """

    # Check L, N, L0 all align
    assert L.norm() == N * L0

    # Find the prime ideal
    if L0 != 1:
        O = L.left_order()
        α = ideal_generator(L)
        L_prime = O * α + O * N
        assert L_prime.norm() == N, "Made the ideal in the wrong way..."
    else:
        L_prime = L

    # For the case of N = l^e, we often need
    # to pick many γ. Remember what we've seen
    # to avoid repeating bad choices.
    seen_γ = set()

    # TODO: how many times to try?
    for _ in range(10):
        # Find a factor L1 such that N*L1 > p to ensure
        # we can compute γ
        L1, facTupdated = compute_L1_KLPT(N, facT.copy())

        # Make sure that L1 isn't too big
        L2_max = prod([p**e for p, e in facTupdated])
        if L2_max < p * N**3:
            print(
                f"DEBUG [strong_approximation_wrapper_KLPT]:"
                "L1 is too big, likely no lattice solutions"
            )

        γ = RepresentIntegerHeuristic(N * L1)
        if γ is None:
            continue

        if γ in seen_γ:
            continue
        seen_γ.add(γ)

        C0, D0 = IdealModConstraint(L_prime, γ)

        if N % 2 == 0:
            # We are going to have to invert the element
            # p*(C^2 + D^2) mod N
            # As a result, we need to remove factors of 2
            # and ultimately reject some valid C0,D0 solutions
            check, C0, D0 = check_ideal_mod_constraint(L_prime, γ, C0, D0)
            # C0, D0 are bad, pick new γ
            if check == False:
                continue

        μ0 = j * (C0 + ω * D0)
        assert μ0 in O0, "μ0 is not contained within O0"
        assert γ * μ0 in L_prime, "The product γ * μ0 is not contained within L"

        ν = StrongApproximationHeuristic(
            N, C0, D0, facTupdated, small_power_of_two=small_power_of_two
        )
        if ν is not None:

            β = γ * ν
            if L0 == 1:
                # Dumb checks...
                assert L == L_prime
                assert β in L, "β is not in prime norm ideal: L"
                assert β.reduced_norm() % N == 0

                # Compute the equivalent ideal
                return chi(β, L)

            # For the near prime case we need to
            # adjust everything so that we have an element
            # contained in L, not L_prime

            # Li_product = L0 * L1 * L2
            Li_product = L0 * (β.reduced_norm() // N)
            δ = ideal_generator(L, coprime_factor=Li_product)
            β = δ * β.conjugate()
            O = L.left_order()
            J = O * β + O * Li_product

            return J
    # Never found a solution, give up so we can pick a new N
    return None


# ==================== #
# End Helper functions #
# ==================== #


def EquivalentSmoothIdealHeuristic(I, T, equivalent_prime_ideal=None, near_prime=False):
    """
    Algorithm 3 (KLPT Algorithm)

    Given an ideal I with left order O0 and a smooth integer T,
    computes an equivalent ideal J with norm dividing T together
    with the quaternion algebra element β.
    """
    # TODO: we could pass in the factors, rather than factoring?
    facT = list(factor(T))

    # Remember N we have found and skip them
    previously_seen = set()

    # Make I as small as possible:
    I = small_equivalent_ideal(I)

    # For case distinction for Inorm = 2^* and otherwise
    nI = ZZ(I.norm())

    small_power_of_two = False
    if nI % 2 == 0 and len(factor(nI)) == 1 and nI < prime_norm_heuristic:
        small_power_of_two = True

    allowed_factors = None
    if near_prime:
        allowed_factors = T

    # TODO: how many times should we try?
    for _ in range(40):
        # Find a prime norm, or prime power norm ideal L with norm N
        L, N, L0, previously_seen = equivalent_prime_ideal_wrapper_KLPT(
            I,
            T,
            previously_seen,
            small_power_of_two=small_power_of_two,
            equivalent_prime_ideal=equivalent_prime_ideal,
            allowed_factors=allowed_factors,
        )

        if L0 != 1:
            print(
                f"DEBUG [EquivalentSmoothIdealHeuristic]: "
                f"Working with a nearprime with L0 = {factor(L0)}"
            )
            # We've allowed non-prime ideals, make sure we made
            # good choices!
            assert N * L0 == L.norm(), "Something went wrong computing L0"
            assert T % L0 == 0, "Allowed factors do not divide the KLPT target norm"

            # Now we remove the factors from the available ones
            T = prod([p**e for p, e in facT])
            facT_updated = list(factor(T // L0))

        else:
            facT_updated = facT.copy()

        # We pass the strong approximation for N,
        # regardless of whether N is prime or 2^k
        # Logic for differences is inside the wrapper.
        J = strong_approximation_wrapper_KLPT(
            L, N, L0, facT_updated, small_power_of_two=small_power_of_two
        )

        if J is None:
            # Could not find a ν, pick a new N or just try again?
            print(
                f"DEBUG [EquivalentSmoothIdeal]:"
                "No solution found, trying again with new starting ideal"
            )
            continue

        # Some last minute checks...
        assert is_integral(J), "Output ideal is not integral"
        assert T % ZZ(J.norm()) == 0, "Ideal does not have target norm"
        assert equivalent_left_ideals(I, J)

        # J is an ideal equivalent to I with norm dividing the target T

        # Ensure that J is cyclic
        J, _ = make_cyclic(J)
        return J

    print(
        f"DEBUG [EquivalentSmoothIdeal]: No Equivalent Smooth Ideal could be computed..."
    )
    return None


# =================================================== #
# Signing KLPT, fixed norm output for SQISign Signing #
# =================================================== #

# ====================== #
# Begin Helper functions #
# ====================== #


def derive_L(I, Iτ, Nτ, O0, O1):
    """
    Given an ideal I with left order O1 and an ideal
    Iτ with left order O0 of norm Nτ computes an ideal
    equivalent to the pullback of I under Iτ with prime
    norm.

    Input: I with left order O1
           Iτ with left order O0 and norm Nτ

    Output L ~ [Iτ]^* I with prime norm
           N = n(L)
           δ such that L = χ(K', δ)
    """
    for _ in range(20):
        # The PoC implementation skips this, but it's
        # not computationally expensive, so we include
        # it anyway
        K = EquivalentRandomEichlerIdeal(I, Nτ)

        # Make K as small as possible
        K = small_equivalent_ideal(K)

        # K' = [Iτ]^* K
        K_prime = pullback_ideal(O0, O1, K, Iτ)

        L, N, δ = EquivalentPrimeIdealHeuristic(K_prime)

        # Bad delta, this will cause EichlerModConstraint to break
        if gcd(δ.reduced_norm(), Nτ) != 1:
            print(f"DEBUG [SigningKLPT]: Not sure why this is happening...")
            print(f"{factor(δ.reduced_norm()) = }")
            print(f"{factor(Nτ.reduced_norm()) = }")
            print(f"{gcd(δ.reduced_norm(), Nτ) = }")
            continue

        if L is not None:
            return L, N, δ

    # If we get here, something is likely broken
    raise ValueError(f"Never found an equivalent prime norm ideal")


def derive_L2_SigningKLPT(γ, L1, e1):
    """
    Given L1 = l^e1 and γ try and compute L2
    so that the output of SigningKLPT has norm
    exactly 2^e
    """
    g = quaternion_basis_gcd(γ, O0)
    extra = 2 * (floor(loglogp / 4) + ZZ(gcd(g, L1).valuation(l)))
    e2 = e - e1 + extra
    return l**e2


def derive_C_and_D_SigningKLPT(L, N, Iτ, Nτ, EichlerIτ, γ, δ, L2):
    """
    Solves IdealModConstraint and EichlerModConstraint
    for a given γ and returns when we find an element
    μ such that L2 / Nrd(μ) is a square mod N*Nτ

    Input: Ideals L, Iτ of prime norm N, Nτ
           EichlerIτ the Eichler order of Iτ
           γ, δ, elements of B
           L2, a divisor of T

    Output C, D such that L2 / p*(C^2 + D^2) is a square mod N*Nτ
    """

    C0, D0 = IdealModConstraint(L, γ)
    C1, D1 = EichlerModConstraint(Iτ, EichlerIτ, γ, δ)

    # Compute CRT
    C = CRT([C0, C1], [N, Nτ])
    D = CRT([D0, D1], [N, Nτ])

    # We need to take a sqrt of this
    μ_norm = p * quadratic_norm(C, D)

    KL = Zmod(N)
    Kτ = Zmod(Nτ)

    # TODO: Sometimes μ_norm is 0 mod N or Kτ
    #       Catch this earlier so this never happens...
    try:
        square_mod_N = (KL(L2) / KL(μ_norm)).is_square()
        square_mod_Nτ = (Kτ(L2) / Kτ(μ_norm)).is_square()
    except:
        return None, None
    # To compute the square root mod N*Nτ
    # Both of these must be true. Will happen
    # about 1/4 of the time.
    if square_mod_N and square_mod_Nτ:
        return C, D

    return None, None


# ====================== #
#  End Helper functions  #
# ====================== #


def SigningKLPT(I, Iτ):
    """
    Algorithm 5 (SQISign paper)

    Input: Iτ a left O0-ideal and right O-ideal of norm Nτ,
           I, a left O-ideal

    Output: J ∼ I of norm l^e, where e is fixed (global param)
    """
    assert is_cyclic(I), "I is not cyclic"
    assert is_cyclic(Iτ), "Iτ is not cyclic"

    # Prime norm ideal
    Nτ = ZZ(Iτ.norm())

    # Make I as small as possible
    I = small_equivalent_ideal(I)

    # Orders needed for pushback and pullforward
    O1 = I.left_order()
    assert Iτ.left_order() == O0
    assert Iτ.right_order() == O1

    # Compute the pullback K of I with left order O0, and
    # find an equivalent prime norm ideal L ~ K.
    L, N, δ = derive_L(I, Iτ, Nτ, O0, O1)

    # We want L1 to be big enough that we sensibly find solutions
    # for RepresentIntegerHeuristic(N * L1) but not so big that we
    # find no solutions to the lattice problem in the strong approx.
    e1 = floor(logp - log(N, l) + 1.74 * loglogp)
    L1 = l**e1

    # EichlerIτ = ℤ + Iτ = OL(I) ∩ OR(I)
    EichlerIτ = eichler_order_from_ideal(Iτ)

    # Store γ which appear to stop checking the element twice
    seen_γ = set()

    # TODO how many times should we try to solve the strong approx?
    for _ in range(2000):
        γ = RepresentIntegerHeuristic(N * L1)
        if γ is None:
            print(f"DEBUG [SigningKLPT]: Unable to compute a γ, trying again.")
            continue

        # No point trying the same element twice
        if γ in seen_γ:
            print(f"DEBUG [SigningKLPT]: Already tried γ, trying again.")
            continue
        seen_γ.add(γ)

        # If this GCD is non-trivial, EichlerModConstraint will break
        if gcd(γ.reduced_norm(), Nτ) != 1:
            continue

        # Given L1 and γ derive the bound L2. We can estimate how non-cyclic
        # the end result will be from the gcd of the elements of γ in the basis
        # of O0, but it's not perfect, so sometimes we need to run SigningKLPT
        # many times to ensure that n(J) = 2^1000
        L2 = derive_L2_SigningKLPT(γ, L1, e1)

        # Look for  Given L1 = l^e1 and γ try and compute L2
        # so that the output of SigningKLPT has norm
        # exactly 2^eμ = j(C + ωD) such that L2 / Nrd(μ) is a square mod N*Nτ
        C, D = derive_C_and_D_SigningKLPT(L, N, Iτ, Nτ, EichlerIτ, γ, δ, L2)
        if C is None:
            print(f"DEBUG [SigningKLPT]: Unable to compute a C,D, given γ.")
            continue

        # Search for a solution to the strong approximation. As the L2 is
        # fixed, we only try once.
        μ = StrongApproximationHeuristic(
            N * Nτ, C, D, factor(L2), composite_factors=(N, Nτ)
        )

        # No solution found, try another γ
        if not μ:
            continue

        # Strong approximation norm check
        assert μ.reduced_norm() == L2
        print(f"INFO [SigningKLPT]: Found a solution to the StrongApproximation!")

        # Now construct the equivalent ideal J
        β = γ * μ
        J_prime = chi(β, L)

        # J = [Iτ]_* J'
        J = pushforward_ideal(O0, O1, J_prime, Iτ)

        # Check things are actually equivalent
        assert equivalent_left_ideals(I, J)

        # Make sure the output is cyclic
        J, _ = make_cyclic(J)

        # Rarely, we will have estimated L2 wrong and
        # in the process of making J cyclic, computed
        # an ideal J with n(J) != l^e
        if J.norm() != l**e:
            print(f"DEBUG [SigningKLPT]: J has the wrong norm, trying again")
            print(f"DEBUG [SigningKLPT]: {factor(J.norm()) = }")
            continue

        return J
