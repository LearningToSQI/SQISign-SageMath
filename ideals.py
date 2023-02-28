"""
Helper functions for various computations associated to
the quaternion algebra, and ideals and orders of the 
quaternion algebra.

Some of these functions could be ported to SageMath. It's a TODO
for when SQISign is being less actively worked on.
"""

# Sage Imports
from sage.all import (
    factor,
    randint,
    ZZ,
    round,
    ceil,
    log,
    gcd,
    Matrix,
    vector,
    prod,
)

# Local imports
from setup import p, O0, ω
from utilities import cornacchia_friendly

# ================================================ #
#  Helpers for elements of the quaternion algebra  #
# ================================================ #

def quadratic_norm(x, y):
    """
    Given two integers x,y, which correspond 
    to the element x + ωy ∈ R[ω], returns
    Nrd(x + ωy)
    

    Note: This function implements the norm
    function f(x,y) in the SQISign papers.

    For SQISign, we have ω = i and i^2 = -1
    so f(x,y) = x**2 + y**2 which is more
    efficient
    """
    if ω**2 != -1:
        raise ValueError(f"quadratic_norm() requires Z[ω] with ω^2 = -1. {ω = }")
    return ZZ(x) ** 2 + ZZ(y) ** 2


def quaternion_change_basis(γ, O):
    """
    Computes the coefficients of a quaternion γ
    in the basis of a given order O
    """
    O_matrix = Matrix([b.coefficient_tuple() for b in O.basis()])
    γ_vector = vector(γ.coefficient_tuple())
    γO_coeffs = γ_vector * O_matrix.inverse()

    assert γ == sum([a * b for a, b in zip(γO_coeffs, O.basis())])
    return γO_coeffs


def quaternion_basis_gcd(γ, O):
    """
    Computes the gcd of the coefficients of a
    quaternion γ in the basis of a given order O
    """
    γO_coeffs = quaternion_change_basis(γ, O)
    return gcd(γO_coeffs)

# ============================================== #
#  Helpers for ideals of the quaternion algebra  #
# ============================================== #

def multiply_ideals(I, J, beta=None):
    """
    Computes I*J when O_R(I) ≃ O_L(J)

    If these orders do not match, we must provide
    an isomorphism which takes O_L(J) to O_R(I)
    """
    if I.right_order() != J.left_order():
        if beta is None:
            raise ValueError(
                "Right and left orders, do not match. Must supply an automorphism, beta"
            )

        J = beta ** (-1) * J * beta
        assert I.right_order() == J.left_order(), "Orders still do not match after applying automorphism"
    return I * J


def is_integral(I):
    """
    Checks whether the input ideal is integral.
    """
    return all([b in I.left_order() for b in I.basis()])


def ideal_basis_gcd(I):
    """
    Computes the gcd of the coefficients of
    the ideal written as a linear combination
    of the basis of its left order.
    """
    I_basis = I.basis_matrix()
    O_basis = I.left_order().unit_ideal().basis_matrix()

    # Write I in the basis of its left order
    M = I_basis * O_basis.inverse()
    return gcd((gcd(M_row) for M_row in M))


def is_cyclic(I):
    """
    Computes whether the input ideal is cyclic,
    all the work is done by the helper function
    `ideal_basis_gcd()`.
    """
    return ideal_basis_gcd(I) == 1


def remove_2_endo(J):
    """
    Helper function for `make_cyclic`. Not currently
    useful, but will be when we need to do SQISign2.
    """
    i, _, _ = J.quaternion_algebra().gens()
    two_endo = O0.left_ideal([i + 1, 2])
    while all(b in two_endo for b in J.basis()):
        J = J * J.right_order([(i - 1) / 2])
    return J


def make_cyclic(I, full=False):
    """
    Given an ideal I, returns a cyclic ideal by dividing
    out the scalar factor g = ideal_basis_gcd(I)
    """
    g = ideal_basis_gcd(I)
    # Ideal was already cyclic
    if g == 1:
        return I, g

    print(f"DEBUG [make_cyclic]: Ideal is not cyclic, removing scalar factor: {g = }")
    J = I.scale(1/g)

    if full:
        # TODO: will remove_2_endo change g?
        # not an issue currently, as we don't
        # use this.
        return remove_2_endo(J), g
    return J, g


def reduced_basis(I, check=False):
    """
    Computes the Minkowski reduced basis of the
    input ideal. Note: this produces the same 
    result for all ideals in the equivalence class
    so corresponds to the reduced basis of the
    smallest equivalent ideal to I
    
    Input: an ideal
    Output: A Minkowski-reduced basis

    Optional: when check is True, the basis is
              checked against the Minkowski bounds
    """

    def _matrix_to_gens(M, B):
        """
        Converts from a matrix to generators in the quat.
        algebra
        """
        return [sum(c * g for c, g in zip(row, B)) for row in M]

    B = I.basis()
    G = I.gram_matrix()
    U = G.LLL_gram().transpose()

    reduced_basis_elements = _matrix_to_gens(U, B)

    if check:
        norm_product = 16 * prod([x.reduced_norm() for x in reduced_basis_elements])
        tmp = p**2 * I.norm() ** 4
        assert norm_product <= 4 * tmp, "Minkowski reduced basis is too large"
        assert norm_product >= tmp, "Minkowski reduced basis is too small"

    return reduced_basis_elements


def small_equivalent_ideal(I, reduced_basis_elements=None):
    """
    Computes the Minkowski reduced basis of the
    ideal and returns the smallest equivalent
    ideal J = Iα ~ I.
    """
    nI = I.norm()

    if not reduced_basis_elements:
        reduced_basis_elements = reduced_basis(I)

    b0 = reduced_basis_elements[0]
    if b0.reduced_norm() == nI**2:
        return I
    return I * I.right_order().left_ideal([b0.conjugate() / nI])


def equivalent_left_ideals(I, J):
    """
    SageMath has this impl. for right ideals
    only. To work around this, we can first
    take the conjugate of I and J.

    TODO: write a function which does this without
          conjugates?
    """
    return I.conjugate().is_equivalent(J.conjugate())


def equivalent_right_ideals(I, J):
    """
    Sage can do this for us already only making a
    wrapped so it matches the above.
    """
    return I.is_equivalent(J)


def invert_ideal(I):
    """
    Computes the inverse of the ideal which is
    the conjugate of I divided by its norm
    """
    return I.conjugate().scale(1 / I.norm())


def left_isomorphism(I, J):
    """
    Given two isomorphic left ideals I, J computes
    α such that J = I*α
    """
    B = I.quaternion_algebra()

    if B != J.quaternion_algebra():
        raise ValueError("Arguments must be ideals in the same algebra.")

    if I.left_order() != J.left_order():
        raise ValueError("Arguments must have the same left order.")

    IJ = I.conjugate() * J
    L = reduced_basis(IJ)
    for t in L:
        α = t / I.norm()
        if J == I * α:
            return α

    raise ValueError("Could not find a left isomorphism...")


def chi(a, I):
    r"""
    From section 3.2 in Antonin's thesis.
    Calculates the equivalent ideal of I, of norm(a)
    Based on the surjection from I \ {0} to the set of equivalent ideals of I
    Obtained by a → I * (a_bar / n(I))
    """
    return I * (a.conjugate() / I.norm())


def chi_inverse(I, J):
    """
    Computes the element α such that

    J = Chi(α, I)
    """
    # Compute the element α
    a = left_isomorphism(I, J)
    assert J == I * a, "Left isomorphism to find the element 'a' failed... why?"
    α = (a * I.norm()).conjugate()
    assert J == chi(α, I), "Something about α is wrong!"
    return α


def scaled_norm(a, I):
    """
    Returns Nrd(a) / n(I), the norm of chi(a, I)
    """
    N = I.norm()
    return a.reduced_norm() / N


def ideal_generator(I, coprime_factor=1):
    """
    Given an ideal I of norm D, finds a generator
    α such that I = O(α,D) = Oα + OD

    Optional: Enure the norm of the generator is coprime 
    to the integer coprime_factor
    """
    OI = I.left_order()
    D = ZZ(I.norm())
    bound = ceil(4 * log(p))

    gcd_norm = coprime_factor * D**2

    # Stop infinite loops.
    for _ in range(1000):
        α = sum([b * randint(-bound, bound) for b in I.basis()])
        if gcd(ZZ(α.reduced_norm()), gcd_norm) == D:
            assert I == OI * α + OI * D
            return α
    raise ValueError(f"Cannot find a good α for D = {D}, I = {I}, n(I) = {D}")


def eichler_order_from_ideal(I):
    """
    The Eichler order is the intersection
    of two orders.

    Given an ideal I, we compute the Eichler
    order ℤ + I from the intersection of the
    left and right orders of I

    Proposition 1 (SQISign paper):

    EichlerOrder = O0 ∩ O = OL(I) ∩ OR(I) = ℤ + I.
    """
    return I.left_order().intersection(I.right_order())


# ========================================= #
#  Pushforward and Pullback of ideals to O0 #
# ========================================= #

def pullback_ideal(O0, O1, I, Iτ):
    """
    Input: Ideal I with left order O1
           Connecting ideal Iτ with left order O0
           and right order O1
    Output The ideal given by the pullback [Iτ]^* I
    """
    assert I.left_order() == O1
    assert Iτ.left_order() == O0
    assert Iτ.right_order() == O1

    N = ZZ(I.norm())
    Nτ = ZZ(Iτ.norm())

    α = ideal_generator(I)
    return O0 * N + O0 * α * Nτ


def pushforward_ideal(O0, O1, I, Iτ):
    """
    Input: Ideal I left order O0
           Connecting ideal Iτ with left order O0
           and right order O1
    Output The ideal given by the pushforward [Iτ]_* I
    """
    assert I.left_order() == O0
    assert Iτ.left_order() == O0
    assert Iτ.right_order() == O1

    N = ZZ(I.norm())
    Nτ = ZZ(Iτ.norm())

    K = I.intersection(O1 * Nτ)
    α = ideal_generator(K)
    return O1 * N + O1 * (α / Nτ)


# ================== #
#  Ideal Filtration  #
# ================== #

def ideal_filtration_steps(I, l, step_size, small_step_first=False):
    """
    Computes the step sizes for the IdealFiltration.

    When an ideal has a left order close to O0 then KLPT can fail as it's
    hard to find equivalent prime norm ideals.

    To remedy this, we can either set the last steps (which may be small)
    near the beginning or end of the filtration.
    """
    # Compute the exp in the ideal's norm
    v = ZZ(I.norm()).valuation(l)

    # We make the first step the shortest
    remaining_steps, small_step = divmod(v, step_size)

    # Edge case when all steps are the same size
    if small_step == 0:
        return [step_size] * remaining_steps

    # When the end is close to E0, we need to make the
    # small step first.
    if small_step_first:
        step_lengths = [small_step] + [step_size] * remaining_steps
    else:
        step_lengths = [step_size] * remaining_steps + [small_step]

    # Make sure we split this up OK
    assert sum(step_lengths) == v

    return step_lengths


def ideal_filtration(I, l, step_size, small_step_first=False):
    """
    Given an ideal I of norm l^* compute
    a chain, length e, of ideals I_i with norm step = l^f
    """
    # ideal_filtration expects a cyclic ideal
    assert is_cyclic(I), "Input ideal is not cyclic!"

    O = I.left_order()
    I_long = O.unit_ideal()
    I_short = O.unit_ideal()

    # Compute appropriate step lengths given I
    step_lengths = ideal_filtration_steps(
        I, l, step_size, small_step_first=small_step_first
    )

    # I = Ik ⊂ ... ⊂ I1 ⊂ I0
    # I_chain is the quotient of the contained ideals
    # I = Ik ⊂ ... ⊂ I1 ⊂ I0
    # with norm step_length

    # Build up a chain of the quotients
    I_chain = [I_short]
    for step_length in step_lengths:
        I_short_i = invert_ideal(I_long) * I + I_long.right_order() * l**step_length
        I_long = I_long * I_short_i
        I_chain.append(I_short_i)

    return I_chain

# ================================================= #
# Helper functions for tests, not used in main code #
# ================================================= #


def non_principal_ideal(O):
    """
    Computes a random non-principal ideal.
    Note: Only used for testing.
    """
    B = O.quaternion_algebra()
    p = -B.invariants()[1]

    # A slightly random ideal (not principal!)
    while True:
        beta = O.random_element(round(p))
        while not cornacchia_friendly(int(beta.reduced_norm())):
            beta = O.random_element(round(p))
        facNorm = factor(beta.reduced_norm())
        if facNorm[-1][1] > 1:
            continue
        N = facNorm[-1][0]
        I = O.left_ideal(
            [N * a for a in O.basis()] + [a * beta.conjugate() for a in O.basis()]
        )
        if I.conjugate().is_equivalent(O.unit_ideal()):
            continue
        return I