"""
Functions which implement the Deuring correspondence specialised for
SQISign.

The main functions which are used are:

EvalEndomorphism(): Given an alg element α ∈ B_{p, ∞} compute the action
                    α(P) for P ∈ E_0 using knowledge of mapping End(E0) and O0

IdealToKernel(): Given an ideal I, compute the kernel generator K ∈ E
                 such that ϕ_I : E / E ⟨K⟩. We follow ia.cr/2023/106
                 for a more efficient algorithm than presented in SQISign, but
                 include SQISign's impl. too.

IdealToIsogenyCoprime(): Given two equivalent ideals J, K with coprime norm and the
                         isogeny ϕK, compute ϕJ

IdealToIsogenyFromKLPT(): Given an ideal I with norm l^* and left order O, the connecting 
                          ideal K with norm l^*, left order O0 and right order O and the 
                          corresponding isogeny ϕK, find the isogeny ϕI.

                          Note: much of this function is handled iteratively by the 
                          function IdealToIsogenySmallFromKLPT() which does something
                          similar, but for input I with norm dividing the available 
                          torsion.

"""

# Sage imports
from sage.all import gcd, ZZ, factor, floor
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

# Local imports
from ideals import (
    equivalent_left_ideals,
    left_isomorphism,
    chi,
    chi_inverse,
    is_cyclic,
    is_integral,
    multiply_ideals,
    ideal_generator,
    quaternion_basis_gcd,
    ideal_filtration
)
from isogenies import (
    torsion_basis,
    dual_isogeny_and_kernel,
    dual_isogeny,
    EllipticCurveIsogenyFactored,
    BiDLP
)
from KLPT import EquivalentSmoothIdealHeuristic
from mitm import meet_in_the_middle_with_kernel
from utilities import has_order_D, print_info

from setup import E0, O0, l, T, f_step_max, Δ, p, sqrt_minus_one, prime_norm_heuristic

# ================================ #
#  Compute the action of End(E0)   #
# ================================ #

def E01(P):
    """
    Identity map, does nothing
    """
    return P


def E0ι(P):
    """
    Returns ι(P) = (-x, √-1 y)
    """
    if P.curve() != E0:
        raise ValueError("The endomorphism ι is defined on the curve E0")

    return E0(-P[0], sqrt_minus_one * P[1], P[2])


def E0π(P):
    """
    Returns π(P) = (X^p, Y^p, Z^p)
    """
    if P.curve() != E0:
        raise ValueError("The endomorphism π is defined on the curve E0")

    return E0(P[0] ** p, P[1] ** p, P[2] ** p)


def E0ιπ(P):
    """
    Returns ιπ(P) = (-X^p, √-1 Y^p, Z^p)
    """
    if P.curve() != E0:
        raise ValueError("The endomorphism ιπ is defined on the curve E0")

    return E0(-P[0] ** p, sqrt_minus_one * P[1] ** p, P[2] ** p)

# Store End(E0) action as an array
EndE0 = [E01, E0ι, E0π, E0ιπ]

# ===================== #
#  Evaluation of End(E) #
# ===================== #

def _check_connecting_isogenies(P, connecting_isogenies):
    """
    Helper function for eval_endomorphism()

    Test whether the curves and isogenies are as expected

    If P ∈ E != E0 then first we map from E → E0 with ϕ_dual.
    This is achieved by supplying the optional argument
    connecting_isogenies = (ϕ, ϕ_dual)

    ϕ : E0 → E, ϕ_dual : E → E0
    
    Note: this feature is unused for SQISign, but may be
    useful in other circumstances. 
    """
    # Initialise empty values to handle when connecting_isogenies
    # is None
    ϕ, ϕ_dual = None, None

    # Curve of the input point
    E = P.curve()

    # Curve has unknown End(E)
    if E != E0:
        if connecting_isogenies is None:
            raise ValueError(
                f"To work on a curve E != E0, a connecting isogeny ϕ : E0 -> E must be known."
            )

        # Check we have both the isogeny and its dual
        if len(connecting_isogenies) != 2:
            raise ValueError(
                "EvalEndomorphism requires both the connecting isogeny ϕ : E0 → E and its dual"
            )

        # Check the domain and codomains line up
        ϕ, ϕ_dual = connecting_isogenies
        if ϕ.domain() != ϕ_dual.codomain() or ϕ.codomain() != ϕ_dual.domain():
            raise ValueError(
                "The connecting isogeny ϕ : E0 → E is incompatible with supplied dual"
            )

        if ϕ.domain() != E0:
            raise ValueError(
                "The connecting isogeny must have to domain of the curve E0 with known End(E0)"
            )

        if ϕ.codomain() != E:
            raise ValueError(
                "The connecting isogeny must have to codomain of the supplied curve E"
            )

        # Now, map the point P so it's on the curve E0
        P = ϕ_dual(P)
    return P, ϕ

def eval_endomorphism(α, P, D, connecting_isogenies=None):
    """
    Evaluates the action of an endomorphism
    f ∈ End(E0) on a point P ∈ E.

    If E is not E0, this can still be done,
    but we need to know the connecting isogeny
    ϕ : E → E0.
    """
    # Verify connecting isogenies are correct, if present
    if connecting_isogenies:
        P, ϕ = _check_connecting_isogenies(P, connecting_isogenies)

    # Unpack the coefficients of the generator α, `d` is the lcm of the denominators
    # of the elements.
    d, *α_coeffs = α.denominator_and_integer_coefficient_tuple()

    # For integral ideals, we expect the denominator of elements to be at most 2
    assert d in (1, 2), "Something is wrong with the input ideal"
    if d == 2:
        # Divide out by two before evaluation if needed
        # TODO: we can avoid this with the Deuring friends paper trick
        P = P.division_points(d)[0]

    # Compute the image of α(P)
    P = sum(c * θ(P) for c, θ in zip(α_coeffs, EndE0))

    # If P ∈ E ≠ E0 then we need to map back
    # from E0 to E using the connecting isogeny
    if connecting_isogenies:
        P = ϕ(P)
    return P


# =============================== #
#    Ideal to Kernel Functions    #
# =============================== #

def derive_cyclic_generator(P, Q, D):
    """
    Given generators <P,Q> of a cyclic group
    of order D, find K such that G = <K>
    
    Heuristically, it seems easy to randomly
    find a K this way, and is about 10x faster
    than the deterministic method as we do not
    need to compute the order of P or Q.
    """
    K = P + Q
    for _ in range(1000):
        if has_order_D(K, D):
            return K
        K += Q
    raise ValueError(f"Never found a cyclic generator!")


def ideal_to_kernel(E, I, connecting_isogenies=None):
    """
    Given a supersingular elliptic curve E
    and a cyclic ideal I produces a generator
    P_I of E[I].

    Optional: If E is not E0, we can still
    find a generator provided a connecting
    isogeny ϕ : E → E0.

    Implementation follows ia.cr/2023/106
    which directly computes the kernel E[I] from the
    action of α_bar, rather than computing the kernel
    via discrete logs from the action of α.

    ker(ϕ) = ⟨α_bar(P), α_bar(Q)⟩ for E[n(I)] = ⟨P,Q⟩ 
    """
    assert is_cyclic(I), "Input ideal is not cyclic"

    # Degree of the isogeny we will to compute
    D = ZZ(I.norm())

    # Compute a generator such that I = O<α, D>
    α = ideal_generator(I)

    # Compute the torsion basis of E[D]
    P, Q = torsion_basis(E, D)

    # Evaluate R = α_bar(P)
    α_bar = α.conjugate()

    # If this has full order, we can stop here as R = α_bar(P)
    # generates the kernel
    R = eval_endomorphism(α_bar, P, D, connecting_isogenies=connecting_isogenies)
    if has_order_D(R, D):
        return R

    # Same again for S = α_bar(Q)
    S = eval_endomorphism(α_bar, Q, D, connecting_isogenies=connecting_isogenies)
    if has_order_D(S, D):
        return S

    # Neither R or S had full order, so we find a
    # linear combination of R, S which has order D
    return derive_cyclic_generator(R, S, D)


# ========================================= #
#     SQISign Ideal to Isogeny Functions    #
# ========================================= #

def IdealToIsogenyCoprime(J, K, ϕK):
    """
    Input:  Two equivalent left ideals J,K of O0
            where: J has norm dividing T^2
                    K has norm l^∙
            The isogeny ϕK : E0 → E0 / <K>

    Output  ϕJ : E0 → E0 / <J>
    """

    # Make sure the left orders are O0
    assert J.left_order() == O0
    assert K.left_order() == O0

    # Ensure the ϕK starts on E0
    assert ϕK.domain() == E0

    # Make sure the norms are correct
    nJ, nK = ZZ(J.norm()), ZZ(K.norm())
    assert gcd(nJ, nK) == 1
    assert nJ.divides(T**2)
    assert nK % l == 0

    # Assert the orders are equivalent
    assert equivalent_left_ideals(J, K)

    # Compute the element α
    α = chi_inverse(K, J)
    assert J == chi(α, K)

    # Compute the ideals Hi
    H1 = J + O0 * T
    H2 = O0 * α + O0 * (nJ / H1.norm())
    assert T**2 % H1.norm() == 0, "Norm of H1 does not divide T^2"
    assert T**2 % H2.norm() == 0, "Norm of H2 does not divide T^2"

    # Compute isogenies from Hi
    # ϕH1 : E0 → E1 = E0 / <H1>

    ϕH1_ker = ideal_to_kernel(E0, H1)
    ϕH1_ker_order = ZZ(H1.norm())
    ϕH1 = EllipticCurveIsogenyFactored(E0, ϕH1_ker, order=ϕH1_ker_order)
    E1 = ϕH1.codomain()
    E1.set_order((p**2 - 1)**2, num_checks=0)

    # We only need the kernel of ϕH2
    ϕH2_ker = ideal_to_kernel(E0, H2)

    # Construct EK, the codomain of ϕK
    EK = ϕK.codomain()
    EK.set_order((p**2 - 1)**2, num_checks=0)

    # ψ: EK → EK / ϕK (ker ϕH2)
    ψ_ker = ϕK(ϕH2_ker)
    ψ_ker_order = ZZ(H2.norm())
    ψ = EllipticCurveIsogenyFactored(EK, ψ_ker, order=ψ_ker_order)

    # Construct the curve Eψ which should be isomorphic to E1
    Eψ = ψ.codomain()
    Eψ.set_order((p**2 - 1)**2, num_checks=0)

    # Check Eψ is isomorphic to E1
    assert Eψ.is_isomorphic(E1)

    # Ensure the codomains match
    iso = Eψ.isomorphism_to(E1)
    ψ = iso * ψ
    ψ_dual = dual_isogeny(ψ, ψ_ker, order=H2.norm())

    # ψ_dual * ϕH1 : E0 → E1 → EK ≃ EJ = E0 / <J>
    ϕJ = ψ_dual * ϕH1

    assert ϕJ.domain() == ϕK.domain(), "ϕJ domain is wrong"
    assert ϕJ.codomain().is_isomorphic(ϕK.codomain()), "ϕJ codomain is wrong"
    assert ϕJ.degree() == nJ, "ϕJ degree is wrong"

    return ϕJ


# ====================================================== #
#  IdealToIsogenySmallFromKLPT. Warning: Big function!   #
# ====================================================== #

def IdealToIsogenySmallFromKLPT(I, J, K, ϕJ, ϕK, equivalent_prime_ideal=None):
    """
    Input: I a left O0-ideal of norm dividing T^2 l^(2f+Δ),
           an O0-ideal in J containing I of norm dividing T^2,
           and an ideal K ∼ J of norm a power of l
           The isogenies ϕJ = E0 / <J> and
                         ϕK = E0 / <K>

           Optional: I_prime allows us to speed up the KLPT
                     algorithm for deriving L by including an
                     ideal with small prime norm equivalent to
                     I

    Output: ϕ  = ϕ2 ◦ θ ◦ ϕ1 : E1 → E2 of degree l^(2f+Δ) such that 
            ϕI = ϕ ◦ ϕJ, L ∼ I of norm dividing T^2, ϕL = E / <L>.

    NOTE: There are some typos in the diagram / algorithm of the SQISign paper
          I hope these corrections help

    - Step 7, Nrd(γ) only has to divide T^2 l^(2f+Δ) n(K)
    - Step 9, using the Figure 4.1 we decompose ϕH2 as ψ2 ∘ ρ_dual_2 as the
      diagram shows ρ2 as the isogeny ρ2 : E5 → Eψ
    - Step 9, φ2 should be replaced with ρ_dual_2
    - Step 11, ψ1' is computed from [ρ_2 η]_* ψ1_dual, not [ϕ2_2 η]_* ψ1_dual


    We want to compute ϕ: E1 → E2
    We'll do this as:
        ϕ = ψ1' ∘ ρ2 ∘ η ∘ ψ1 ∘ ϕ1

    Figure 4 from SQISign

               ψ1'
    ┌────>Eψ ───────>E2
    │     ^          ^
    │     │ρ2        │ϕ2
    │     │          │
    │     E6         E4
    │     ^          ^
    │ψ2   ┊η         ┊θ
    │     ┊    ψ1    ┊
    │     E5<────────E3
    │                ^
    │                │ϕ1
    │       ϕJ       │
    E0══════════════>E1
            ϕK
    """

    # ==================================== #
    # Ensure that the input is as expected
    # ==================================== #

    # Check I is as expected
    assert I.left_order() == O0
    assert T**2 * l ** (2 * f_step_max + Δ) % I.norm() == 0

    # Check J, ϕJ are as expected
    assert J.left_order() == O0
    assert ZZ(J.norm()).divides(T**2)
    assert ϕJ.degree() == J.norm()
    assert ϕJ.domain() == E0

    # Check K, ϕK are as expected
    assert K.left_order() == O0
    assert K.norm() % l == 0 or K.norm() == 1
    assert ϕK.degree() == K.norm()
    assert ϕK.domain() == E0

    # Make sure ϕJ and ϕK start and end up on isomorphic curves
    assert ϕJ.domain().is_isomorphic(ϕK.domain())
    assert ϕJ.codomain().is_isomorphic(ϕK.codomain())

    # ==================================================== #
    #   Helper Functions for IdealToIsogenySmallFromKLPT   #
    # ==================================================== #
    r"""
    derive_ϕ1(): Used to compute ϕ1 from I
                 the isogeny ϕ1 : E1 → E3

    derive_L_and_gamma(): Given the ideals I,K,L
                 compute the ideal L and γ. Runs
                 until gcd(γ) = 2^k for k >= 0

    derive_isogenies_from_H1(),
    derive_isogenies_from_H2():
                 Given the ideals H1,H2 find the
                 (factored) isogenies with Hi as
                 the isogenies kernel.

    derive_final_pushforwards():
                 Given ψ1, ρ2_dual, and, η compute
                 ψ1_prime, λ1 and λ3 used to derive
                 the final isogenies.
    """

    def derive_ϕ1(I, ϕJ, E1, step_size):
        """
        Given an ideal I, compute an ideal I1 of
        norm l^f and then using the pushforward
        from ϕJ, compute a degree l^f isogeny 
        ϕ1: E1 → E3

        """
        # Compute the ideal I1
        I1 = I + O0 * l**step_size
        assert l**step_size % I1.norm() == 0

        # Compute the isogeny from ϕ1': E0 -> E0 / <I1>
        ϕ1_prime_ker = ideal_to_kernel(E0, I1)
        assert (
            ϕ1_prime_ker.order() == l**step_size
        ), "The degree of the kernel ϕ1_prime_ker is incorrect"

        # Now we compute the push forward to compute
        # ϕ1 : E1 → E3
        #
        # Compute the isogeny ϕ1 : E1 → E3
        ϕ1_ker = ϕJ(ϕ1_prime_ker)
        ϕ1_ker_order = ZZ(I1.norm())
        ϕ1 = EllipticCurveIsogenyFactored(E1, ϕ1_ker, order=ϕ1_ker_order)
        E3 = ϕ1.codomain()
        E3.set_order((p**2 - 1)**2, num_checks=0)
        assert (
            ϕ1.degree().divides(l**f_step)
        ), f"Degree of {factor(ϕ1.degree()) = } does not divide {l^f_step}"
        return ϕ1, E1, E3

    def derive_L_and_gamma(I, J, K, step_size, equivalent_prime_ideal=None):
        """
        Given ideals I,J,K find a prime norm ideal L equivalent to 
        I and then compute γ, used to generate the ideals H1 and H2

        Optional: equivalent_prime_ideal allows us to supply an ideal
        I' equivalent to I with prime norm.
        """
        # g keeps track of backtracking in γ
        g = -1
        while g != 1 and not (len(factor(g)) == 1 and g % l == 0):

            # If the ideal K has small l-valuation, we can send in M as a
            # power of two to KLPT and skip prime ideal generation
            if l ** (step_size + ZZ(K.norm()).valuation(l)) < prime_norm_heuristic:
                α = left_isomorphism(J, K)
                assert J * α == K

                # Send in an ideal with norm a power of two
                M = I * α
                L = EquivalentSmoothIdealHeuristic(
                    M,
                    T**2,
                    equivalent_prime_ideal=equivalent_prime_ideal,
                    near_prime=True,
                )

            else:
                L = EquivalentSmoothIdealHeuristic(
                    I,
                    T**2,
                    equivalent_prime_ideal=equivalent_prime_ideal,
                    near_prime=True,
                )
            if L is None:
                print(
                    f"DEBUG: [IdealToIsogenySmallFromKLPT]"
                    "EquivalentSmoothIdealHeuristic failed... trying again"
                )
                continue
            nL = ZZ(L.norm())
            assert nL.divides(T**2), "The norm of the ideal L does not divide T^2"

            # Compute the elements α,β

            # If K has norm 1 then ChiInverse(K, J)
            # will return i resulting in the kernel being
            # twisted by the automorphism E0ι.
            # To avoid this, manually set α = 1
            if K.norm() == 1:
                α = 1
            else:
                α = chi_inverse(K, J)
            β = chi_inverse(I, L)

            # Check we have the correct elements α, β
            assert J == chi(α, K) and α in K
            assert L == chi(β, I) and β in I

            # Compute gamma and check its properties
            nJ = J.norm()
            γ = (β * α) / nJ
            g = quaternion_basis_gcd(γ, O0)

            if g != 1:
                # g is a power of l, so we can easily correct
                # for backtracking by extending the mitm
                if g % l == 0 and len(factor(g)) == 1:
                    # Make gamma primitive
                    γ = γ / g
                else:
                    # The gcd is bad, so we need to compute L again
                    print(
                        (f"DEBUG: [IdealToIsogenySmallFromKLPT]:"
                         "gcd(γ) = {g}, edge case not currently supported, generating L again...")
                    )
                    continue

            return L, nL, γ, g

    def derive_isogenies_from_H1(H1_odd, ϕ1, ϕK, E3):
        """
        ϕH1 : E0 -> E5
        ϕH1 = ψ1 ◦ ϕ1 ◦ ϕK where ψ1 has degree T
                                 ϕK has degree nK
                                 ϕ1 has degree l^f_step

        ϕH1 = ϕodd ∘ ϕeven where ϕodd has degree T
                                 ϕeven has degree nK*l^f_step

        Our goal is to compute ψ1

                   ϕeven
              E0 ─────────> Eeven
              │             │
              │             │
        ϕ1∘ϕk │             │ ϕodd
              │             │
              │             │
              v             v
              E3 ─────────> E5
                    ψ1

        """
        print(f"INFO [IdealToIsogenySmallFromKLPT]: Computing isogenies from H1...")

        # First compute the kernel from the ideal
        ϕH1_ker_odd = ideal_to_kernel(E0, H1_odd)

        # Compute the pushforward
        ψ1_ker = ϕ1(ϕK(ϕH1_ker_odd))
        ψ1_order = ZZ(H1_odd.norm())

        ψ1 = EllipticCurveIsogenyFactored(E3, ψ1_ker, order=ψ1_order)
        E5 = ψ1.codomain()
        E5.set_order((p**2 - 1)**2, num_checks=0)
        return ψ1, ψ1_ker, ψ1_order, E5

    def derive_isogenies_from_H2(H2):
        """
        Want to derive isogenies ψ2 : E0 → Eψ
        and ρ̂2  : Eψ → E6 from ϕH2 : E0 -> E6 where
        ϕH2 = ρ̂2 ◦ ψ2 where ψ2 has degree dividing T
                            ρ̂2 has degree 2^f

        """
        print(f"INFO [IdealToIsogenySmallFromKLPT]: Computing isogenies from H2...")
        # First compute the kernel from H2
        ϕH2_ker = ideal_to_kernel(E0, H2)
        H2_l_val = (H2.norm()).valuation(l)

        # Compute even and odd kernel orders
        ρ2_dual_order = l**H2_l_val
        ψ2_order = H2.norm() // ρ2_dual_order

        # Find the subgroup of order (dividing) T within ϕH2_ker
        ψ2_ker = ρ2_dual_order * ϕH2_ker

        # Compute the isogeny ψ2 and the curve Eψ
        # Note: Eψ is not named in the SQISign paper
        ψ2 = EllipticCurveIsogenyFactored(E0, ψ2_ker, order=ψ2_order)
        Eψ = ψ2.codomain()
        Eψ.set_order((p**2 - 1)**2, num_checks=0)

        # Compute the subgroup of order 2^f
        # on the curve Eψ by pushing through
        # a point of order 2^f on E0
        ρ2_dual_ker = ψ2_order * ϕH2_ker
        ρ2_dual_ker = ψ2(ρ2_dual_ker)

        # Compute the isogeny ρ̂2 : Eψ -> E6
        ρ2_dual = EllipticCurveIsogenyFactored(Eψ, ρ2_dual_ker, order=ρ2_dual_order)
        E6 = ρ2_dual.codomain()
        E6.set_order((p**2 - 1)**2, num_checks=0)

        # Check the end points all match up
        assert ψ2.domain() == E0
        assert ψ2.codomain() == ρ2_dual.domain()

        return ψ2, Eψ, ρ2_dual, ρ2_dual_ker, ρ2_dual_order, E6

    def derive_final_pushforwards(
        ψ1, ψ1_ker, ψ1_order, ρ2_dual, ρ2_dual_ker, ρ2_dual_order, η, η_ker, η_order
    ):
        """
        Compute the pushforwards

        We derive two isogenies from the following
        isogeny square

        ψ1' = [ρ2 ∘ η]_* ψ1_dual
            = (ρ2 ∘ η) ker(ψ1_dual)

        ϕ2 ∘ θ = [ψ1_dual]_* ρ2 ∘ η
               = ψ1_dual(ker(ρ2 ∘ η))

        We do not know ker(ρ2 ∘ η) so instead

        ϕ2 ∘ θ = λ3 * λ1

        Where:

        λ1 : E3 -> Eλ has kernel ψ1_dual(ker(η))
        λ2 : E6 -> Eλ has kernel η(ker(ψ1_dual))
        λ3 : Eλ -> E2 has kernel λ2(ker(ρ2))

                   ψ1'
            Eψ ───────────>E2<────┐
            ^              ^      │
            │              │      │
            │              │      │
        ρ2  │              │ λ3   │ ϕ2
            │              │      │
            │      λ2      │      │
            E6 ───────────>Eλ     E4
            ^              ^      ^
            ┊              │      ┊
         η  ┊              │ λ1   ┊ θ
            ┊              │      ┊
            ┊              │      ┊
            E5 ──────────> E3 ╌╌╌╌┘
                ψ1_dual
        """

        # Compute the dual of ψ1, ψ2, ρ2 together with a point generating
        # the dual isogeny's kernel
        print(
            f"INFO [IdealToIsogenySmallFromKLPT]: Computing the duals of ψ1 and ρ2_dual"
        )
        ψ1_dual, ψ1_dual_ker = dual_isogeny_and_kernel(ψ1, ψ1_ker, order=ψ1_order)
        ρ2, ρ2_ker = dual_isogeny_and_kernel(ρ2_dual, ρ2_dual_ker, order=ρ2_dual_order)

        # As we know ψ1_dual_ker, this is easy
        # ψ1' = [ρ2 ∘ η]_* ψ1_dual
        print(f"INFO [IdealToIsogenySmallFromKLPT]: Computing the isogeny ψ1_prime")
        ψ1_prime_ker = ρ2(η(ψ1_dual_ker))
        ψ1_prime = EllipticCurveIsogenyFactored(Eψ, ψ1_prime_ker, order=ψ1_order)

        # ϕ2 ∘ θ = [ψ1_dual]_* ρ2 ∘ η
        # We do not know the kernel of
        # (ρ2 ∘ η), so here's a work around

        # λ1 : E3 → Eλ
        λ1_ker = ψ1_dual(η_ker)
        print(f"INFO [IdealToIsogenySmallFromKLPT]: Computing the isogeny λ1")
        λ1 = EllipticCurveIsogenyFactored(E3, λ1_ker, order=η_order)
        Eλ = λ1.codomain()
        Eλ.set_order((p**2 - 1)**2, num_checks=0)

        # λ2 : E6 → Eλ
        λ2_ker = η(ψ1_dual_ker)
        print(f"INFO [IdealToIsogenySmallFromKLPT]: Computing the isogeny λ2")
        λ2 = EllipticCurveIsogenyFactored(E6, λ2_ker, order=ψ1_order)
        Eλ2 = λ2.codomain()
        iso = Eλ2.isomorphism_to(Eλ)
        λ2 = iso * λ2
        assert λ2.codomain() == Eλ

        # λ3 : Eλ → E2
        λ3_ker = λ2(ρ2_ker)
        print(f"INFO [IdealToIsogenySmallFromKLPT]: Computing the isogeny λ3")
        λ3 = EllipticCurveIsogenyFactored(Eλ, λ3_ker, order=ρ2_dual_order)

        return ψ1_prime, λ1, λ3

    # ======================================================== #
    #   End Helper Functions for IdealToIsogenySmallFromKLPT   #
    # ======================================================== #

    # Set ϕK to have the same codomain as ϕJ
    iso = ϕK.codomain().isomorphism_to(ϕJ.codomain())
    ϕK = iso * ϕK
    assert ϕK.codomain() == ϕJ.codomain()

    # Accounts for last step where Δ may be smaller
    nI = I.norm()
    step_size = nI.valuation(l)

    f_step = min(f_step_max, floor(step_size / 2))
    Δ_actual = max(step_size - 2 * f_step, 0)
    assert 2 * f_step + Δ_actual == step_size

    # Norms will be useful later
    nJ = ZZ(J.norm())
    nK = ZZ(K.norm())

    # First, find the domain of the isogeny ϕ1
    E1 = ϕJ.codomain()
    E1.set_order((p**2 - 1)**2, num_checks=0)
    assert E1.is_isomorphic(ϕK.codomain()), "ϕJ and ϕK do not end on the same curve!"

    # When the step size is small enough, we can skip the complicated
    # steps and directly derive the needed isogeny. This assumes that
    # this is performed as the last step in `IdealToIsogenyFromKLPT()`
    if step_size < f_step:
        ϕ1, _, _ = derive_ϕ1(I, ϕJ, E1, step_size)
        L = EquivalentSmoothIdealHeuristic(I, T**2)
        nL = ZZ(L.norm())
        assert T**2 % nL == 0, "The norm of the ideal L does not divide T^2"

        # Early return
        return ϕ1, L, None

    # Derive ϕ1
    ϕ1, E1, E3 = derive_ϕ1(I, ϕJ, E1, f_step)

    # To continue, we first need to do some KLPT magic
    # Compute the ideal L equivalent to I with norm dividing T^2
    L, nL, γ, g = derive_L_and_gamma(
        I, J, K, step_size, equivalent_prime_ideal=equivalent_prime_ideal
    )

    # Check γ is in the correct ideals
    assert g * γ in K and g * γ.conjugate() in L
    # Check γ has correct reduced norm
    assert γ.reduced_norm() == (nI * nL * nK // g**2) / (nJ)

    # Check γ has reduced norm that divides (nI T^2 nK) / nJ
    assert T**2 * l ** (2 * f_step + Δ_actual) * nK % γ.reduced_norm() == 0
    
    # Compute the ideals H1, H2
    # TODO: we can remove this, but we use
    # n(H1) when computing H2
    H1 = O0 * γ + O0 * (nK * l**f_step * T)

    # We only will need the odd part to compute ψ1
    H1_odd = O0 * γ + O0 * T

    # Note:
    # The algorithm in the paper states:
    # H2 = O0*γ.conjugate() + O0*l^f_step*T
    # But this is only true when
    # Nrd(γ) = T^2*l^(2*f_step+Δ)*nK
    # As γ will only divide this, we use the following:
    H2 = O0 * γ.conjugate() + O0 * (γ.reduced_norm() / ((l**Δ_actual * H1.norm())))

    assert is_cyclic(H1), "H1 is not cyclic"
    # TODO: this bug sometimes appears. What is the cause?
    if not is_cyclic(H2):
        print(f"{is_integral(H1) = }")
        print(f"{is_integral(H2) = }")
        print(f"{factor(H1.norm()) = }")
        print(f"{factor(H2.norm()) = }")
        print(f"{factor(γ.conjugate().reduced_norm()) = }")

    ψ1, ψ1_ker, ψ1_order, E5 = derive_isogenies_from_H1(H1_odd, ϕ1, ϕK, E3)
    ψ2, Eψ, ρ2_dual, ρ2_dual_ker, ρ2_dual_order, E6 = derive_isogenies_from_H2(H2)

    # =============================== #
    #   Meet in the middle step (😬)  #
    # =============================== #

    # We expect there to be an isogeny degree l^Δ linking these.
    # However, if we had a non-trivial gcd, we have a little
    # further to brute-force
    gap = Δ_actual + 2 * ZZ(g.valuation(l))
    η_order = l**gap
    print(
        f"DEBUG [IdealToIsogenySmallFromKLPT]: Attempting a mitm with: {factor(η_order)}"
    )
    η, η_ker = meet_in_the_middle_with_kernel(E5, E6, l, gap)
    print(
        f"INFO [IdealToIsogenySmallFromKLPT]: Found the meet in the middle isogeny"
    )

    # ================================= #
    #   Compute the final pushforwards  #
    # ================================= #

    ψ1_prime, λ1, λ3 = derive_final_pushforwards(
        ψ1, ψ1_ker, ψ1_order, ρ2_dual, ρ2_dual_ker, ρ2_dual_order, η, η_ker, η_order
    )

    # Compute ϕ2θ : E3 → E2
    ϕ2θ = λ3 * λ1

    # ϕ : E1 → E2
    ϕ = ϕ2θ * ϕ1

    # ψ : E0 → E2
    ψ = ψ1_prime * ψ2

    # Do the curves start/end in the right spot
    assert ϕ.codomain().is_isomorphic(
        ψ.codomain()
    ), "ϕ and ψ do not end on isomorphic curves"
    assert ϕ.domain().is_isomorphic(ϕ1.domain()), "ϕ does not start on E1"
    assert ψ.domain().is_isomorphic(E0), "ψ does not start on E0"

    # Do the degrees / norms line up
    assert ψ.degree() == nL, "degree of ψ is incorrect"

    return ϕ, L, ψ

# ======================= #
#  IdealToIsogenyFromKLPT #
# ======================= #

def IdealToIsogenyFromKLPT(I, K, ϕK, I_prime=None, K_prime=None, end_close_to_E0=False):
    """
    Computes the isogeny ϕI whose kernel corresponds to the ideal I.

    Input: A left O-ideal I of norm a power of l,
           K a left O0-ideal and right O-ideal of norm l^•
           The corresponding ϕK : E / <K>.

           Optional:
           I_prime: an ideal equivalent to I with small prime norm
           K_prime: an ideal equivalent to K with small prime norm

           Explanation:
           If we know an ideal with small prime norm which is equivalent
           to I or K we can speed up this algorithm by skipping part of
           the KLPT step inside IdealToIsogenySmallFromKLPT or
           IdealToIsogenyCoprime respectively.

    Output: ϕI : E / <I>
    """
    # Ensure the norms are as expected
    assert I.norm() % l == 0 or I.norm() == 1
    assert K.norm() % l == 0 or K.norm() == 1

    # Ensure the orders are as expected
    assert I.left_order() == K.right_order()
    assert K.left_order() == O0

    # If we supply equivalent_prime_ideal, make sure it
    # is of the right form
    if I_prime:
        assert equivalent_left_ideals(
            I, I_prime
        ), "Input I_prime is not equivalent to I"
        assert ZZ(I_prime.norm()).is_prime(), "Input I_prime does not have prime order"

    if K_prime:
        assert equivalent_left_ideals(
            K, K_prime
        ), "Input K_prime is not equivalent to I"
        assert ZZ(K_prime.norm()).is_prime(), "Input K_prime does not have prime order"

    # =============================================== #
    #   Helper Functions for IdealToIsogenyFromKLPT   #
    # =============================================== #

    def derive_J_and_phi_J(K, ϕK, K_prime=None):
        """
        Given a connecting ideal K and corresponding isogeny
        compute an equivalent ideal J with norm coprime to
        K and the equivalent isogeny ϕJ

        Optional: K_prime is an equivalent ideal to K with 
        prime norm
        """
        # In keygen we send
        #  - K  = O0.unit_ideal()
        #  - ϕK = E0.isogeny(E0(0))
        # Which allows us to set 'J' to be the unit 
        # ideal, and ϕJ to be a trivial isogeny too
        if K.norm() == 1:
            J = O0.unit_ideal()
            ϕJ = E0.isogeny(E0(0))
            return J, ϕJ

        # Sometimes we already know an equivalent prime
        # norm ideal, so we use this and save some time
        # within KLPT
        if K_prime:
            # TODO: is there any point in lopping this??
            for _ in range(10):
                J = EquivalentSmoothIdealHeuristic(
                    K, T**2, equivalent_prime_ideal=K_prime
                )
                if J:
                    break
            ϕJ = IdealToIsogenyCoprime(J, K, ϕK)
            return J, ϕJ

        # Generic case
        # Compute a smooth norm ideal J
        J = None
        for _ in range(10):
            J = EquivalentSmoothIdealHeuristic(K, T**2)
            if J is not None:
                break

        if J is None:
            exit("Was unable to compute an equivalent ideal J ~ K")

        # Compute the isogeny ϕJ : E0 / <J>
        ϕJ = IdealToIsogenyCoprime(J, K, ϕK)

        assert equivalent_left_ideals(J, K)
        return J, ϕJ

    # ==================== #
    # End Helper Functions #
    # ==================== #

    J, ϕJ = derive_J_and_phi_J(K, ϕK, K_prime=K_prime)

    # Compute a chain:
    # I = Iv ⊂ ... ⊂ I1 ⊂ I0
    # I_short is the quotient of these elements
    I_short = ideal_filtration(
        I, l, (2 * f_step_max + Δ), small_step_first=end_close_to_E0
    )

    # Create a list to store the isogeny factors
    ϕi_factors = []

    # For the last step, we can use JIi_prime = I_prime
    # To avoid problems computing equivalent prime norm ideals
    JIi_prime = None

    # First element in Iv is the unit ideal of O
    # We don't need this.
    for ind, Ii in enumerate(I_short[1:], start=1):
        print_info(
            f"STARTING WITH ELEMENT {ind}/{len(I_short) - 1} OF FILTRATION CHAIN...", 
            banner="-"
        )
        # On last loop, use the trick that we know an equivalent
        # prime norm ideal
        if ind == len(I_short) - 1:
            JIi_prime = I_prime

        # J * Ii
        alpha = left_isomorphism(K, J)  # K*α = J
        JIi = multiply_ideals(J, Ii, beta=alpha)

        # Compute
        # ϕ = ϕ2 ◦ θ ◦ ϕ1 : E1 → E2 of degree l^(2f+Δ) such that ϕJIi = ϕ ◦ ϕJ,
        # J ∼ JIi of norm dividing T^2,
        # ϕJ = E / <J> (output J, not input J)

        ϕi, J, ϕJ = IdealToIsogenySmallFromKLPT(
            JIi, J, K, ϕJ, ϕK, equivalent_prime_ideal=JIi_prime
        )

        # Add the isogeny ϕi to the list of factors
        # May need to correct it with an isomorphism
        if ind > 1:
            ϕi_prev = ϕi_factors[ind - 2]
            ι = ((ϕi_prev).codomain()).isomorphism_to(ϕi.domain())
            ϕi = ϕi * ι
        ϕi_factors.append(ϕi)

        # If we're not in the last loop
        if ind != len(I_short) - 1:
            # Update the ideal K and the isogeny ϕK
            K = K * Ii
            ϕK = ϕi * ϕK

    ϕI = EllipticCurveHom_composite.from_factors(ϕi_factors)

    return ϕI


# ============================================ #
#   Functions for Kernel to Isogeny Functions  #
# ============================================ #


def compute_coprime_basis(D):
    """
    Start with basis <1, i, (i + j) / 2, (1 + k) / 2>
    and find a new basis such that the norm of each basis
    element is coprime to `D`.

    TODO: is this the best method?
    """
    O0_basis = O0.basis()
    θs = []
    for f in O0_basis:
        while True:
            if gcd(f.reduced_norm(), D) == 1:
                θs.append(f)
                break
            f += O0_basis[0] + O0_basis[1]
    return θs


def find_torsion_basis_EndE(θPs, D):
    """
    Looks for θi, θj such that θi(P), θj(P) generates E[D]
    """
    for i in range(4):
        for j in range(i + 1, 4):
            eθiθj = θPs[i].weil_pairing(θPs[j], D, algorithm="pari")
            if has_order_D(eθiθj, D, multiplicative=True):
                return i, j
    raise ValueError(f"No basis for E[D] found with given point")


def kernel_to_ideal(P, D, connecting_isogenies=None):
    """
    Given a point P ∈ E[D] compute the
    ideal I(<P>))

    Optional: If E is not E0, we can still
    find an ideal provided a connecting
    isogeny ϕ : E → E0.
    """
    # Compute a basis β1,β2,β3,β4 of O0
    # with norm coprime to D
    βs = compute_coprime_basis(D)

    # Compute the image of all the points
    # β(P) by acting with θ ≅ β
    θs = [eval_endomorphism(β, P, D, connecting_isogenies=connecting_isogenies) for β in βs]

    # Find θi, θj which generates E[D]
    i, j = find_torsion_basis_EndE(θs, D)
    θi, θj = θs[i], θs[j]

    # Pick k ≠ i,j such that
    k = set([0, 1, 2, 3]).difference([i, j]).pop()
    θk = θs[k]

    # Solve the discrete log
    a, b = BiDLP(θk, θi, θj, D)
    assert a * θi + b * θj == θk

    # Create the Quaternion Algebra element
    α = βs[k] - a * βs[i] - b * βs[j]
    return O0 * α + O0 * D
