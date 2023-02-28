import unittest
from KLPT import *
from ideals import equivalent_left_ideals, non_principal_ideal
from setup import *

I = non_principal_ideal(O0)

def CyclicIdeal(l, e):
    """
    Finds a cyclic Ideal I with norm l^e

    Note: I have no idea if this is a good method,
    I just grabbed it / adapted it from:
    https://hxp.io/blog/34/hxp-CTF-2017-crypto600-categorical-writeup/
    """
    while True:
        γ = RepresentIntegerHeuristic(l**e)
        I = O0.unit_ideal().scale(γ)
        if is_cyclic(I):
            return I

class TestKLPT(unittest.TestCase):
    @staticmethod
    def compute_J_and_gamma():
        """
        Helper for tests
        """
        previously_seen = set()
        facT = list(factor(T))

        for _ in range(100):
            J, N, _ = EquivalentPrimeIdealHeuristic(I, previous=previously_seen)
            previously_seen.add(N)
            if not J: continue
            if kronecker(l, N) == -1:
                ex, _ = FindExtra(N, facT.copy())
                γ = RepresentIntegerHeuristic(N*ex)
                if γ is not None:
                    return J, N, γ
        raise ValueError("Never found a prime norm with a suitable γ...")

    def test_RepresentIntegerHeuristic(self):
        γ = None
        for _ in range(100):
            M = randint(p+2, p^2)
            γ = RepresentIntegerHeuristic(M)
            if γ: break

        if not γ:
            print(f"DEBUG: [test_RepresentIntegerHeuristic] Could not find a γ.")
            self.assertTrue(False)
        # Found a gamma, make sure it works!
        self.assertEqual(γ.reduced_norm(), M)

    def test_EquivalentPrimeIdealHeuristic(self):
        J = None
        for _ in range(100):
            J, N, α = EquivalentPrimeIdealHeuristic(I)
            if J: break

        if not J:
            print(f"DEBUG: [test_EquivalentPrimeIdealHeuristic] Could not find an ideal")
            self.assertTrue(False)

        # Found a prime ideal, make sure it has the right
        # properties. 
        self.assertTrue(Chi(α, I) == J)
        self.assertTrue(ZZ(N).is_prime())
        self.assertTrue(equivalent_left_ideals(I, J))

    def test_IdealModConstraint(self):
        J, _, γ = self.compute_J_and_gamma()
        C0, D0 = IdealModConstraint(J, γ)
        # Ensure that γ*μ0 is an element of J
        μ0 = j*(C0 + ω*D0)
        self.assertTrue(γ*μ0, J)

    def test_StrongApproximationHeuristic(self):
        # Ensure that for a given l,
        # that l is a non-quadratic residue
        # mod N
        J, N, γ = self.compute_J_and_gamma()
        C0, D0 = IdealModConstraint(J, γ)
        ν = StrongApproximationHeuristic(N, C0, D0, factor(l^e))
        if ν is None:
            print("Heuristic algorithm `StrongApproximationHeuristic` failed")
            self.assertTrue(False)
        ν_norm = ZZ(ν.reduced_norm())
        self.assertTrue(ν_norm % l == 0)

    def test_EquivalentSmoothIdealHeuristic_le(self):
        norm_l = l^e
        J = EquivalentSmoothIdealHeuristic(I, norm_l)
        self.assertTrue(equivalent_left_ideals(I,J))
        self.assertTrue(norm_l % ZZ(J.norm()) == 0)
        print(f"Norm for l^e: {factor(J.norm())}")

    def test_EquivalentSmoothIdealHeuristic_T(self):
        norm_T = T^2
        J = EquivalentSmoothIdealHeuristic(I, norm_T)
        print(f"Norm: {factor(J.norm())}")
        self.assertTrue(equivalent_left_ideals(I,J))
        self.assertTrue(norm_T % ZZ(J.norm()) == 0)

    def test_EquivalentSmoothIdealHeuristic_small_norm(self):
        norm_T = T^2
        I = CyclicIdeal(2, 80)
        I_small = I + O0*2^10
        assert is_cyclic(I_small)

        J = EquivalentSmoothIdealHeuristic(I_small, norm_T)
        print(f"Norm: {factor(J.norm())}")
        self.assertTrue(equivalent_left_ideals(I_small,J))
        self.assertTrue(norm_T % ZZ(J.norm()) == 0)

if __name__ == '__main__' and '__file__' in globals():
    unittest.main()