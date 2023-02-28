import unittest
from KLPT import RepresentIntegerHeuristic
from deuring import *
from setup import *

# Speedy and still (mostly) correct
proof.all(False)

class Testideal_to_kernel(unittest.TestCase):
    def helper(self, E, I, connecting_isogenies=None):
        print(f"Testing an ideal of norm n(I) = {factor(I.norm())}")

        D = ZZ(I.norm())
        K = ideal_to_kernel(E, I, connecting_isogenies=connecting_isogenies)
        self.assertEqual(K.order(), D), "Recovered kernel has the wrong order"

        ϕ = E.isogeny(K, algorithm="factored")
        self.assertEqual(ϕ.degree(), D), "ϕ has the wrong degree"

        J = kernel_to_ideal(K, D, connecting_isogenies=connecting_isogenies)
        self.assertEqual(J.norm(), D), "Recovered ideal has the wrong degree"
        self.assertEqual(I, J), "Ideals are not equal"

    def test_ideal_to_kernel(self):
        even_part = 2^10
        odd_part = T_prime

        γ = None
        while γ is None:
            γ = RepresentIntegerHeuristic(even_part * odd_part, parity=True)

        I_even = O0*γ + O0*even_part
        I_odd = O0*γ.conjugate() + O0*odd_part
        assert I_even.norm() == even_part
        assert I_odd.norm() == odd_part

        print(f"Simple test, using the curve E0")
        self.helper(E0, I_even)
        self.helper(E0, I_odd)

        print(f"Testing again, but using a curve != E0")
        P,Q = torsion_basis(E0, Dc)
        K = P + Q
        ϕ = E0.isogeny(K, algorithm="factored")
        E = ϕ.codomain()
        E.set_order((p^2 - 1)^2)
        ϕ_dual = dual_isogeny(ϕ, K, Dc)
        connecting_isogenies = (ϕ, ϕ_dual)

        self.helper(E, I_even, connecting_isogenies=connecting_isogenies)
        self.helper(E, I_odd, connecting_isogenies=connecting_isogenies)


class Testkernel_to_ideal(unittest.TestCase):
    def helper(self, K, D, connecting_isogenies=None):
        assert K.order() == D

        E = K.curve()

        print(f"Testing kernel of degree {factor(D)}")
        ϕ = E.isogeny(K, algorithm="factored")
        self.assertEqual(ϕ.degree(), D), "ϕ has the wrong degree"

        I = kernel_to_ideal(K, D, connecting_isogenies=connecting_isogenies)
        self.assertEqual(I.norm(), D), "Ideal has the wrong degree"

        _K = ideal_to_kernel(E, I, connecting_isogenies=connecting_isogenies)
        self.assertEqual(_K.order(), D), "The ideal derived has the wrong degree"

        _ϕ = E.isogeny(_K, algorithm="factored")
        self.assertEqual(_ϕ.degree(), D), "The isogeny created from the kernel created from the ideal is wrong"
        self.assertEqual(ϕ.codomain().j_invariant(), _ϕ.codomain().j_invariant()), "The recovered isogeny does not have the same codomain up to isomorphism..."
        self.assertEqual(ϕ.codomain(), _ϕ.codomain()), "The recovered isogeny does not have the same codomain up..."


    def test_kernel_to_idealE0(self):
        even_part = 2^10
        odd_part = T_prime

        P_even, Q_even = torsion_basis(E0, even_part)
        P_odd, Q_odd   = torsion_basis(E0, odd_part)

        self.helper(P_even + Q_even, even_part)
        self.helper(P_odd + Q_odd, odd_part)

    def test_kernel_to_ideal(self):
        # Same test, but with a connecting isogeny
        Pc,Qc = torsion_basis(E0, Dc)
        Kc = Pc + Qc
        ϕc = E0.isogeny(Kc, algorithm="factored")
        Ec = ϕc.codomain()
        Ec.set_order((p^2 - 1)^2)
        ϕc_dual = dual_isogeny(ϕc, Kc, Dc)
        connecting_isogenies = (ϕc, ϕc_dual)

        even_part = 2^10
        odd_part = T_prime

        P_even, Q_even = torsion_basis(Ec, even_part)
        P_odd, Q_odd   = torsion_basis(Ec, odd_part)

        self.helper(P_even + Q_even, even_part, connecting_isogenies=connecting_isogenies)
        self.helper(P_odd + Q_odd, odd_part, connecting_isogenies=connecting_isogenies)


if __name__ == '__main__' and '__file__' in globals():
    unittest.main()

