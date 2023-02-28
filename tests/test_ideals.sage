import unittest
from ideals import *
from setup import *

class TestIdealHelpers(unittest.TestCase):
    def test_equivalent_right_ideals(self):
        I1 = O0.right_ideal([b for b in O0.basis()])
        I2 = B.random_element() * I1
        self.assertTrue(equivalent_right_ideals(I1, I2))

    def test_equivalent_left_ideals(self):
        I1 = O0.left_ideal([b for b in O0.basis()])
        I2 = I1 * B.random_element()
        self.assertTrue(equivalent_left_ideals(I1, I2))

    def test_LeftIsomorphism(self):
        I = non_principal_ideal(O0)
        for _ in range(10):
            β = B.random_element()
            J = I*β
            α = left_isomorphism(I,J)
            self.assertEqual(J, I*α)

if __name__ == '__main__' and '__file__' in globals():
    unittest.main()