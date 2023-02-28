import unittest
from lattices import *

class TestCloseVectors(unittest.TestCase):
    def test_SolveCVP(self):
        """
        Solution found from
        IntegerLattice(A).closest_vector()

        Which does what we need, but is slower.
        SolveCVP uses fpylll and is about 35x faster
        (2ms compared to 70ms)
        """
        A = [(-8, 0, 0, -1), (10, -1, 0, 1), (14, 1, 2, -1), (-2, -1, 2, 2)]
        A = Matrix(ZZ, A)
        t = vector([1,2,3,4])

        c = vector([0, 2, 2, 3])
        self.assertEqual(SolveCVP(A,t), c)

    def test_generate_close_vectors(self):
        # Stupid example
        A = diagonal_matrix(ZZ, 4, [1,2,3,4])
        assert A.is_positive_definite()
        bound = 100
        short_vectors = generate_close_vectors(A, bound)
        self.assertTrue(all([v.dot_product(v) <= bound for v in short_vectors]))

if __name__ == '__main__' and '__file__' in globals():
    unittest.main()