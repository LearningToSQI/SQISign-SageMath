"""
Functions which solve lattice problems for use in subalgorithms of KLPT.

For KLPT there are a few spots where we need to enumerate short vectors
of lattices to ensure the smallest possible solutions to Diophantine equations.

Namely, we need close vectors to a special lattice for the strong approximation
to ensure the output bound is ~pN^3. This is accomplished with GenerateCloseVectors.
We use FPYLLL for the underlying lattice computations, which seem to outperform
Pari. We also have the ability to enumerate rather than precompute all vectors, 
which is better than Pari's qfminim.

For the generation of equivalent prime norm ideals, we have an ideal basis and 
we find short norm vectors of this and immediately output algebra elements. 
There's probably ways to reuse the GenerateShortVectors, but there's a few 
things about the prime norm elements which require special treatment so we 
chose to suffer code duplication for clearer specific functions.
"""

# Sage Imports
from sage.all import vector, floor, ZZ, Matrix, randint

# fpylll imports
import fpylll
from fpylll import IntegerMatrix, CVP
from fpylll.fplll.gso import MatGSO


def solve_closest_vector_problem(lattice_basis, target):
    """
    Use the fpylll library to solve the CVP problem for a given
    lattice basis and target vector
    """
    L = IntegerMatrix.from_matrix(lattice_basis.LLL())
    v = CVP.closest_vector(L, target)
    # fpylll returns a type `tuple` object
    return vector(v)


def generate_short_vectors_fpyll(L, bound, count=2000):
    """
    Helper function for GenerateShortVectors and 
    generate_small_norm_quat which builds an iterator
    for short norm vectors of an LLL reduced lattice
    basis.
    """
    # # Move from Sage world to Fypll world
    A = IntegerMatrix.from_matrix(L)

    # Gram-Schmidt Othogonalization
    G = MatGSO(A)
    _ = G.update_gso()

    # Enumeration class on G with `count`` solutions
    # BEST_N_SOLUTIONS:
    # Starting with the nr_solutions-th solution, every time a new solution is found
    # the enumeration bound is updated to the length of the longest solution. If more
    # than nr_solutions were found, the longest is dropped.
    E = fpylll.Enumeration(
        G, nr_solutions=count, strategy=fpylll.EvaluatorStrategy.BEST_N_SOLUTIONS
    )

    # We need the row count when we call enumerate
    r = L.nrows()

    # If enumerate finds no solutions it raises an error, so we
    # wrap it in a try block
    try:
        # The arguments of enumerate are:
        # E.enumerate(first_row, last_row, max_dist, max_dist_expo)
        short_vectors = E.enumerate(0, r, bound, 0)
    except Exception as e:
        print(f"DEBUG [generate_short_vectors_fpyll]: No short vectors could be found...")
        print(f"{e}")
        short_vectors = []
        
    return short_vectors

def generate_short_vectors(lattice_basis, bound, count=2000):
    """
    Generate a generator of short vectors with norm <= `bound`
    returns at most `count` vectors.
    
    Most of the heavy lifting of this function is done by 
    generate_short_vectors_fpyll
    """
    L = lattice_basis.LLL()
    short_vectors = generate_short_vectors_fpyll(L, bound, count=count)
    for _, xis in short_vectors:
        # Returns values x1,x2,...xr such that
        # x0*row[0] + ... + xr*row[r] = short vector
        v3 = vector([ZZ(xi) for xi in xis])
        v = v3 * L
        yield v


def generate_close_vectors(lattice_basis, target, p, L, count=2000):
    """
    Generate a generator of vectors which are close, without
    bound determined by N to the `target`. The first
    element of the list is the solution of the CVP.
    """
    # Compute the closest element
    closest = solve_closest_vector_problem(lattice_basis, target)
    yield closest

    # Now use short vectors below a bound to find
    # close enough vectors

    # Set the distance
    diff = target - closest
    distance = diff.dot_product(diff)

    # Compute the bound from L
    b0 = L // p
    bound = floor((b0 + distance) + (2 * (b0 * distance).sqrt()))

    short_vectors = generate_short_vectors(lattice_basis, bound, count=count)

    for v in short_vectors:
        yield closest + v


def generate_small_norm_quat(Ibasis, bound, count=2000):
    """
    Given an ideal I and an upper bound for the scaled
    norm Nrd(a) / n(I), finds elements a ∈ B such that
    a has small norm.
    """
    # Before starting anything, just send out the basis
    # sometimes this works, and much quicker.
    for bi in Ibasis:
        yield bi

    # Recover Quaternion algebra from IBasis for use later
    B = Ibasis[0].parent()

    # Write Ibasis as a matrix
    Ibasis_matrix = Matrix([x.coefficient_tuple() for x in Ibasis]).transpose()

    # Can't do LLL in QQ, so we move to ZZ by clearing
    # the denominator
    lattice_basis, _ = Ibasis_matrix._clear_denom()
    L = lattice_basis.LLL()

    # Move from Sage world to Fypll world
    short_vectors = generate_short_vectors_fpyll(L, bound, count=count)

    for _, xis in short_vectors:
        # Returns values x1,x2,...xr such that
        # x0*row[0] + ... + xr*row[r] = short vector
        v3 = vector([ZZ(xi) for xi in xis])

        # Often the very shortest are all composite
        # this forces some of the latter element?
        # Decide with a coin-flip?
        if randint(0, 1) == 1:
            v3[2] = 1

        v = Ibasis_matrix * v3

        yield B(v)

    print(
        f"WARNING [generate_small_norm_quat]: "
        "Exhausted all short vectors, if you're seeing this SQISign is likely about to fail."
    )
