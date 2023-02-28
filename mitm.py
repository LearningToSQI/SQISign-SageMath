"""
Algorithms to derive an unknown isogeny between two elliptic curves
of known degree 2^k. 

This is implemented by computing two isogeny graphs from each 
elliptic curve and looking for a collision in the leaves of the 
respective tree graphs.

The graph is efficiently built by deriving j-invariants from the
roots of modular polynomials and once a path through the graph is
known, the isogeny is computed by checking all degree 2 isogenies
from each graph.
"""

# Sage imports
from sage.all import (
    floor,
    PolynomialRing,
    factor,
    EllipticCurveIsogeny,
    EllipticCurve,
)
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

# Local imports
from isogenies import (
    generate_kernels_division_polynomial,
    kernel_from_isogeny_prime_power,
    generate_kernels_prime_power,
)
from setup import Fp2

# ============================================= #
# Compute j-invariants from modular polynomials #
# ============================================= #

# For faster quadratic root computation
Fp2_inv_2 = Fp2(1) / 2

def sqrt_Fp2(a):
    """
    Efficiently computes the sqrt
    of an element in Fp2 using that
    we always have a prime p such that
    p ≡ 3 mod 4.
    """
    p = Fp2.characteristic()
    i = Fp2.gens()[0]  # i = √-1

    a1 = a ** ((p - 3) // 4)
    x0 = a1 * a
    α = a1 * x0

    if α == -1:
        x = i * x0
    else:
        b = (1 + α) ** ((p - 1) // 2)
        x = b * x0

    return x


def quadratic_roots(b, c):
    """
    Computes roots to the quadratic polynomial

        f = x^2 + b * x + c

    Using the quadratic formula

    Just like in school!
    """
    d2 = b**2 - 4 * c
    d = sqrt_Fp2(d2)
    return ((-b + d) * Fp2_inv_2, -(b + d) * Fp2_inv_2)


def generic_modular_polynomial_roots(j1):
    """
    Compute the roots to the Modular polynomial
    Φ2, setting x to be the input j-invariant.

    When only one j-invariant is known, we
    find up to three new j-invariant values.

    This is fairly slow, but is only done
    once per graph.
    """
    R = PolynomialRing(j1.parent(), "y")
    y = R.gens()[0]
    Φ2 = (
        j1**3
        - j1**2 * y**2
        + 1488 * j1**2 * y
        - 162000 * j1**2
        + 1488 * j1 * y**2
        + 40773375 * j1 * y
        + 8748000000 * j1
        + y**3
        - 162000 * y**2
        + 8748000000 * y
        - 157464000000000
    )

    return Φ2.roots(multiplicities=False)


def quadratic_modular_polynomial_roots(jc, jp):
    """
    When we have the current node's value as
    well as the parent node value then we can
    find the remaining roots by solving a
    quadratic polynomial following
    
    https://ia.cr/2021/1488
    """
    jc_sqr = jc**2
    α = -jc_sqr + 1488 * jc + jp - 162000
    β = (
        jp**2
        - jc_sqr * jp
        + 1488 * (jc_sqr + jc * jp)
        + 40773375 * jc
        - 162000 * jp
        + 8748000000
    )
    # Find roots to x^2 + αx + β
    return quadratic_roots(α, β)


def find_j_invs(j1, l, j_prev=None):
    """
    Compute the j-invariants of the elliptic
    curves 2-isogenous to the elliptic curve
    with j(E) = j1

    The optional param j_prev is used to
    pass through the parent node in the
    graph.

    This is important to stop backtracking,
    but we can also use it to reduce the degree
    of the modular polynomial to make root
    derivation more efficient.
    """
    # TODO: only l=2 supported
    if l != 2:
        raise ValueError("Currently, `find_j_invs` is only implemented for l=2")

    if j_prev:
        roots = quadratic_modular_polynomial_roots(j1, j_prev)

    else:
        roots = generic_modular_polynomial_roots(j1)

    # Dont include the the previous node to avoid backtracking
    return [j for j in roots if j != j_prev]


# ============================================= #
# Construct and find paths in an isogeny graph  #
# ============================================= #

def j_invariant_isogeny_graph(j1, l, e, middle_j_vals=None):
    """
    A depth-first search of the isogeny graph.

    For G1 where we compute the whole graph, dfs or bfs is
    the same. But when we are looking for a collision in
    the leaves, we can supply the end values from G1 as
    middle_j_vals, and stop as soon as we find a collision.

    This means for a isogeny of degree 2^k, the best case for
    G2 would be computing only (k/2) nodes rather than the whole
    graph!

    Actual graph is a list of dictionaries, with levels as elements
    of the list and nodes in the dictionary with node : parent_node
    as data pairs.
    """
    isogeny_graph = [{} for _ in range(e + 1)]

    # Set the root of the graph
    isogeny_graph[0][j1] = None

    # Set stack for DFS
    stack = [(j1, 0)]

    while stack:
        # Get a new element from the stack
        node, level = stack.pop()

        # Parent of node (for backtracking)
        node_parent = isogeny_graph[level][node]

        # Where we will store the next nodes
        child_level = level + 1

        # Set a bool to see if we should check the middle
        check_middle = child_level == e and middle_j_vals is not None

        # Compute all child nodes, except the node which
        # goes back down the tree
        node_children = find_j_invs(node, l, j_prev=node_parent)

        # Compute a dictionary of child : parent nodes for this
        # level in the tree.
        for child in node_children:
            # Add children node to graph
            isogeny_graph[child_level][child] = node

            # Return early when DFS finds the middle j_invariant
            if check_middle and child in middle_j_vals:
                return isogeny_graph, child

            # New nodes to search through
            if child_level != e:
                stack.append((child, child_level))

    return isogeny_graph, None


def j_invariant_path(isogeny_graph, j1, j2, e, reversed_path=False):
    """
    Compute a path through a graph with root j1 and
    last child j2. This is efficient because of our
    data structure for the graph (simply e look ups
    in a dictionary).
    """
    # Make sure the end node is where we expect
    assert j1 in isogeny_graph[0]
    assert j2 in isogeny_graph[e]

    j_path = [j2]
    j = j2
    for k in reversed(range(1, e + 1)):
        j = isogeny_graph[k][j]
        j_path.append(j)

    if not reversed_path:
        j_path.reverse()
    return j_path


def isogeny_from_j_invariant_path(E1, j_invariant_path, l):
    """
    Given a starting curve E1 and a path of j-invariants
    of elliptic curves l-isogenous to its neighbour, compute
    an isogeny ϕ with domain E1 and codomain En with
    j(En) equal to the last element of the path
    """
    # Check we're starting correctly
    assert E1.j_invariant() == j_invariant_path[0]

    # We will compute isogenies linking
    # Ei, Ej step by step
    ϕ_factors = []
    Ei = E1

    for j_step in j_invariant_path[1:]:
        # Compute the isogeny between nodes
        ϕij = brute_force_isogeny_jinv(Ei, j_step, l, 1)

        # Store the factor
        ϕ_factors.append(ϕij)

        # Update the curve Ei
        Ei = ϕij.codomain()

    # Composite isogeny from factors
    ϕ = EllipticCurveHom_composite.from_factors(ϕ_factors)
    return ϕ

# ======================================== #
# Brute force ell isogenies between nodes  #
# ======================================== #

def brute_force_isogeny_jinv(E1, j2, l, e):
    """
    Finds an isogeny of degree l^e, between
    two curves E1, and an unknown curve
    with j-invariant j2

    TODO: write this to be combined with
    `BruteForceSearch` so we don't have so
    much code duplication?
    """

    # Compute the j-invariants of the end
    # points
    jE1 = E1.j_invariant()

    # Degree of isogeny we search for
    D = l**e

    # Handle case when the curves are
    # isomorphic
    if D == 1:
        if jE1 == j2:
            F = E1.base()
            E2 = EllipticCurve(F, j=j2)
            return E1.isomorphism_to(E2)
        else:
            raise ValueError(
                f"A degree 1 isogeny cannot be found, as the curves are not isomorphic"
            )

    # Enumerate through kernel generators
    if e == 1:
        kernels = generate_kernels_division_polynomial(E1, l)
    else:
        kernels = generate_kernels_prime_power(E1, l, e)

    for K in kernels:
        ϕ = EllipticCurveIsogeny(E1, K, degree=D, check=False)
        Eϕ = ϕ.codomain()
        jEϕ = Eϕ.j_invariant()
        if jEϕ == j2:
            return ϕ

    raise ValueError(
        f"No degree {D} isogeny found linking the curves E1 and curve with invariant j2"
    )


def brute_force_isogeny(E1, E2, l, e):
    """
    Finds an isogeny of degree l^e, between
    two curves E1, and E2

    TODO:
    Currently only implemented for prime
    power degrees
    """
    # Degree of isogeny we search for
    D = l**e

    # Handle case when the curves are
    # isomorphic
    if D == 1:
        if E1.is_isomorphic(E2):
            return E1.isomorphism_to(E2), E1(0)
        else:
            raise ValueError(
                f"A degree 1 isogeny cannot be found, as the curves are not isomorphic"
            )

    # Enumerate through kernel generators
    if e == 1:
        kernels = generate_kernels_division_polynomial(E1, l)
    else:
        kernels = generate_kernels_prime_power(E1, l, e)

    for K in kernels:
        ϕ = EllipticCurveIsogeny(E1, K, degree=D, check=False)
        Eϕ = ϕ.codomain()
        if Eϕ.is_isomorphic(E2):
            iso = Eϕ.isomorphism_to(E2)
            ϕ = iso * ϕ
            return ϕ, K

    raise ValueError(
        f"[-] BruteForceSearch failed. No degree {D} isogeny found linking the curves E1 and E2"
    )


# ============================================ #
# Claw-finding attack to recover mitm isogeny  #
# ============================================ #

def claw_finding_attack(E1, E2, l, e):
    """
    Finds a meet in the middle isogeny of
    degree D, between two curves E1, E2
    by searching along two halves of a
    graph of j invariants
    """

    # We'll search sqrt(l**e) from E1
    # and E2 and find a curve in the
    # middle.

    # Let e2 >= e1, as dfs approach means
    # we'll probably not compute all of
    # the second graph.
    e1 = floor(e / 2)
    e2 = e - e1

    # Coerce the elements from Fp4 to Fp2
    # as Ei are supersingular, this is always
    # possible.
    j1 = Fp2(E1.j_invariant())
    j2 = Fp2(E2.j_invariant())

    # Compute the isogeny graph with root j1
    isogeny_graph_e1, _ = j_invariant_isogeny_graph(j1, l, e1)

    # The top level of the graph are the middle
    # j-invariants we want to find a collision with
    e1_middle_vals = set(isogeny_graph_e1[e1].keys())

    # Compute the isogeny graph with root j1
    # `j_invariant_isogeny_graph` will terminate as soon as
    # it finds a j invariant with j_middle in e1_middle_vals
    isogeny_graph_e2, j_middle = j_invariant_isogeny_graph(
        j2, l, e2, middle_j_vals=e1_middle_vals
    )

    # If we didn't find a value in the middle, then there's
    # no isogeny, or we have a bug!
    if not j_middle:
        raise ValueError(f"No isogeny of degree {factor(l**e)} linking E1 and E2.")

    # Given the middle j-inv, compute paths from Ei → Middle
    j1_path = j_invariant_path(isogeny_graph_e1, j1, j_middle, e1)
    j2_path = j_invariant_path(isogeny_graph_e2, j2, j_middle, e2, reversed_path=True)

    # Make sure both paths end in the middle
    assert j1_path[-1] == j2_path[0]

    # Construct a path from j(E1) to j(E2)
    j_path = j1_path + j2_path[1:]

    # Compute the isogeny from the j-inv path
    ϕ = isogeny_from_j_invariant_path(E1, j_path, l)

    # Fix the end of the isogeny with an isomorphism
    E2ϕ = ϕ.codomain()
    iso = E2ϕ.isomorphism_to(E2)

    return iso * ϕ


def meet_in_the_middle_with_kernel(E1, E2, l, e):
    """
    Wrapper for the Claw finding algorithm.
    Additionally computes the kernel of the
    recovered isogeny, which is needed in
    SQISign for other computations.

    Edge cases:

    When e = 0 the curves are isomorphic,
    so just compute the isomorphism and exit.

    When e = 1, we can't do a mitm, so just
    brute force all l-isogenies from one end.
    """
    # Deal with isomorphisms
    if e == 0:
        ϕ = E1.isomorphism_to(E2)
        return ϕ, E1(0)

    # Mitm doesn't work if there's no middle to meet in
    if e == 1:
        ϕ, K = brute_force_isogeny(E1, E2, l, e)
        return ϕ, K

    ϕ = claw_finding_attack(E1, E2, l, e)
    if ϕ == None:
        raise ValueError(
            f"[-] ClawFindingAttack failed. No isogeny of degree {factor(l**e)} linking E1 and E2."
        )
    K = kernel_from_isogeny_prime_power(ϕ)
    return ϕ, K
