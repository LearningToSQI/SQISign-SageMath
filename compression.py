"""
Functions for the compression and decompression of the 
response isogeny for SQISign.

---

Note:

For SQISign, the suggested compressed representation is

        S = S1 || s2 || S2 || ... || sv || Sv

where Si are bit strings representing solutions to the
kernels Ki = Ri + [Si]Qi for Ri,Qi ∈ Ei[D] and si are 
integers which hint to computing the Ri orthogonal to Qi 
in some deterministic manner.

In order for there to always be kernels P + xQ, we need that
the image of σ1(Q1) has full order. As σi have degrees of
prime power, then we are guaranteed that either σ1(P1) or 
σ1(Q2) has full order, but not both.

When σ1(Q1) has order != 2^f we perform a swap P,Q = Q,P
ensuring that the kernel can be written as P + [x]Q. However, 
this has to be communicated to decompression, so we increase the 
size of the compression by 1-bit and append a "1" if we swapped
the first torsion basis and "0" otherwise:

S = swap_bit || S1 || s2 || S2 || ... || sv || Sv

Note, as this fixed σ1(Q1) to have full order, and the kernel
has K = P + [x]Q, then Qi is always orthogonal to Ki and σi(Qi)
will always have full order. So we only need this swap bit for
the first run, when we compute a random (but deterministic) torsion
basis for EA[2^f].
"""

# Sage imports
from sage.all import ZZ
from sage.schemes.elliptic_curves.hom_composite import EllipticCurveHom_composite

# Local imports
from isogenies import (
    torsion_basis,
    generate_random_point,
    EllipticCurveIsogenyFactored,
    DLP
)
from utilities import has_order_D

# ========================================= #
#  Functions to pack an isogeny into blocks #
# ========================================= #

def isogeny_into_blocks(σ, block_size):
    """
    This is a bit of a hack to deal with nested
    composite isogenies.

    The problem is, σ is of type EllipticCurveHom_composite
    which has as elements, other isogenies of type EllipticCurveHom_composite
    which contains prime degree isogenies and could also contain isomorphisms
    of degree 1.

    This function recursively looks for prime isogenies / morphisms
    and then puts them into blocks of at size block_size, then as
    output gives an array of isogenies of degree dividing block_size
    (Only the last element should have degree < block_size).
    """

    def update_blocks(σ_block, σ_blocks):
        """
        Take an array of isogenies of combined degree
        block_size and make a single composite isogeny
        """
        σi = EllipticCurveHom_composite.from_factors(σ_block)
        σ_blocks.append(σi)
        return 1, []

    def rec_prime_degree(σ, block_degree, σ_block, σ_blocks):
        """
        Recursively look for prime degree isogenies and morphisms
        and update blocks when appropriate
        """
        # Is the factor also a composite?
        for σi in σ.factors():
            if isinstance(σi, EllipticCurveHom_composite):
                block_degree, σ_block, σ_blocks = rec_prime_degree(
                    σi, block_degree, σ_block, σ_blocks
                )

            else:
                σ_block.append(σi)
                block_degree *= σi.degree()
                if block_degree == block_size:
                    block_degree, σ_block = update_blocks(σ_block, σ_blocks)

        return block_degree, σ_block, σ_blocks

    σ_blocks = []
    σ_block = []
    block_degree = 1

    block_degree, σ_block, σ_blocks = rec_prime_degree(
        σ, block_degree, σ_block, σ_blocks
    )

    # Update for the last block
    if σ_block:
        update_blocks(σ_block, σ_blocks)

    # Make sure the blocks line up with the input
    assert σ_blocks[0].domain() == σ.domain()
    assert σ_blocks[-1].codomain() == σ.codomain()

    # Make sure the chain σv ∘ ... ∘ σ1 connects
    for v in range(len(σ_block) - 1):
        assert σ_blocks[v].codomain() == σ_blocks[v + 1].domain()

    return σ_blocks

# ========================================= #
#    Convert between bit strings and data   #
# ========================================= #

def bitstring_to_data(σ_compressed, f):
    """
    Given a bit string of data in the following form:

    S = swap_bit || S1 || s2 || S2 || ... || sv || Sv

    Compute the swap_bit as an integer and two arrays 
    of integers [S1 ... Sv] and [s2 ... sv] 
    """
    # Extract out the swap_bit and first dlog
    swap_and_S1, σ_compressed = σ_compressed[: f + 1], σ_compressed[f + 1 :]
    swap_bit, S1 = swap_and_S1[0], swap_and_S1[1:]

    # Arrays to store values as Integers
    swap_bit = ZZ(swap_bit, 2)
    dlogs = [ZZ(S1, 2)]
    hints = []

    # Each chunk si || Si has 4 + f bits
    bit_count = 4 + f

    # TODO
    # I'm parsing the string twice here, which is annoying,
    # This could be made more efficient...
    # Split the bitstring into blocks of size (4 + f)
    siSi_bitstrings = [
        σ_compressed[i : i + bit_count] for i in range(0, len(σ_compressed), bit_count)
    ]
    # Extract out si, Si and parse as integers
    for siSi_bits in siSi_bitstrings:
        si_bits, Si_bits = siSi_bits[:4], siSi_bits[4:]

        # Parse bitstrings to integers
        dlogs.append(ZZ(Si_bits, 2))
        hints.append(ZZ(si_bits, 2))

    return swap_bit, dlogs, hints


def int_to_bitstring(x, n):
    """
    Returns a bit string representing x
    of length exactly n by left padding
    with zeros.

    TODO:
    Assumes that x < 2^n. We could add a check
    for this.
    """
    return bin(x)[2:].zfill(n)


def hint_to_bitstring(hint):
    """
    If the hint is smaller than 15, store
    the binary representation of hint. Else
    set the bit string to "1111", which
    communicates that the first 15
    generated points can be skipped.
    """
    # Only use 4 bits, so the max hint is
    # 15
    hint = min(hint, 15)
    return int_to_bitstring(hint, 4)


# ========================================= #
# Helpers for compression and decompression #
# ========================================= #

def compute_R(E, Q, D, hint=None):
    """
    Deterministically generate a point R linearly
    independent from Q. If `hint` is given, we can
    skip known bad points.

    Input: E an elliptic curves
           a point Q ∈ E[D]
           D the order of Q

    Output: A point R ∈ E[D] such that
            E[D] = <R, Q>
            A hint stating the number of
            iterations taken to find R
    """
    # We need R to have order D, so
    # compute the cofactor n
    p = E.base().characteristic()
    n = (p**2 - 1) // D

    # If the hint is smaller than 2^k
    # then we know the `hint` point is
    # correct
    if hint is not None and hint < 15:
        return n * generate_random_point(E, seed=hint), None

    # If hint is not none, then we know we can
    # skip the first 15 values
    k_start = 0
    if hint is not None:
        k_start = 15

    for k in range(k_start, 2000):
        # Find a random point
        R = n * generate_random_point(E, seed=k)

        # Point is not in E[D]
        if R.is_zero() or not has_order_D(R, D):
            continue

        # We now have a point in E[D]
        # check it's linearly independent to Q
        pair = R.weil_pairing(Q, D, algorithm="pari")
        if has_order_D(pair, D, multiplicative=True):
            R._order = ZZ(D)
            return R, k

    raise ValueError(f"Never found a point Q")


def compute_S_and_next_Q(σ, P, Q, f, first_step=False):
    """
    Given a torsion basis P, Q finds x
    such that the kernel of σ is P + xQ

    Additionally computes S, which is
    the bitstring of fixed length f
    representing x and the image σ(Q)
    """
    # Recover degree of isogeny
    D = σ.degree()

    # Map through points
    imP, imQ = σ(P), σ(Q)

    # For the first step, imQ may not have order D
    # which means there will be no solution to the
    # dlog.
    #
    # We deal with this by swapping P,Q (as D is a prime
    # power, imP has order D when imQ doesn't).
    # We must communicate to decompression that this swap
    # occurred, which we do with the `swap_bit`.
    swap_bit = 0

    if first_step and not has_order_D(imQ, D):
        print(f"DEBUG [compute_S_and_next_Q]: swapped the basis around")
        P, Q = Q, P
        imP, imQ = imQ, imP
        swap_bit = 1

    x = DLP(-imP, imQ, D)
    S = int_to_bitstring(x, f)

    # The isogeny with kernel <K> will only be the same as
    # σ up to isomorphism. To ensure that Q_next is the
    # correct point, we must make σ_new from <K> and then
    # re-evaluate σ_new(Q) to compute the image of a point
    # linearly independent from ker(σ_new).
    K = P + x * Q
    Eσ = K.curve()
    p = Eσ.base().characteristic()
    Eσ.set_order((p**2 - 1)**2, num_checks=0)
    σ_new = EllipticCurveIsogenyFactored(Eσ, K, order=D)

    # Now compute the new point Q_(i+1)
    Q_next = σ_new(Q)
    return S, Q_next, swap_bit

# ============================= #
# Compression and Decompression #
# ============================= #

def compression(E, σ, l, f):
    """
    Given an isogeny σ of degree l^e = l^vf compute a
    compressed representation of the isogeny as a bit string
    in the form:

    σ_compressed = swap_bit || S1 || s2 || S2 || ... || sv || Sv

    Note: swap_bit denotes when the torsion basis of the first
    step must be swapped: R, Q = Q, R to ensure a solution to
    the discrete log. This bit is needed to communicate to decom.
    to do the same swap.
    """
    σ_compressed = ""

    # Split σ into v isogenies of degree f
    σ_chain = isogeny_into_blocks(σ, l**f)
    assert E == σ_chain[0].domain(), "Supplied curve is incorrect"

    # Step 1, we need to compute the ker(σ1) and
    # a point linearly dependent ker(σ1)
    # To do this we need the torsion basis
    # E0[D].
    σ1 = σ_chain[0]
    D = σ1.degree()

    R1, Q1 = torsion_basis(E, D, canonical=True)
    S1, Qi, swap_bit = compute_S_and_next_Q(σ1, R1, Q1, f, first_step=True)

    # Update compressed bitstring
    σ_compressed += str(swap_bit)
    σ_compressed += S1

    # For the remaining steps, we can use Qi as one
    # basis element so we only need to compute some
    # Ri linearly independent to Qi. There will always
    # be a solution to the dlog as Qi has order D.
    for σi in σ_chain[1:]:
        Ei = Qi.curve()

        # We need to align the next step in the chain
        # with σ_new computed previously, which means
        # mapping the domain of the next step with the
        # codomain of the last step.
        assert Ei.is_isomorphic(σi.domain())
        iso = Ei.isomorphism_to(σi.domain())
        σi = σi * iso
        assert Ei == σi.domain()

        D = σi.degree()
        # The last element of the chain has degree D | l^f
        # So we can ensure Qi ∈ E[D] by multiplying by a
        # cofactor
        if D != l**f:
            cofactor = l**f // D
            Qi = cofactor * Qi

        # Add hint to compression to help decompression
        # recover Ri
        Ri, hint = compute_R(Ei, Qi, D)
        σ_compressed += hint_to_bitstring(hint)

        # Add dlog to compression and derive next Qi
        Si, Qi, _ = compute_S_and_next_Q(σi, Ri, Qi, f)
        σ_compressed += Si

    # Our compressed rep. is 1 bit longer as we encode the swap bit
    v = len(σ_chain)
    assert len(σ_compressed) == (v - 1) * (f + 4) + f + 1

    return σ_compressed


def decompression(E_start, E_end, σ_compressed, l, f, σ_length):
    """
    Given a bit string:

    σ_compressed = swap_bit || S1 || s2 || S2 || ... || sv || Sv

    Compute the isogeny σ : E_start → E_end of degree l^σ_length isogeny
    """
    # Extract integers from the encoded bitstring
    swap_bit, dlogs, hints = bitstring_to_data(σ_compressed, f)

    # Compute canonical torsion basis E[D]
    D = l**f
    Ri, Qi = torsion_basis(E_start, D, canonical=True)

    # In compression, if im(Q) does not have full order,
    # we swap R,Q so we can solve the discrete log. If
    # this happened, we need to also swap R,Q in decom.
    if swap_bit == 1:
        Qi, Ri = Ri, Qi

    σ_factors = []
    Ei = E_start
    for Si, hint in zip(dlogs, hints):
        Ki = Ri + Si * Qi
        σi = EllipticCurveIsogenyFactored(Ei, Ki, order=D)
        σ_factors.append(σi)
        Ei = σi.codomain()
        Qi = σi(Qi)
        Ri, _ = compute_R(Ei, Qi, D, hint=hint)

    # The last step has length D | 2^f.
    # I can't see a way to derive it, but as the response length
    # is public, I think we just have to include it as a param...
    # Anyway...
    # When the last step != D, we need to clear the cofactor
    # to make sure Ri, Qi have order = σv.degree()
    last_step_length = σ_length % f
    if last_step_length != 0:
        cofactor = D // l**last_step_length
        Ri = cofactor * Ri
        Qi = cofactor * Qi
        D = l**last_step_length

    Si = dlogs[-1]
    Ki = Ri + Si * Qi
    σi = EllipticCurveIsogenyFactored(Ei, Ki, order=D)
    σ_factors.append(σi)

    σ = EllipticCurveHom_composite.from_factors(σ_factors)
    Eσ = σ.codomain()

    assert Eσ.is_isomorphic(
        E_end
    ), "The isogeny σ does not end at a curve isomorphic to E_end"
    # Use an isomorphism to ensure the codomain
    # of σ is E_end
    iso = Eσ.isomorphism_to(E_end)
    return iso * σ
