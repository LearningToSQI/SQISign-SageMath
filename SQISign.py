"""
SageMath implementation of SQISign

    SQISign: compact post-quantum signatures from quaternions and isogenies
    Luca De Feo, David Kohel, Antonin Leroux, Christophe Petit, and Benjamin Wesolowski,
    https://ia.cr/2020/1240

# ========================== #
# Proof of knowledge example #
# ========================== #

sage: from SQISign import SQISign
sage: prover, verifier = SQISign(), SQISign()
sage: prover.keygen()
sage: EA = prover.export_public_key()
sage: E1 = prover.commitment()
sage: ϕ_ker = verifier.challenge(E1)
sage: S = prover.response(ϕ_ker)
sage: assert verifier.verify_response(EA, E1, S, ϕ_ker)

# ========================== #
#       Signing example      #
# ========================== #

sage: from SQISign import SQISign
sage: signer, verifier = SQISign(), SQISign()
sage: msg = b"Learning to SQI!"
sage: signer.keygen()
sage: EA = signer.export_public_key()
sage: sig = signer.sign(msg)
sage: assert verifier.verify(EA, sig, msg)

# ========================== #
#      SQISign Functions     #
# ========================== #

keygen(): Generates two equivalent ideals Iτ (prime norm)
    and Jτ (smooth norm). Computes the isogeny 
    τ_prime : E0 → EA = E0 / <Jτ>. 

    (τ_prime, Iτ, Jτ) are secret values,
    EA is the public value.

    All values are stored in `self`

export_public_key(): Returns the public key EA from self.
    Requires that .keygen() has been run.

commitment(): Computes a secret isogeny ψ : E0 → E1 of degree
    T_prime, together with the corresponding ideal Iψ. 
    (ψ, Iψ) are stored in self.
    Returns the public value E1 for use in generating a 
    challenge.

challenge(): Computes a public isogeny ϕ : E1 → E2 of degree
    Dc. Returns ϕ_ker.

challenge_from_message(): Given a message `msg` and curve E1, 
    computes an isogeny ϕ : E1 → E2 deterministically of degree 
    Dc and returns ϕ_ker

response(): Given an isogeny ϕ : E1 → EA and the secret values
    from .keygen() and .commitment() computes an isogeny 
    σ : EA → E2 of fixed degree l^e. 

sign(): Given a message `msg` computes a random commitment ψ and
    then generates a challenge from the commitment and the 
    message using `challenge_from_message()`. A response to
    the generated challenge is computed and returned.

verify_response(): Given the public key EA and the response σ
    check whether σ has the right degree and codomains and 
    whether ϕ_dual ∘ σ is a cyclic isogeny : EA → E1.

verify(): Given a message, and the signature (E1, σ) generates a 
    challenge ϕ and runs `verify_response()` to verify the
    signature.
"""

# Python imports
from hashlib import shake_128

# SageMath imports
from sage.all import randint, ZZ, factor, proof

# Local imports
from ideals import (
    is_integral,
    is_cyclic,
    multiply_ideals,
    equivalent_left_ideals,
    left_isomorphism,
)
from isogenies import torsion_basis, dual_isogeny, EllipticCurveIsogenyFactored
from deuring import IdealToIsogenyFromKLPT, kernel_to_ideal
from KLPT import RepresentIntegerHeuristic, SigningKLPT
from compression import compression, decompression
from utilities import inert_prime, has_order_D
from setup import E0, O0, Bτ, eτ, p, l, Dc, T_prime, ω, e, f_step_max

proof.all(False)


class SQISign:
    def __init__(self):
        """
        TODO

        We only use the object to store intermediate values,
        we could also init this with parameters and stop using
        globals.
        """
        # Key pair
        # pk = EA
        # sk = (τ_prime, Iτ, Jτ)
        self.pk = None
        self.sk = None

        # Secret commitment values
        # commitment_secrets = (ψ_ker, ψ, Iψ)
        self.commitment_secrets = None

    def keygen(self):
        """
        Efficient keygen as described in Appendix D
        of the SQISign paper

        Input: None
        Output: None

        Stores to self:
            self.pk = EA
            self.sk = (τ_prime, Iτ, Jτ)

            EA: the codomain of the isogeny τ_prime
            τ_prime: the secret isogeny from E0 → EA
            Iτ: ideal with prime norm equivalent to Jτ
            Jτ: ideal with smooth norm, equivalent to Iτ and to
                τ_prime under the Deuring correspondence

        Note:
            To send the public key use the function
            `self.export_public_key()`
        """
        # Compute a random prime ≤ Bτ which is inert
        # in R[ω].
        # Note: this is the same as picking p ≡ 3 mod 4
        Nl = l**eτ
        
        # Stop infinite loops
        for _ in range(1000):
            Nτ = inert_prime(Bτ, -ZZ(ω**2))
            # We need the product to be large enough for
            # RepresentIntegerHeuristic.
            if Nτ * Nl > 2 * p:
                break

        # Compute an endomorphism γ of norm Nτ l^eτ
        # Nτ < Bτ
        γ = None

        # Stop infinite loops
        for _ in range(1000):
            γ = RepresentIntegerHeuristic(Nτ * Nl, parity=True)
            if γ is not None:
                break
        
        if γ is None:
            exit("Never found an alg element with norm (Nτ * Nl), Exiting...")

        Iτ = O0 * γ + O0 * Nτ
        Jτ = O0 * γ.conjugate() + O0 * Nl

        # Iτ has prime norm
        assert Iτ.norm() == Nτ, f"{Iτ.norm() = }, {Nτ = }, {Nl = }"
        # Jτ has smooth norm l^e
        assert Jτ.norm() == Nl, f"{Jτ.norm() = }, {Nτ = }, {Nl = }"

        # Iτ is an integral ideal: Iτ ⊂ O0
        assert is_integral(Iτ), "Iτ is not integral"

        # Jτ is an integral ideal: Iτ ⊂ O0
        assert is_integral(Jτ), "Jτ is not integral"

        # Jτ is a cyclic isogeny
        assert is_cyclic(Jτ), "Jτ is not cyclic"

        # Compute the secret isogeny τ
        I_trivial = O0.unit_ideal()
        ϕ_trivial = E0.isogeny(E0(0))
        τ_prime = IdealToIsogenyFromKLPT(
            Jτ, I_trivial, ϕ_trivial, I_prime=Iτ, end_close_to_E0=True
        )
        EA = τ_prime.codomain()

        # The isogeny τ_prime should have degree = n(Jτ)
        assert (
            τ_prime.degree() == Jτ.norm()
        ), f"{factor(τ_prime.degree()) = } {factor(Jτ.norm()) = }"

        self.pk = EA
        self.sk = (τ_prime, Iτ, Jτ)

        return None

    def export_public_key(self):
        """
        Helper function to return the public key

        TODO: this could be compressed, probably.
        """
        if self.pk is None:
            raise ValueError(f"Must first generate a keypair with `self.keygen()`")
        return self.pk

    def commitment(self):
        """
        Compute the challenge isogeny and corresponding ideal
        of degree / norm T_prime

        Input: None
        Output: E1: the codomain of the commitment isogeny

        Stores to self:
        self.commitment_secrets = (ψ_ker, ψ, Iψ))
            ψ_ker: the kernel of ψ
            ψ: the secret commitment isogeny ψ : E0 → E1
            Iψ: the ideal equivalent to ψ.
        """
        # Generate a random kernel
        # of order T_prime
        P, Q = torsion_basis(E0, T_prime)
        x = randint(1, T_prime)
        ψ_ker = P + x * Q

        # Generate the ideal Iψ from ψ_ker
        Iψ = kernel_to_ideal(ψ_ker, T_prime)
        assert Iψ.norm() == T_prime, "Iψ has the wrong norm"

        # Generate the ideal ψ and its codomain
        ψ = EllipticCurveIsogenyFactored(E0, ψ_ker, order=T_prime)
        E1 = ψ.codomain()

        # Store secret results
        self.commitment_secrets = (ψ_ker, ψ, Iψ)

        return E1

    @staticmethod
    def challenge(E1, x=None):
        """
        Compute the challenge isogenys

        Input: E1 the codomain of the commitment and domain
               of the challenge isogeny
        Output: ϕ_ker: The kernel isogeny ϕ : E1 → E2 of degree Dc
        """
        # Generate a random kernel ∈ E1[Dc]
        E1.set_order((p**2 - 1) ** 2)
        P, Q = torsion_basis(E1, Dc, canonical=True)

        # If x isn't supplied, generate a random x
        if x is None:
            x = randint(1, Dc)

        # Compute the kernel of the challenge isogeny
        ϕ_ker = P + x * Q

        return ϕ_ker

    def challenge_from_message(self, E1, msg):
        """
        Compute a challenge deterministically from a
        message

        Input: E1: the codomain of the commitment and domain
               of the challenge isogeny
               msg: the message to be signed

        Output: ϕ_ker: The kernel isogeny ϕ : E1 → E2 of degree Dc

        TODO: this was just thrown together, almost certainly not
              what we should be doing here.
        """
        # Compute a scalar from the message
        h = shake_128(msg).digest(128)
        x = int.from_bytes(h, "big")

        # Reduce modulo Dc
        x = ZZ(x % Dc)

        # Compute a challenge using a kernel
        # K = P + [x]Q ∈ E1[Dc]
        return self.challenge(E1, x=x)

    def response(self, ϕ_ker):
        """
        Compute the isogeny σ : EA → E2 of degree l^e where
        e is a SQISign parameter. Does this by via the Deuring
        correspondence from an ideal of norm l^e.

        Input:  ϕ_ker: The kernel isogeny ϕ : E1 → E2 of degree Dc
        Output: S: a bitstring corresponding to an isogeny σ : EA → E2
        """
        if self.pk is None or self.sk is None:
            raise ValueError(f"Must first generate a keypair with `self.keygen()`")

        if self.commitment_secrets is None:
            raise ValueError(
                f"Must first generate a commitment with `self.commitment()`"
            )

        # Extract secret values from keygen
        EA = self.pk
        τ_prime, Iτ, Jτ = self.sk

        # Extract values from commitment
        ψ_ker, ψ, Iψ = self.commitment_secrets

        # Recover the dual of ψ from ψ and its kernel
        ψ_dual = dual_isogeny(ψ, ψ_ker, order=T_prime)

        # Deviation from paper time!
        # We are asked to first compute Iϕ
        # Then compute: Iτ_bar * Iψ * Iϕ
        # But we don't actually do this.
        # Instead, we directly compute
        # Iψ * Iϕ = Iψ ∩ I_([ψ]^* ϕ)
        #         = Iψ ∩ I_([ψ_dual]_* ϕ)
        #

        # First compute the ideal from the pullback
        # I_([ψ_dual]_* ϕ)
        Iϕ_pullback = kernel_to_ideal(ψ_dual(ϕ_ker), Dc)
        IψIϕ = Iψ.intersection(Iϕ_pullback)
        assert IψIϕ.norm() == Iψ.norm() * Iϕ_pullback.norm()

        # Compute the product of ideals
        # I = Iτ_bar * Iψ * Iϕ
        Iτ_bar = Iτ.conjugate()
        I = multiply_ideals(Iτ_bar, IψIϕ)
        assert I.norm() == Iτ_bar.norm() * IψIϕ.norm()

        print(f"INFO [SQISign Response]: Running SigningKLPT")
        J = SigningKLPT(I, Iτ)
        assert J.norm() == l**e, "SigningKLPT produced an ideal with incorrect norm"
        print(f"INFO [SQISign Response]: Finished SigningKLPT")

        assert equivalent_left_ideals(
            I, J
        ), "Signing KLPT did not produce an equivalent ideal!"
        assert is_cyclic(J), "SigningKLPT produced a non-cyclic ideal"

        # Ensure that the left and right orders match
        α = left_isomorphism(Iτ, Jτ)
        J = α ** (-1) * J * α
        assert J.left_order() == Jτ.right_order()

        print(f"INFO [SQISign Response]: Computing the corresponding isogeny")
        σ = IdealToIsogenyFromKLPT(J, Jτ, τ_prime, K_prime=Iτ)
        print(f"INFO [SQISign Response]: Computed the isogeny EA → E2")

        print(f"INFO [SQISign Response]: Compressing the isogeny σ to a bitstring")
        S = compression(EA, σ, l, f_step_max)
        print(
            f"INFO [SQISign Response]:"
            f"Compressed the isogeny σ to a bitstring of length {len(S)}"
        )

        return S

    def sign(self, msg):
        """
        Use SQISign to sign a message by creating a challenge
        isogeny from the message and generating a response S
        from the challenge.

        Input: msg: the message to be signed

        Output: sig: a signature tuple (E1, S)
                    E1 : the codomain of the commitment
                    S: a compressed bitstring of the response isogeny EA → E2
        """
        # Make a commitment
        E1 = self.commitment()

        # Use the message to find a challenge
        ϕ_ker = self.challenge_from_message(E1, msg)

        # Compute a response for the challenge
        S = self.response(ϕ_ker)

        return (E1, S)

    def verify_response(self, EA, E1, S, ϕ_ker):
        """
        Verify that the compressed bitstring S corresponds to
        an isogeny σ EA → E2 of degree l^e such that ϕ_dual ∘ σ
        is cyclic

        Input: EA: the public key, and codomain of the secret isogeny τ_prime
               E1: the codomain of the secret commitment ψ : E0 → E1
               S: a compressed bitstring of the response isogeny EA → E2
               ϕ_ker: the kernel of the challenge isogeny ϕ : E1 → E2
        Output: True if the response is value, False otherwise
        """
        # Compute the challenge isogeny from the challenge kernel
        ϕ = EllipticCurveIsogenyFactored(E1, ϕ_ker, order=Dc)
        E2 = ϕ.codomain()
        E2.set_order((p**2 - 1) ** 2)

        # Decompress σ
        print(f"INFO [SQISign Verify]: Decompressing the isogeny σ from a bitstring")
        σ = decompression(EA, E2, S, l, f_step_max, e)

        print(f"INFO [SQISign Verify]: Verifying the degree and (co)domains of σ")
        # Ensure that the domain of σ is EA
        if not σ.domain() == EA:
            print(f"DEBUG [SQISign Verify]: The domain of σ is not EA")
            return False

        if not σ.codomain() == E2:
            print(f"DEBUG [SQISign Verify]: The codomain of σ is not E2")
            return False

        # Check the degree of σ is as expected
        if ZZ(σ.degree()) != l**e:
            print(
                f"DEBUG [SQISign Verify]:"
                f"The degree σ is {factor(σ.degree())}, expected {l}^{e}"
            )
            return False

        # Check that the isogeny ϕ_dual ∘ σ is cyclic
        print(f"INFO [SQISign Verify]: Verifying that ϕ_dual * σ is cyclic")

        # Compute torsion basis EA[2^f]
        D = l**f_step_max
        P, Q = torsion_basis(EA, D)
        ϕ_dual = dual_isogeny(ϕ, ϕ_ker)

        # Compute ϕ_dual ∘ σ : EA → E1
        ϕ_dual_σ = ϕ_dual * σ
        imP = ϕ_dual_σ(P)
        assert imP.curve() == E1, "Mapping is incorrect"

        # Check if ϕ_dual ∘ σ is cyclic
        if has_order_D(imP, D):
            return True

        print(
            f"DEBUG [SQISign Verify]: ϕ_dual_σ(P) does not have full order, checking Q"
        )

        imQ = ϕ_dual_σ(Q)
        assert imQ.curve() == E1, "Mapping is incorrect"
        if has_order_D(imQ, D):
            return True

        print(f"DEBUG [SQISign Verify]: ϕ_dual_σ is not cyclic!")
        return False

    def verify(self, EA, sig, msg):
        """
        Wrapper for verify for when the challenge must be
        generated from a message

        Input: EA: the public key, and codomain of the secret isogeny τ_prime
               sig: a signature tuple (E1, S)
                   E1: the codomain of the secret commitment ψ : E0 → E1
                   S: a compressed bitstring of the response isogeny EA → E2
               msg: the message which has been signed
        Output: True if the response is value, False otherwise
        """
        # Extract pieces from signature
        E1, S = sig

        # Generate ϕ_ker from the message
        ϕ_ker = self.challenge_from_message(E1, msg)

        # Verify signature
        return self.verify_response(EA, E1, S, ϕ_ker)
