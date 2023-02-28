# SQISign-SageMath

A SageMath implementation of SQISign following the paper 
[SQISign: compact post-quantum signatures from quaternions and isogenies](https://eprint.iacr.org/2020/1240), 
by Luca De Feo, David Kohel, Antonin Leroux, Christophe Petit, and Benjamin Wesolowski (2020).

## Learning to SQI

Accompanying our code, we have written a blog, [Learning to SQI](https://learningtosqi.github.io), which includes detailed 
write-ups of many of the (sub-)algorithms of SQISign. It has not been written to be self-contained, but rather as
supplementary material to the original [SQISign paper](https://eprint.iacr.org/2020/1240).

We do not claim novelty in any of the posts, but rather hope that this resource will be valuable to whoever is thinking about implementing SQISign, or some other 
protocol which relies on the Deuring correspondence.

## Example Usage

Example of SQISign as a one-round interactive identification protocol between two parties:

```python
sage: from SQISign import SQISign
sage: prover, verifier = SQISign(), SQISign()
sage: prover.keygen()
sage: EA = prover.export_public_key()
sage: E1 = prover.commitment()
sage: phi_ker = verifier.challenge(E1)
sage: S = prover.response(phi_ker)
sage: assert verifier.verify_response(EA, E1, S, phi_ker)
```

Example of signing a message `msg` with SQISign:

```python
sage: from SQISign import SQISign
sage: signer, verifier = SQISign(), SQISign()
sage: msg = b"Learning to SQI!"
sage: signer.keygen()
sage: EA = signer.export_public_key()
sage: sig = signer.sign(msg)
sage: assert verifier.verify(EA, sig, msg)
```

## Project Overview

SQISign itself is implemented in the file [`SQISign.py`](SQISign.py). Reading this code is enough to see
a high-level description of how the main functions: `keygen()`, `commitment()`, `challenge()`, `response()`
and `verify_response()` have been implemented.

If you wish to run SQISign yourself, we include two examples following the above snippets with additional comments:

- [`example_SQISign.sage`](example_SQISign.sage)
- [`example_signing.sage`](example_signing.sage)

### Performance

For the SQISign prime $p_{6983}$, it is expected that our implementation of SQISign will successfully run in about 15 minutes. 
If you want to see everything work more quickly, we also include $p_{\text{toy}}$, which allows SQISign to run in approximately 
30 seconds. By default, SQISign will run with the intended 256-bit prime.

If you wish to run SQISign with the smaller toy parameters, then simply change the last lines of [`parameters.py`](parameters.py)
such that `params = p_toy`.

**A Note on Slowness**:
Our implementation of SQISign is particularly slow as we have not optimised the isogeny computations.
For a full run of the SQISign protocol, we find that about 75% of the total running time is spent computing odd-degree isogenies. 
To use the in-built SageMath functions to compute isogenies, (either with with Vélu's formula or the optimised $\sqrt{elu}$), 
we must work with projective coordinates. SQISign is intended to be implemented with $x$-only arithmetic, allowing simultaneous 
access to the torsion sub-groups of both an elliptic curve and it's quadratic twist when performing isogeny computations.
Rather than implement $x$-only isogenies, we instead work with supersingular elliptic curves over the extension field
$E / \mathbb{F}_{p^4}$. This allows us to work with curves with maximal available torsion $(p+1)(p-1)$, and avoid the need
for computing over the twist.
This has the benefit of allowing us to directly work with the SageMath class, `EllipticCurveIsogeny`, but results
in relatively slow isogenies due to the extension field.

We discuss this, and plans for future work on the page [Future Work](https://learningtosqi.github.io/posts/future-work/) on the
Learning to SQI blog.

### Helper functions

- The file [`ideals.py`](ideals.py) contains helper functions for working with the quaternion algebra $\mathcal{B}_{p, \infty}$, 
  and ideals and orders of the quaternion algebra.
- Similarly, [`isogenies.py`](isogenies.py) contains helper functions for computing and working with isogenies between
  supersingular elliptic curves $E / \mathbb{F}_{p^4}$.
- One step in SQISign requires computing a brute force search for a degree $2^\Delta$ isogeny. 
  We implemented this using a meet-in-the-middle algorithm in [`mitm.py`](mitm.py).
- We implement both the naïve KLPT and SigningKLPT algorithms in [`KLPT.py`](KLPT.py). Although this generally follows
  the seminal paper [On the quaternion l-isogeny path problem](https://arxiv.org/abs/1406.0981), by David Kohel, Kristin Lauter,
  Christophe Petit, Jean-Pierre Tignol, we include many adjustments and edge cases such that our algorithms work as is needed 
  for SQISign.
- At certain points in the KLPT algorithm, we need to enumerate short vectors of a lattice to find good solutions to
  various problems. We use [`fpylll`](https://github.com/fplll/fpylll) 
  to achieve this and the code implementing it is contained in [`lattices.py`](lattices.py).
- The functions responsible for translating between ideals and isogenies following the Deuring correspondence is contained
  in [`deuring.py`](deuring.py). This includes the standard `ideal_to_isogeny()`, as well as the specialised algorithms
  introduced in the SQISign paper: `IdealToIsogenyFromKLPT()`.
- Suitable SQISign parameters are stored as dictionaries in [`parameters.py`](parameters.py). These are then imported
  into [`setup.py`](setup.py), which then computes all the global parameters which are needed for various sub-algorithms 
  of SQISign.
- SQISign computes a large $\ell^e$ degree isogeny during `response()`. Before sending this to the verifier, the isogeny is
  compressed and then similarly decompressed by the verifier. The file [`compression.py`](compression.py) file handles these
  functions.
- Anything else we used a lot but didn't seem to belong anywhere else is stored in [`utilities.py`](utilities.py).


## Future Work

- Implement the new algorithms from [New algorithms for the Deuring correspondence: toward practical and secure SQISign signatures](https://eprint.iacr.org/2022/234), 
  by Luca De Feo, Antonin Leroux, Patrick Longa and Benjamin Wesolowski.
- Once we have the new SQISign algorithms, we can start benchmarking various SQISign parameter sets, such as the ones recently suggested in 
  [Cryptographic Smooth Neighbors](https://eprint.iacr.org/2022/1439), by 
  Giacomo Bruno, Maria Corte-Real Santos, Craig Costello, Jonathan Komada Eriksen, Michael Naehrig, Michael Meyer and Bruno Sterner.
- Currently all isogenies are computed between curves $E / \mathbb{F}\_{p^4}$. We can make our implementation more
  efficient by using $x$-only arithmetic, such that we can access the $(p+1)$ torsion of $E / \mathbb{F}\_{p^2}$ and 
  the $(p-1)$ torsion of its quadratic twist without needing to perform expensive computations on the extension 
  field $\mathbb{F}\_{p^4}$.
- As computing and evaluating isogenies consume most of the computation time, we could write specialised isogeny methods using 
  Montgomery and Edwards curves. These are not suitable for all curves, but should lead to significant performance improvements
  for our SQISign implementation
- Finish implementing the more efficient variation of `ideal_to_kernel()` from 
  [Deuring for the People: Supersingular Elliptic Curves with Prescribed Endomorphism Ring in General Characteristic](https://ia.cr/2023/106), 
  by Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni.

### References

- [SQISign: compact post-quantum signatures from quaternions and isogenies](https://eprint.iacr.org/2020/1240), Luca De Feo, David Kohel, Antonin Leroux, Christophe Petit, and Benjamin Wesolowski (2020).
- [New algorithms for the Deuring correspondence: toward practical and secure SQISign signatures](https://eprint.iacr.org/2022/234), Luca De Feo, Antonin Leroux, Patrick Longa and Benjamin Wesolowski (2022).
- [Quaternion algebras and isogeny-based cryptography](https://www.lix.polytechnique.fr/Labo/Antonin.LEROUX/manuscrit_these.pdf), Antonin Leroux (PhD Thesis) (2022).
- [On the quaternion $\ell$-isogeny path problem](https://arxiv.org/abs/1406.0981), David Kohel, Kristin Lauter, Christophe Petit, Jean-Pierre Tignol (2014).
- [Supersingular isogeny graphs and endomorphism rings: reductions and solutions](https://eprint.iacr.org/2018/371), Kirsten Eisenträger, Sean Hallgren, Kristin Lauter, Travis Morrison Christophe Petit (2018).
- [An improvement to the quaternion analogue of the l-isogeny path problem](https://crypto.iacr.org/2018/affevents/mathcrypt/page.html), Christophe Petit and Spike Smith (2018).
- [Deuring for the People: Supersingular Elliptic Curves with Prescribed Endomorphism Ring in General Characteristic](https://ia.cr/2023/106) Jonathan Komada Eriksen, Lorenz Panny, Jana Sotáková, and Mattia Veroni (2023).
