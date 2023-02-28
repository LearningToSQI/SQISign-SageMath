"""
Pari a.log(b) doesn't allow supplying the order of the base `b`. 
We write a new wrapper around `pari.fflog` to speed up dlogs in Fp4

If Pari outperforms in other areas, we can add more functions here :)
"""

# Import pari interface
import cypari2

# SageMath imports
from sage.all import ZZ

# Make instance of Pari
pari = cypari2.Pari()

def discrete_log_pari(a, base, order):
    """
    Wrapper around pari discrete log.
    Works like a.log(b), but allows
    us to use the optional argument
    order.
    """
    x = pari.fflog(a, base, order)
    return ZZ(x)
