
import time

from setup import *
from isogenies import EllipticCurveIsogenyFactored

P, Q = E0.gens()

n = (p^2 - 1) // T_prime
K = n * P
K.set_order(T_prime)

t0 = time.time()
phi = E0.isogeny(K, algorithm="factored")
t1 = time.time() - t0
print(f"Boring: {t1}")

t0 = time.time()
phi = EllipticCurveIsogenyFactored(E0, K, order=T_prime, velu_bound=10000000000)
t1 = time.time() - t0
print(f"Inf bound: {t1}")

for bound in [200, 400, 800]:
    t0 = time.time()
    EllipticCurveIsogenyFactored(E0, K, order=T_prime, velu_bound=bound)
    t1 = time.time() - t0
    print(f"Bound: {bound}, Time: {t1}")