# Python imports
import cProfile
import pstats
import time

# Local imports
from mitm import *
from isogenies import *
from setup import *

# Speedy and still (mostly) correct
proof.all(False)


l, e = 2, 14
D = l^e
P, Q = torsion_basis(E0, D)

total_time = 0
iter_count = 100
for _ in range(iter_count):
    K = P + randint(0, D) * Q 
    K.set_order(D)
    ϕA = E0.isogeny(K, algorithm="factored")
    EA = ϕA.codomain()

    t0 = time.time()
    ϕ = claw_finding_attack(E0, EA, l, e)
    total_time += time.time() - t0
    assert ϕ.codomain() == EA

print(f"Mitm time taken: {(total_time/iter_count):.5f}")


# K = P + randint(0, D) * Q
# ϕA = E0.isogeny(K, algorithm="factored")
# EA = ϕA.codomain()

# cProfile.run("ClawFindingAttack(E0, EA, l, e)", 'mitm.cProfile')
# p = pstats.Stats('mitm.cProfile')
# p.strip_dirs().sort_stats("cumtime").print_stats(int(25))