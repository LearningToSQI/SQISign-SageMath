"""
Parameter sets for SQISign stored as dictionaries to be imported
into setup.py.

p_toy: A toy 54-bit prime which is suitable for running SQISign and testing 
         functions.

p6983:   The prime from the SQISign paper. 

Parameter set is decided by setting `params` on line 76.
"""

# Sage imports
from sage.all import ZZ

# A 54-bit toy value:
# (p^2 - 1) = 2^17 * 3^2 * 17^2 * 29^2 * 37 * 41 * 43 * 59 * 71 * 79 * 97 * 101 * 103 * 107 * 137
p_toy = {
    "p": ZZ(9568331647090687),
    "q": 1,  # For constructing B(-q, -p)
    "l": 2,
    # All torsion is available for p_toy
    "available_torsion": ZZ(91552970508717179193151202131968),
    # We use all odd torsion p_toy
    # T = available_torsion // 2^17
    "T": ZZ(698493732518899377389154069),
    "e": 213,  # Fixed signing exp. ~ log(p^4)
    "Δ": 5,  # Meet in the middle exp
    "T_prime": 43 * 59 * 71 * 79 * 97 * 101 * 103 * 107 * 137,  # Commitment torsion
    "Dc": 3**2 * 17**2 * 29**2 * 37 * 41,  # Challenge torsion
    "f_step_max": 11,  # Maximum step in IdealToIsogenyFromKLPT
}

p6983 = {
    "p": ZZ(
        73743043621499797449074820543863456997944695372324032511999999999999999999999
    ),
    "q": 1,  # For constructing B(-q, -p)
    "l": 2,
    # 2^34 * 3^53 * 5^21 * 7^2 * 11 * 31 * 43 * 83 * 103^2 * 107 * 109 * 137 * 199 * 227 * 
    # 419 * 491 * 569 * 631 * 677 * 751 * 827 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 
    # 3691 * 4019 * 4283 * 6983
    "available_torsion": ZZ(
        395060348595898919675269756303091126739265710380751097355414818834918473504507691894531931983730526183038976000000000000000000000
    ),
    # 3^53 * 5^21 * 7^2 * 11 * 31 * 43 * 83 * 103^2 * 107 * 109 * 137 * 199 * 227 * 419 * 491 * 
    # 569 * 631 * 677 * 751 * 827 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 3691 * 4019 * 
    # 4283 * 6983
    "T": ZZ(
        22995538811426314040743149601801480558890948909145727366873460055498783098563301138028731335870137332417011260986328125
    ),
    "e": 1000,  # Fixed signing exp. ~ log(p^4)
    "Δ": 15,  # Meet in the middle exp
    # Commitment torsion
    # 7^2 * 11 * 31 * 43 * 83 * 103^2 * 107 * 109 * 137 * 199 * 227 * 419 * 491 * 569 * 631 * 677 * 
    # 751 * 827 * 857 * 859 * 883 * 1019 * 1171 * 1879 * 2713 * 3691 * 4019 * 4283 * 6983
    "T_prime": ZZ(
        2487980652789837077620845135662275571903620147320463972934655988236926594121011
    ),
    "Dc": 3**53 * 5**21,  # Challenge torsion, gcd(Dc, T_prime) = 1
    # Maximum step in IdealToIsogenyFromKLPT
    "f_step_max": 31,
}

# ************************* #
# Pick parameter dictionary #
# ************************* #

# Picking params from the above constants will
# select the parameters used in SQISign

# Expect SQISign to take 30 seconds for p_toy
# and about 15 minutes for p6983

# params = p_toy
params = p6983
