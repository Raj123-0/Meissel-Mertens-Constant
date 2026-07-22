===============================================================================
PROJECT: Meissel-Mertens Constant Computation Engine
===============================================================================

OVERVIEW:
Calculates the Meissel-Mertens constant (M_1 ≈ 0.26149721284764278375...) to 
arbitrary precision (N digits). It is fundamental to number theory and prime 
reciprocal sums (Mertens' second theorem).

ALGORITHM & MATHEMATICS:
- Möbius Inversion Series:
    M_1 = gamma + sum_{k=2}^{infinity} (mu(k) / k) * ln(zeta(k))
- Evaluated with Euler-Mascheroni constant and logarithmic Riemann Zeta values.
