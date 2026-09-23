#!/usr/bin/env python3
"""
Strengthened tail-asymptotic audit for No.20.

This file does NOT change the paper or theorem.  It supplements the original
tail-constant audit by checking actual PMF/asymptotic ratios from exact
finite-n coefficients.

The theorem is asymptotic; practical parameter choices can converge slowly
because neighboring poles are close and large alternating coefficients cancel.
Accordingly, the audit:
  1. checks a small exact case where convergence is visible,
  2. checks monotone approach toward 1 over selected n,
  3. reports (but does not require near-1 convergence for) the practical
     m=64,k=4,b=60,s=2 case.
Standard library only.
"""
from fractions import Fraction
from math import comb, factorial
import importlib.util
from pathlib import Path

HERE = Path(__file__).resolve().parent
spec = importlib.util.spec_from_file_location(
    "core", HERE / "NO20_RCOTBF_EXACT_V0_3.py"
)
core = importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)

def tail_parameters(m,k,b,s):
    if b < 1:
        raise ValueError("tail theorem assumes b>=1")
    if s < 1:
        raise ValueError("s>=1")
    mu = core.moments_L(m,k,b)[0]
    r = Fraction(b**k, m**k)
    c = comb(m,b)
    C = (
        Fraction(c**s,1)
        * (1-r)**(s-1)
        / (mu * r**(s-1) * factorial(s-1))
    )
    return mu, r, c, C

def exact_residence_pmf_to(m,k,b,s,N):
    """Exact PMF coefficients P(R_s=n), 0<=n<=N, using Fraction."""
    mu = core.moments_L(m,k,b)[0]
    q = [core.survival_geom(m,k,b,n) for n in range(N+1)]
    A = [x/mu for x in q]
    L = [Fraction(0)] + [q[n-1]-q[n] for n in range(1,N+1)]
    R = A
    for _ in range(s-1):
        out = [Fraction(0)] * (N+1)
        for i, ri in enumerate(R):
            if not ri:
                continue
            for j in range(1, N-i+1):
                if L[j]:
                    out[i+j] += ri * L[j]
        R = out
    return R

def ratio_at(R,n,C,r,s):
    approx = float(C) * (n ** (s-1)) * (float(r) ** n)
    return float(R[n]) / approx

def main():
    # Small case: convergence is visible and numerically stable.
    small = (3,1,1,2)
    m,k,b,s = small
    _,r,_,C = tail_parameters(*small)
    R = exact_residence_pmf_to(*small,100)
    ns = [10,20,50,100]
    ratios = [ratio_at(R,n,C,r,s) for n in ns]
    print("SMALL_CASE", *small)
    for n,x in zip(ns,ratios):
        print(f"ratio_n{n}={x:.15f}")
    assert all(0 < x < 1.01 for x in ratios)
    assert all(ratios[i] < ratios[i+1] for i in range(len(ratios)-1))
    assert ratios[-1] > 0.98

    # Practical case: deliberately only checks directional convergence.
    # Neighboring poles are close, so convergence is much slower.
    practical = (64,4,60,2)
    m,k,b,s = practical
    _,r,_,C = tail_parameters(*practical)
    R = exact_residence_pmf_to(*practical,200)
    ns2 = [100,150,200]
    ratios2 = [ratio_at(R,n,C,r,s) for n in ns2]
    print("PRACTICAL_CASE", *practical)
    for n,x in zip(ns2,ratios2):
        print(f"ratio_n{n}={x:.15f}")
    assert all(x > 0 for x in ratios2)
    assert ratios2[0] < ratios2[1] < ratios2[2]

    print("TAIL_RATIO_AUDIT_PASS")

if __name__ == "__main__":
    main()
