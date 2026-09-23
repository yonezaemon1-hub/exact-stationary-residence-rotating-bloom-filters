#!/usr/bin/env python3
"""Tail constant audit for No.20 V0.3. Standard library only."""
from fractions import Fraction
from math import comb, factorial
import importlib.util
from pathlib import Path

HERE=Path(__file__).resolve().parent
spec=importlib.util.spec_from_file_location("core", HERE/"NO20_RCOTBF_EXACT_V0_3.py")
core=importlib.util.module_from_spec(spec)
spec.loader.exec_module(core)

def tail_parameters(m,k,b,s):
    if b < 1:
        raise ValueError("tail theorem assumes b>=1")
    if s < 1:
        raise ValueError("s>=1")
    mu=core.moments_L(m,k,b)[0]
    r=Fraction(b**k,m**k)
    c=comb(m,b)
    K=Fraction(c**s,1)*(1-r)**(s-1)/(mu*r**(s-1))
    C=K/Fraction(factorial(s-1),1)
    # p_n ~ C n^(s-1) r^n
    S=C*r/(1-r)
    # P(R_s>n) ~ S n^(s-1) r^n
    return mu,r,c,C,S

def self_test():
    for m,k,b,s in [(8,1,4,1),(8,1,4,2),(16,2,8,3),(64,4,32,2)]:
        mu,r,c,C,S=tail_parameters(m,k,b,s)
        assert mu>0 and 0<r<1 and c>0 and C>0 and S>0
    print("TAIL_FORMULA_SELF_TEST_PASS")

if __name__=="__main__":
    self_test()
    for p in [(64,4,32,2),(64,4,60,2),(64,4,60,4)]:
        mu,r,c,C,S=tail_parameters(*p)
        print("PARAM",*p)
        print("mu=",float(mu))
        print("r=",float(r))
        print("c=",c)
        print("pmf_constant=",float(C))
        print("survival_constant=",float(S))
