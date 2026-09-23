#!/usr/bin/env python3
"""
NO20 V0.2
Exact stationary residence-time laws for occupancy-triggered rotating Bloom-filter rings.

Model:
- m bits per sub-filter
- k iid uniform bit placements per inserted distinct item, with replacement
- active sub-filter freezes/rotates at the first insertion whose post-insertion
  occupancy exceeds b, 0 <= b < m
- s queryable sub-filters in a rotating ring
- a tagged insertion is sampled at a stationary insertion epoch

This script implements:
(1) exact cycle survival via Stirling occupancy probabilities,
(2) finite-geometric representation,
(3) moments through order 3,
(4) discrete equilibrium residual law,
(5) exact finite-truncation residence PMF and quantiles,
(6) exact mean/variance of residence count,
(7) Poisson-arrival Laplace transform,
(8) independent occupancy-DP and Markov-recursion checks.

No third-party packages required.
"""

from __future__ import annotations
from fractions import Fraction
from functools import lru_cache
from math import comb
import argparse
import random


def check(m: int, k: int, b: int, s: int | None = None) -> None:
    if m < 1 or k < 1 or not (0 <= b < m):
        raise ValueError("require m>=1, k>=1, 0<=b<m")
    if s is not None and s < 1:
        raise ValueError("require s>=1")


@lru_cache(maxsize=None)
def S2(n: int, j: int) -> int:
    if n == 0:
        return 1 if j == 0 else 0
    if j == 0 or j > n:
        return 0
    if j == 1 or j == n:
        return 1
    return j * S2(n - 1, j) + S2(n - 1, j - 1)


def falling(m: int, j: int) -> int:
    p = 1
    for r in range(j):
        p *= m-r
    return p


def coeffs(m: int, b: int):
    return [(-1)**(b-i) * comb(m, i) * comb(m-i-1, b-i)
            for i in range(b+1)]


def survival_stirling(m: int, k: int, b: int, n: int) -> Fraction:
    check(m,k,b)
    if n < 0:
        raise ValueError
    t = k*n
    return Fraction(sum(falling(m,j)*S2(t,j) for j in range(b+1)), m**t)


def survival_geom(m: int, k: int, b: int, n: int) -> Fraction:
    check(m,k,b)
    if n < 0:
        raise ValueError
    ans = Fraction(0)
    for i,c in enumerate(coeffs(m,b)):
        r = Fraction(i**k, m**k)
        ans += c * r**n
    return ans


def pmf_L(m: int, k: int, b: int, n: int) -> Fraction:
    if n < 1:
        return Fraction(0)
    return survival_geom(m,k,b,n-1)-survival_geom(m,k,b,n)


def moments_L(m: int, k: int, b: int):
    mu=e2=e3=Fraction(0)
    for i,c in enumerate(coeffs(m,b)):
        r = Fraction(i**k, m**k)
        one = 1-r
        mu += Fraction(c,1)/one
        e2 += Fraction(c,1)*(1+r)/(one**2)
        e3 += Fraction(c,1)*(1+4*r+r*r)/(one**3)
    var=e2-mu*mu
    return mu,e2,e3,var


def H(m:int,k:int,b:int,z:Fraction)->Fraction:
    ans=Fraction(0)
    for i,c in enumerate(coeffs(m,b)):
        r=Fraction(i**k,m**k)
        ans += Fraction(c,1)/(1-r*z)
    return ans


def G_L(m:int,k:int,b:int,z:Fraction)->Fraction:
    return 1-(1-z)*H(m,k,b,z)


def residual_pmf(m:int,k:int,b:int,a:int)->Fraction:
    if a<0: return Fraction(0)
    mu,_,_,_=moments_L(m,k,b)
    return survival_geom(m,k,b,a)/mu


def residual_moments(m:int,k:int,b:int):
    mu,e2,e3,_ = moments_L(m,k,b)
    ea=(e2-mu)/(2*mu)
    ea2=(2*e3-3*e2+mu)/(6*mu)
    va=ea2-ea*ea
    return ea,ea2,va


def residence_moments(m:int,k:int,b:int,s:int):
    check(m,k,b,s)
    mu,e2,e3,varL=moments_L(m,k,b)
    ea,ea2,varA=residual_moments(m,k,b)
    er=ea+(s-1)*mu
    vr=varA+(s-1)*varL
    return er,vr


def G_R(m:int,k:int,b:int,s:int,z:Fraction)->Fraction:
    check(m,k,b,s)
    mu,_,_,_=moments_L(m,k,b)
    return H(m,k,b,z)/mu * G_L(m,k,b,z)**(s-1)


def poisson_laplace(m:int,k:int,b:int,s:int,lam:Fraction,u:Fraction)->Fraction:
    if lam <= 0 or u < 0:
        raise ValueError
    return G_R(m,k,b,s,lam/(lam+u))


def placement_dp_survival(m:int,k:int,b:int,n:int)->Fraction:
    p=[Fraction(0)]*(m+1)
    p[0]=Fraction(1)
    for _ in range(k*n):
        q=[Fraction(0)]*(m+1)
        for j,x in enumerate(p):
            if x==0: continue
            q[j] += x*Fraction(j,m)
            if j<m:
                q[j+1] += x*Fraction(m-j,m)
        p=q
    return sum(p[:b+1],Fraction(0))


def item_transition(m:int,k:int,j:int):
    p={j:Fraction(1)}
    for _ in range(k):
        q={}
        for x,prob in p.items():
            q[x]=q.get(x,Fraction(0))+prob*Fraction(x,m)
            if x<m:
                q[x+1]=q.get(x+1,Fraction(0))+prob*Fraction(m-x,m)
        p=q
    return p


def markov_mean_cycle(m:int,k:int,b:int)->Fraction:
    n=b+1
    A=[[Fraction(0) for _ in range(n+1)] for __ in range(n)]
    for j in range(n):
        A[j][j]=1
        for h,p in item_transition(m,k,j).items():
            if h<=b:
                A[j][h]-=p
        A[j][n]=1

    for col in range(n):
        piv=next(r for r in range(col,n) if A[r][col])
        A[col],A[piv]=A[piv],A[col]
        d=A[col][col]
        A[col]=[x/d for x in A[col]]
        for r in range(n):
            if r==col: continue
            f=A[r][col]
            if f:
                A[r]=[A[r][c]-f*A[col][c] for c in range(n+1)]
    return A[0][n]


def truncate_distribution_L(m:int,k:int,b:int,eps=1e-14):
    out=[]
    n=1
    tail=1.0
    while tail>eps:
        p=float(pmf_L(m,k,b,n))
        out.append(p)
        tail=float(survival_geom(m,k,b,n))
        n+=1
        if n>1_000_000:
            raise RuntimeError("tail truncation too long")
    return out


def truncate_distribution_A(m:int,k:int,b:int,eps=1e-14):
    mu=float(moments_L(m,k,b)[0])
    out=[]
    a=0
    mass=0.0
    while 1-mass>eps:
        p=float(survival_geom(m,k,b,a))/mu
        out.append(p)
        mass+=p
        a+=1
        if a>1_000_000:
            raise RuntimeError("tail truncation too long")
    return out


def conv(a,b):
    out=[0.0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        if x==0: continue
        for j,y in enumerate(b):
            if y: out[i+j]+=x*y
    return out


def residence_distribution(m:int,k:int,b:int,s:int,eps=1e-14):
    A=truncate_distribution_A(m,k,b,eps)
    L=[0.0]+truncate_distribution_L(m,k,b,eps)
    R=A
    for _ in range(s-1):
        R=conv(R,L)
    return R


def quantile(pmf,q):
    acc=0.0
    for i,p in enumerate(pmf):
        acc+=p
        if acc>=q:
            return i
    return len(pmf)-1


def monte_carlo(m:int,k:int,b:int,s:int,reps:int,seed:int=12345):
    rng=random.Random(seed)

    def one_cycle():
        bits=set()
        n=0
        while True:
            n+=1
            for _ in range(k):
                bits.add(rng.randrange(m))
            if len(bits)>b:
                return n

    pool=[one_cycle() for _ in range(max(5000,reps//2))]
    total=sum(pool)
    cumsum=[]
    z=0
    for L in pool:
        z+=L
        cumsum.append(z)

    import bisect
    vals=[]
    for _ in range(reps):
        slot=rng.randrange(total)
        idx=bisect.bisect_right(cumsum,slot)
        L0=pool[idx]
        start=0 if idx==0 else cumsum[idx-1]
        pos=slot-start
        A=L0-1-pos
        R=A
        for __ in range(s-1):
            R+=one_cycle()
        vals.append(R)
    mean=sum(vals)/len(vals)
    var=sum((x-mean)**2 for x in vals)/len(vals)
    return mean,var


def self_test():
    cases=[(3,1,1,2),(5,2,2,2),(8,3,4,3),(16,3,8,2)]
    for m,k,b,s in cases:
        for n in range(8):
            a=survival_stirling(m,k,b,n)
            g=survival_geom(m,k,b,n)
            d=placement_dp_survival(m,k,b,n)
            assert a==g==d,(m,k,b,n,a,g,d)

        mu,_,_,_=moments_L(m,k,b)
        assert mu==markov_mean_cycle(m,k,b),(m,k,b,mu,markov_mean_cycle(m,k,b))

        R=residence_distribution(m,k,b,s,1e-13)
        mass=sum(R)
        assert abs(mass-1)<1e-10,(m,k,b,s,mass)
        er=float(residence_moments(m,k,b,s)[0])
        er_num=sum(i*p for i,p in enumerate(R))
        assert abs(er-er_num)<1e-8,(m,k,b,s,er,er_num)

        assert G_L(m,k,b,Fraction(1))==1
        assert G_R(m,k,b,s,Fraction(1))==1

    print("SELF_TEST_PASS")


def report(m,k,b,s):
    self_test()
    mu,e2,e3,varL=moments_L(m,k,b)
    ea,_,varA=residual_moments(m,k,b)
    er,vr=residence_moments(m,k,b,s)
    R=residence_distribution(m,k,b,s)
    print(f"PARAM m={m} k={k} b={b} s={s}")
    print(f"E[L]={float(mu):.12g}")
    print(f"Var(L)={float(varL):.12g}")
    print(f"E[A]={float(ea):.12g}")
    print(f"Var(A)={float(varA):.12g}")
    print(f"E[R_s]={float(er):.12g}")
    print(f"Var(R_s)={float(vr):.12g}")
    for q in (0.05,0.50,0.95,0.99):
        print(f"Q{int(q*100):02d}(R_s)={quantile(R,q)}")


def main():
    ap=argparse.ArgumentParser()
    ap.add_argument("--m",type=int,default=64)
    ap.add_argument("--k",type=int,default=4)
    ap.add_argument("--b",type=int,default=32)
    ap.add_argument("--s",type=int,default=2)
    ap.add_argument("--self-test",action="store_true")
    ap.add_argument("--mc",type=int,default=0)
    a=ap.parse_args()
    if a.self_test:
        self_test()
        return
    report(a.m,a.k,a.b,a.s)
    if a.mc:
        th=residence_moments(a.m,a.k,a.b,a.s)
        mm,mv=monte_carlo(a.m,a.k,a.b,a.s,a.mc)
        print(f"MC mean={mm:.8g}, theory={float(th[0]):.8g}")
        print(f"MC var ={mv:.8g}, theory={float(th[1]):.8g}")


if __name__=="__main__":
    main()
