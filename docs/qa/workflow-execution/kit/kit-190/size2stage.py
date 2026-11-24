# KIT-190 bound-declaration.md section 8: simulated size of the two-stage procedure at the margin (stdlib Python only;
# indicative; no package code). Run: python3 -I size2stage.py. The figures in the declaration came from the functions
# below with 20,000 (effect) and 6,000 (null) simulated cells per setting, seeds as printed there.
# Size of the two-stage procedure at the margin (stdlib only; indicative).
import random, math
from statistics import NormalDist, mean, stdev
N=NormalDist(); z=N.inv_cdf(0.95); T=0.80**(1/8)
def r_effect(s, rmin, rmax=1000, m=0.05):
    for R in range(rmin, rmax+1, 50):
        if N.cdf(m*math.sqrt(R)/s - z) >= T: return R
    return rmax
def r_null(d, rmin, rmax=1000, m=0.10):
    d=max(d, 1e-6)
    for R in range(rmin, rmax+1, 50):
        if N.cdf(m*math.sqrt(R)/math.sqrt(d) - z) >= T: return R
    return rmax
def qhat(b,c,n,d):
    lo,hi=1e-12,(1-d)/2-1e-12
    f=lambda q: b/(q+d)+(c/q if c else 0)-2*(n-b-c)/(1-2*q-d)
    for _ in range(100):
        mid=(lo+hi)/2
        if f(mid)>0: lo=mid
        else: hi=mid
    return (lo+hi)/2
def tango(b,c,n,m=0.10):
    q=qhat(b,c,n,m); return (b-c-n*m)/math.sqrt(n*(2*q+m*(1-m))) < -z
def size_null(delta, sims, seed, rmin=150):
    rnd=random.Random(seed); p10=(delta+0.10)/2; p01=(delta-0.10)/2; hits=0; Rs=[]
    def draw(n):
        b=c=0
        for _ in range(n):
            u=rnd.random()
            if u<p10: b+=1
            elif u<p10+p01: c+=1
        return b,c
    for _ in range(sims):
        b1,c1=draw(50); R=r_null((b1+c1)/50, rmin); b2,c2=draw(R-50); Rs.append(R)
        hits+=tango(b1+b2,c1+c2,R)
    return hits/sims, mean(Rs)
def size_effect(s, sims, seed, rmin=200, skew=True):
    rnd=random.Random(seed); hits=0; Rs=[]
    def draw(n):
        # skewed: centred exponential scaled to SD s, plus the margin 0.05 (true mean exactly at the margin)
        return [0.05 + s*((rnd.expovariate(1.0)-1.0) if skew else rnd.gauss(0,1)) for _ in range(n)]
    for _ in range(sims):
        x1=draw(50); R=r_effect(stdev(x1), rmin); x=x1+draw(R-50); Rs.append(R)
        e=[v-0.05 for v in x]; t=mean(e)/(stdev(e)/math.sqrt(len(e)))
        hits+= t < -z   # t reference approximates the studentised sign-flip test
    return hits/sims, mean(Rs)
for d in (0.12,0.25,0.40):
    print('null delta',d, size_null(d, 3000, int(d*1000)))
for s in (0.20,0.40):
    for sk in (False,True):
        print('effect sd',s,'skewed' if sk else 'normal', size_effect(s, 4000, int(s*100)+sk))
