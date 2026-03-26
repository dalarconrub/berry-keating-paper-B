#!/usr/bin/env python3
"""
E4 Anclas multi-L: R₃_CS en 5 puntos × 8 valores de L.
Ajuste 3p robusto: δR₃ = a/L + b/L² + d/L³

Tiempo estimado: ~25 min (5 pts × 8 L × ~38s/eval)
Salida: scripts/e4_anclas_multi_L_results.json

Ejecutar:
  python scripts/e4_anclas_multi_L.py
"""
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(__file__), '..', 'src'))

import mpmath, numpy as np, json, time
from scipy.optimize import curve_fit

mpmath.mp.dps = 20

def sieve(N):
    is_p = np.ones(N+1, dtype=bool); is_p[0]=is_p[1]=False
    for i in range(2, int(N**0.5)+1):
        if is_p[i]: is_p[i*i::i] = False
    return np.where(is_p)[0]

PRIMES = sieve(50000).astype(float)
NP = 3000

def A_f(x):
    r = mpmath.mpf(1)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p); p1x = mpmath.power(p, 1+x)
        r *= (1 - 1/p1x) * (1 - 2/p + 1/p1x) / (1 - 1/p)**2
    return r

def B_f(x):
    r = mpmath.mpf(0)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p)
        r += (mpmath.log(p) / (mpmath.power(p, 1+x) - 1))**2
    return r

def Q_f(x, y):
    r = mpmath.mpf(0)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p)
        r -= mpmath.log(p)**3 / (mpmath.power(p, 2+x+y) *
             (1 - 1/mpmath.power(p, 1+x)) * (1 - 1/mpmath.power(p, 1+y)))
    return r

def P_f(x, y):
    Ax = A_f(x); S = mpmath.mpf(0)
    for p in PRIMES[:NP]:
        p = mpmath.mpf(p); lp = mpmath.log(p)
        px = mpmath.power(p, x); py = mpmath.power(p, y)
        p1x = mpmath.power(p, 1+x); p1y = mpmath.power(p, 1+y)
        num = (1 - 1/px) * (1 - 1/px - 1/py + 1/p1y) * (-lp)
        den = (1 - 1/mpmath.power(p, 1-x+y)) * (1 - 1/p1y) * \
              (1 - 2/p + 1/p1x) * mpmath.power(p, 2-x+y)
        if abs(den) > 1e-30: S += num / den
    return Ax * S

def zpz(s): return mpmath.zeta(s, derivative=1) / mpmath.zeta(s)
def zpz_prime(s):
    z = mpmath.zeta(s); zp = mpmath.zeta(s, derivative=1)
    return mpmath.zeta(s, derivative=2)/z - (zp/z)**2

def I_closed(a1, a2, beta, T):
    L = mpmath.log(T/(2*mpmath.pi)); s1 = a1+beta; s2 = a2+beta
    r = mpmath.mpc(0)
    for sa, sb, aa, ab in [(s1,s2,a1,a2), (s2,s1,a2,a1)]:
        if abs(sa) < 1e-20: continue
        ph = T*mpmath.exp(-sa*L)/(1-sa)
        zz = mpmath.zeta(1-sa)*mpmath.zeta(1+sa)
        Qv = Q_f(sa, sb); Av = A_f(sa); Pv = P_f(sa, sb)
        if abs(ab-aa) > 1e-20 and abs(ab+beta) > 1e-20:
            zd = zpz(1+ab-aa) - zpz(1+ab+beta)
        elif abs(ab-aa) > 1e-20:
            zd = zpz(1+ab-aa) + 1/(ab+beta)
        else: zd = mpmath.mpf(0)
        r += ph*zz*(Qv + Av*zd + Pv)
    return r

def I1_corrected(a, beta, T):
    L = mpmath.log(T/(2*mpmath.pi)); s = a+beta
    if abs(s) < 1e-20: s = mpmath.mpf('1e-10')
    zpzp = zpz_prime(1+s); il = T*(L-1)
    esL = mpmath.exp(-s*L)
    ilp = T*esL*(L/(1-s) + 1/(1-s)**2)
    zz = mpmath.zeta(1+s)*mpmath.zeta(1-s)
    return il*zpzp + ilp*(zz*A_f(s) - B_f(s))

def R3_GUE(v1, v2):
    S = np.sinc; s01 = S(v1); s02 = S(v2); s12 = S(v1-v2)
    return 1 - s01**2 - s02**2 - s12**2 + 2*s01*s02*s12

def R3_CS(v1, v2, Lv):
    Tv = mpmath.exp(Lv)*2*mpmath.pi; rho = Lv/(2*np.pi)
    a1 = mpmath.mpc(0, 2*np.pi*v1/Lv); a2 = mpmath.mpc(0, 2*np.pi*v2/Lv); z = mpmath.mpf(0)
    Is = mpmath.re(I_closed(a1,a2,z,Tv)+I_closed(z,a1,-a2,Tv)+I_closed(z,a2,-a1,Tv)+
                   I_closed(-a1,-a2,z,Tv)+I_closed(z,-a2,a1,Tv)+I_closed(z,-a1,a2,Tv))
    I1s = mpmath.re(I1_corrected(z,a2,Tv)+I1_corrected(z,a1,Tv)+I1_corrected(-a2,a1,Tv)+
                    I1_corrected(-a2,z,Tv)+I1_corrected(-a1,a2,Tv)+I1_corrected(-a1,z,Tv))
    log3 = float(Tv)*(Lv**3 - 3*Lv**2 + 6*Lv - 6)
    return (log3+float(Is)+float(I1s))/((2*np.pi)**3*float(Tv)*rho**3)

# ============================================================
# Configuracion
# ============================================================
POINTS = [(0.7,1.8), (1.0,2.0), (1.0,1.5), (0.8,1.5), (1.0,2.5)]
L_VALS = [15, 20, 25, 30, 40, 50, 60, 80]
OUTPUT = os.path.join(os.path.dirname(__file__), 'e4_anclas_multi_L_results.json')

def m2(L, a, b): return a/L + b/L**2
def m3(L, a, b, d): return a/L + b/L**2 + d/L**3

# ============================================================
# Computo
# ============================================================
print(f'Anclas multi-L: {len(POINTS)} puntos x {len(L_VALS)} L values')
print(f'L = {L_VALS}')
print(f'Estimado: ~{len(POINTS)*len(L_VALS)*38/60:.0f} min')
print('='*60)
sys.stdout.flush()

all_results = []
t_total = time.time()

for idx, (v1, v2) in enumerate(POINTS):
    R3g = R3_GUE(v1, v2)
    dR3_vals = []
    t0 = time.time()

    for Li, Lv in enumerate(L_VALS):
        r3 = R3_CS(v1, v2, float(Lv))
        dR3_vals.append(float(r3) - R3g)
        elapsed = time.time() - t0
        print(f'  ({v1},{v2}) L={Lv}: dR3={dR3_vals[-1]:.6f}  dR3*L2={dR3_vals[-1]*Lv**2:+.2f}  [{elapsed:.0f}s]')
        sys.stdout.flush()

    dt = time.time() - t0
    Ls = np.array(L_VALS, dtype=float)
    dR3 = np.array(dR3_vals)

    p2, cov2 = curve_fit(m2, Ls, dR3)
    p3, cov3 = curve_fit(m3, Ls, dR3, p0=[p2[0], p2[1], 0])
    errs3 = np.sqrt(np.diag(cov3))

    rec = {
        'v1': v1, 'v2': v2, 'R3_GUE': float(R3g),
        'L_vals': L_VALS,
        'dR3': [float(x) for x in dR3_vals],
        'fit_2p': {'a': float(p2[0]), 'b': float(p2[1])},
        'fit_3p': {'a': float(p3[0]), 'b': float(p3[1]), 'd': float(p3[2]),
                   'a_err': float(errs3[0]), 'b_err': float(errs3[1]), 'd_err': float(errs3[2])},
    }
    all_results.append(rec)

    print(f'  -> 2p: a={p2[0]:+.4f} b={p2[1]:+.2f}')
    print(f'  -> 3p: a={p3[0]:+.4f}+-{errs3[0]:.4f}  b={p3[1]:+.2f}+-{errs3[1]:.2f}  d={p3[2]:+.1f}+-{errs3[2]:.1f}')
    print(f'  [{idx+1}/{len(POINTS)}] {dt:.0f}s')
    print()
    sys.stdout.flush()

    # Guardar parcial
    with open(OUTPUT, 'w') as f:
        json.dump(all_results, f, indent=2)

dt_total = time.time() - t_total
print('='*60)
print(f'Total: {dt_total/60:.1f} min')
print(f'Guardado: {OUTPUT}')
print()
print('RESUMEN b (3p):')
for r in all_results:
    f3 = r['fit_3p']
    print(f'  ({r["v1"]},{r["v2"]}): b={f3["b"]:+.2f}+-{f3["b_err"]:.2f}  a={f3["a"]:+.4f}  d={f3["d"]:+.1f}')
b_all = [r['fit_3p']['b'] for r in all_results]
print(f'  mean(b) = {np.mean(b_all):+.2f}  std = {np.std(b_all):.2f}')
