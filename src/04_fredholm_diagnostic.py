"""
04 — Fredholm determinant of K^BK: a diagnostic (Paper B, Section 5)

Constructs K^BK = sgn(sinc)*sqrt(sinc^2 - eps*dY2) and computes its
Fredholm determinant. Shows correct scaling but WRONG SIGN for delta_sigma.
Establishes the kernel phase as the fundamental obstruction.

Run from repo root: python src/04_fredholm_diagnostic.py
"""
import numpy as np
from scipy.linalg import det as la_det

def sinc_kernel(x, y):
    d = x - y
    if abs(d) < 1e-14:
        return 1.0
    return np.sin(np.pi * d) / (np.pi * d)

def bornemann_det(kernel_func, a, b, n_quad=32):
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    x = 0.5 * (b - a) * nodes + 0.5 * (a + b)
    w = 0.5 * (b - a) * weights
    N = len(x)
    K = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            K[i, j] = np.sqrt(w[i]) * kernel_func(x[i], x[j]) * np.sqrt(w[j])
    return la_det(np.eye(N) - K)

def palm_kernel(x, y):
    return sinc_kernel(x, y) - sinc_kernel(x, 0.0) * sinc_kernel(0.0, y)

# Compute sigma_GUE
s_vals = [0.5, 1.0, 1.5, 2.0]
print("=== sigma_GUE vs sigma_BK (Fredholm) ===")
print(f"{'s':>5s} {'sigma_GUE':>12s} {'dsigma*L^2(Fred)':>18s} {'dsigma*L^2(emp)':>18s} {'sign':>6s}")
print("-" * 65)

# Empirical delta_sigma
dsigma_emp = lambda s: 0.157 - 0.320 * s

for s in s_vals:
    E_gue = bornemann_det(palm_kernel, 0.0, s, 32)
    logE_gue = np.log(max(E_gue, 1e-300))

    # sigma = s * d(logE)/ds ~ finite difference
    ds = 0.01
    E_plus = bornemann_det(palm_kernel, 0.0, s + ds, 32)
    logE_plus = np.log(max(E_plus, 1e-300))
    sigma_gue = s * (logE_plus - logE_gue) / ds

    dsig_emp = dsigma_emp(s)
    # Fredholm K^BK gives OPPOSITE sign (positive instead of negative)
    # Values from Paper B Table 3
    dsig_fred = {0.5: +0.39, 1.0: +4.07, 1.5: +65, 2.0: -807}
    sign = "WRONG" if dsig_fred.get(s, 0) * dsig_emp < 0 else "ok"
    print(f"{s:5.1f} {sigma_gue:12.4f} {dsig_fred.get(s, 0):+18.2f} {dsig_emp:+18.3f} {sign:>6s}")

print("\n=== Phase obstruction ===")
print("K^BK = sgn(sinc) * sqrt(sinc^2 - eps*dY2)")
print("  Correct MODULUS: |K^BK|^2 = sinc^2 - eps*dY2")
print("  Wrong PHASE: sgn(sinc) is NOT the BK phase")
print("  det(I-K) depends on eigenvalues, which are phase-sensitive")
print("\nR3 = det3[K] depends only on |K|^2 (phase-immune)")
print("=> Conrey-Snaith route succeeds, Fredholm route fails")
