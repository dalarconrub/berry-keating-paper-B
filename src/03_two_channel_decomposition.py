"""
03 — Two-channel decomposition (Paper B, Section 4)

Computes c_R3 using the parameter-free formula with E_cond (Schur complement).
No threshold needed: E_cond naturally suppresses rare gaps.

Run from repo root: python src/03_two_channel_decomposition.py
"""
import json
import numpy as np
from scipy.linalg import det as la_det

def sinc_kernel(x, y):
    d = x - y
    if abs(d) < 1e-14:
        return 1.0
    return np.sin(np.pi * d) / (np.pi * d)

def p2_components(s1, s2, n_quad=24):
    """det3 and E_cond for consecutive gaps (s1, s2)."""
    S = s1 + s2
    pts = np.array([0.0, s1, S])
    M3 = np.array([[sinc_kernel(pts[i], pts[j]) for j in range(3)] for i in range(3)])
    det3 = np.linalg.det(M3)
    M3_inv = np.linalg.inv(M3)
    nodes, weights = np.polynomial.legendre.leggauss(n_quad)
    x = 0.5 * S * nodes + 0.5 * S
    w = 0.5 * S * weights
    N = len(x)
    K_mat = np.zeros((N, N))
    for i in range(N):
        for j in range(N):
            kxy = sinc_kernel(x[i], x[j])
            kx_pts = np.array([sinc_kernel(x[i], pts[k]) for k in range(3)])
            kpts_y = np.array([sinc_kernel(pts[k], x[j]) for k in range(3)])
            K_mat[i, j] = np.sqrt(w[i]) * (kxy - kx_pts @ M3_inv @ kpts_y) * np.sqrt(w[j])
    return det3, la_det(np.eye(N) - K_mat)

# Load fine grid with precomputed E_cond
with open("data/e4_fine_grid_with_Econd.json") as f:
    data = json.load(f)
pts = data["points"]

R_GUE = 0.59971
dv = 0.05
c_emp = 1.245

# Compute c_R3 for different R3 thresholds
print("=== c_R3 sensitivity to grid coverage ===")
print(f"{'R3_min':>8s} {'N':>5s} {'Z_p2':>8s} {'c_R3':>8s} {'c_E':>8s}")
print("-" * 42)
for R3_min in [0.0, 0.10, 0.20, 0.30, 0.50, 0.80]:
    num, den, n = 0.0, 0.0, 0
    for p in pts:
        if p["R3_GUE"] < R3_min:
            continue
        s1, s2 = p["v1"], p["v2"] - p["v1"]
        r_val = min(s1, s2) / max(s1, s2) if max(s1, s2) > 0 else 0
        b = p.get("b_rich_1k_3k", p.get("b_rich", 0))
        Ec = p["E_cond"]
        p2 = p["R3_GUE"] * Ec
        num += (r_val - R_GUE) * b * Ec * dv**2
        den += p2 * dv**2
        n += 1
    c_R3 = num / den if den > 0 else 0
    print(f"{R3_min:8.2f} {n:5d} {den:8.4f} {c_R3:+8.4f} {c_emp - c_R3:+8.4f}")

print(f"\nc_total = {c_emp} (invariant by construction)")
print("The formula with E_cond has NO free parameters (no threshold needed)")
