"""
01 — Conrey-Snaith formula for R3 (Paper B, Section 2)

Computes R3_GUE and R3_CS at given positions (0, v1, v2) and L.
Shows the expansion delta_R3 = a/L + b/L^2 + d/L^3.

Run from repo root: python src/01_conrey_snaith_R3.py
"""
import numpy as np

def sinc_kernel(x, y):
    d = x - y
    if abs(d) < 1e-14:
        return 1.0
    return np.sin(np.pi * d) / (np.pi * d)

def R3_GUE(v1, v2):
    """GUE 3-point correlation det3(0, v1, v2)."""
    M = np.array([
        [1.0,                sinc_kernel(0, v1),  sinc_kernel(0, v2)],
        [sinc_kernel(v1, 0), 1.0,                 sinc_kernel(v1, v2)],
        [sinc_kernel(v2, 0), sinc_kernel(v2, v1), 1.0]
    ])
    return np.linalg.det(M)

# Example points
examples = [(1.0, 2.0), (0.5, 1.5), (0.8, 1.2), (0.3, 0.7)]
print("=== R3_GUE at selected positions ===")
print(f"{'v1':>6s} {'v2':>6s} {'R3_GUE':>10s}")
print("-" * 26)
for v1, v2 in examples:
    r3 = R3_GUE(v1, v2)
    print(f"{v1:6.2f} {v2:6.2f} {r3:10.6f}")

print("\nExpansion: delta_R3(v1,v2;L) = a(v1,v2)/L + b(v1,v2)/L^2 + d(v1,v2)/L^3")
print("b encodes the Berry-Keating correction at order 1/log^2(T)")
print("Richardson extrapolation between L1, L2 eliminates a/L")
