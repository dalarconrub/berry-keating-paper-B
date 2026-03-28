"""
02 — Richardson extrapolation at large L (Paper B, Section 3)

Loads precomputed R3_CS at L=1000,3000 (fine grid, 1082 pts),
performs Richardson extrapolation, and shows convergence of b.
Generates Figure 1 (heatmap of b) and Figure 4 (b by pair).

Run from repo root: python src/02_richardson_extrapolation.py
"""
import json
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Load fine grid data
with open("data/e4_fine_grid_with_Econd.json") as f:
    data = json.load(f)
pts = data["points"]

print(f"Loaded {len(pts)} fine grid points (dv=0.05)")

# Extract arrays
v1s = np.array([p["v1"] for p in pts])
v2s = np.array([p["v2"] for p in pts])
bs = np.array([p.get("b_rich_1k_3k", p.get("b_rich", 0)) for p in pts])
r3s = np.array([p["R3_GUE"] for p in pts])

# Mirror for symmetry (data only has v1 < v2)
v1_all, v2_all, b_all = [], [], []
for v1, v2, b in zip(v1s, v2s, bs):
    v1_all.extend([v1, v2])
    v2_all.extend([v2, v1])
    b_all.extend([b, b])
v1_all = np.array(v1_all)
v2_all = np.array(v2_all)
b_all = np.array(b_all)

# Statistics by R3 region
print(f"\n{'R3 range':>12s} {'N':>5s} {'<b>':>8s} {'std(b)':>8s}")
print("-" * 38)
for lo, hi in [(0.0, 0.30), (0.30, 0.50), (0.50, 0.80), (0.80, 1.01)]:
    mask = (r3s >= lo) & (r3s < hi)
    if mask.sum() == 0:
        continue
    print(f"[{lo:.2f},{hi:.2f}) {mask.sum():5d} {bs[mask].mean():+8.2f} {bs[mask].std():8.2f}")

# --- Figure 1: heatmap ---
v1_unique = np.sort(np.unique(v1_all))
v2_unique = np.sort(np.unique(v2_all))
dv = np.median(np.diff(v1_unique))
b_grid = np.full((len(v2_unique), len(v1_unique)), np.nan)
for i in range(len(v1_all)):
    ix = np.argmin(np.abs(v1_unique - v1_all[i]))
    iy = np.argmin(np.abs(v2_unique - v2_all[i]))
    b_grid[iy, ix] = b_all[i]

x_edges = np.append(v1_unique - dv/2, v1_unique[-1] + dv/2)
y_edges = np.append(v2_unique - dv/2, v2_unique[-1] + dv/2)

fig, ax = plt.subplots(figsize=(7, 6))
im = ax.pcolormesh(x_edges, y_edges, b_grid, cmap="RdBu_r", vmin=-20, vmax=20, shading="flat")
plt.colorbar(im, ax=ax, shrink=0.85, label="$b_{rich}$ [Rich(1000, 3000)]")
ax.set_xlabel("$v_1$", fontsize=13)
ax.set_ylabel("$v_2$", fontsize=13)
ax.set_title("Richardson coefficient $b(v_1, v_2)$", fontsize=13)
ax.set_aspect("equal")
lim = [min(v1_unique[0], v2_unique[0]), max(v1_unique[-1], v2_unique[-1])]
ax.plot(lim, lim, "k--", alpha=0.3, lw=0.8)
plt.tight_layout()
plt.savefig("figures/fig1_b_rich_map.pdf", dpi=300, bbox_inches="tight")
plt.savefig("figures/fig1_b_rich_map.png", dpi=150, bbox_inches="tight")
print("\nSaved figures/fig1_b_rich_map.{pdf,png}")
