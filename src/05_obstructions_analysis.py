"""
05 — Why the Conrey-Snaith route is necessary (Paper B, Section 6)

Summarizes the 5 layers of obstruction that block perturbative approaches
and explains why R3 = |K|^2 at large L bypasses all of them.

Run from repo root: python src/05_obstructions_analysis.py
"""

print("=== 5 layers of obstruction ===\n")

layers = [
    ("Layer 1: sinc(1) = 0",
     "Any kernel perturbation dK gives dY2(r) = 2*sinc(r)*dK(r).\n"
     "At r=1 (typical GUE gap): sinc(1) = sin(pi)/pi = 0 exactly.\n"
     "Decouples perturbation from most relevant scale.\n"
     "NOTE: artifact of linearization. Exact K^BK = sqrt(sinc^2 + eps*dY2) is smooth."),

    ("Layer 2: Sign obstruction",
     "For antisymmetric dK (preserves A1=0): |K+iepsA|^2 >= sinc^2.\n"
     "=> Y2 increases => triples more unequal => <r> decreases => c < 0.\n"
     "Physical correction needs c > 0. Narrowing of p(s) is non-perturbative."),

    ("Layer 3: Born approximation — partial success",
     "Born (1st order in dY2): c_Born = +0.32 (26% of c_emp, correct sign).\n"
     "Higher orders diverge: dY2(1) = -0.273 != 0 => dK ~ 1/(r-1) singular.\n"
     "Born = 26% is the MAXIMUM of any perturbative expansion."),

    ("Layer 4: O(sqrt(T)) cancellation",
     "Prime sum f(r;T) = sum_p (logp/sqrt(p))*sin(2pi*tau_p*r) ~ O(sqrt(T)/logT).\n"
     "Physical c/log^2(T) emerges from perfect cancellation of O(sqrt(T)) terms.\n"
     "No numerical truncation at P_max < T captures this."),

    ("Layer 5: d/L^3 contamination",
     "CS expansion has d ~ -253. Richardson at (80,1000) retains d/80 ~ -3.2.\n"
     "Comparable to b ~ 3.9, suppresses result by 50-70%.\n"
     "Only at L >= 1000: d/1000 ~ -0.25 (< 1% contamination).")
]

for title, desc in layers:
    print(f"--- {title} ---")
    print(desc)
    print()

print("=== Resolution ===")
print("Conrey-Snaith route bypasses ALL 5 layers simultaneously:")
print("  R3 = |K|^2 is phase-immune          (Layers 1-3)")
print("  R3_CS evaluated at exact L           (Layer 4)")
print("  Richardson at (1000, 3000)            (Layer 5)")
