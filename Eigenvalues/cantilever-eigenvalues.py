# =============================================================================
#  Eigenvalue Analysis of a Cantilever Steel Pipe
#  Using OpenSeesPy
#
#  Model:   Vertical cantilever, fixed at base, free at top
#  Section: Hollow circular pipe (steel)
#  Method:  Consistent mass matrix, eigensolver
# =============================================================================

import openseespy.opensees as ops
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.patches as mpatches

# =============================================================================
# 1. INPUT PARAMETERS  —  modify these freely
# =============================================================================

# --- Geometry ---
L      = 10.0        # Total length of the cantilever [m]
n_elem = 10          # Number of beam elements along the length

# --- Pipe cross-section [m] ---
D_out  = 0.2         # Outer diameter [m]
D_in   = 0.18        # Inner diameter [m]  (wall thickness = (D_out-D_in)/2 = 10 mm)

# --- Steel material ---
E      = 210e9       # Young's modulus [Pa]
rho    = 7850.0      # Mass density [kg/m³]

# --- Analysis ---
n_modes = 3          # Number of modes to compute

# =============================================================================
# 2. DERIVED CROSS-SECTION PROPERTIES
# =============================================================================

A  = (np.pi / 4.0) * (D_out**2 - D_in**2)   # Cross-sectional area [m²]
Iz = (np.pi / 64.0) * (D_out**4 - D_in**4)  # Second moment of area [m⁴]
m  = rho * A                                 # Mass per unit length [kg/m]

print("=" * 55)
print("  PIPE CROSS-SECTION PROPERTIES")
print("=" * 55)
print(f"  Outer diameter       D_out = {D_out*1000:.1f} mm")
print(f"  Inner diameter       D_in  = {D_in*1000:.1f} mm")
print(f"  Wall thickness       t     = {(D_out-D_in)/2*1000:.1f} mm")
print(f"  Area                 A     = {A*1e4:.4f} cm²")
print(f"  2nd moment of area   Iz    = {Iz*1e8:.4f} cm⁴")
print(f"  Linear mass          m     = {m:.4f} kg/m")

# =============================================================================
# 3. THEORETICAL REFERENCE (Euler-Bernoulli cantilever)
#    f_n = (beta_n * L)^2 / (2*pi*L^2) * sqrt(EI / (rho*A))
#    beta_n * L : 1.8751, 4.6941, 7.8548  for modes 1, 2, 3
# =============================================================================

beta_L = [1.87510, 4.69409, 7.85476]
f_theory = [(bL**2 / (2 * np.pi * L**2)) * np.sqrt(E * Iz / (rho * A))
            for bL in beta_L]

print("\n" + "=" * 55)
print("  THEORETICAL FREQUENCIES (Euler-Bernoulli)")
print("=" * 55)
for i, f in enumerate(f_theory):
    print(f"  Mode {i+1}:  f = {f:.4f} Hz    T = {1/f:.4f} s")

# =============================================================================
# 4. OPENSEESPY MODEL
# =============================================================================

ops.wipe()

# 2D model — 2 spatial dimensions, 3 DOFs per node (ux, uy, rotation)
ops.model('basic', '-ndm', 2, '-ndf', 3)

# --- Nodes ---
# Node i is placed at height y_i = i * (L / n_elem)
n_nodes = n_elem + 1
node_spacing = L / n_elem

for i in range(n_nodes):
    ops.node(i + 1, 0.0, i * node_spacing)

# --- Boundary condition: fix base node (node 1) ---
ops.fix(1, 1, 1, 1)   # ux=fixed, uy=fixed, rz=fixed

# --- Elastic material section (tag = 1) ---
# elasticBeamColumn uses E, A, Iz directly
# We also pass the mass per unit length for consistent mass matrix
ops.geomTransf('Linear', 1)

for i in range(n_elem):
    ops.element(
        'elasticBeamColumn',
        i + 1,          # element tag
        i + 1,          # node i  (start)
        i + 2,          # node i+1 (end)
        A,              # cross-section area
        E,              # Young's modulus
        Iz,             # moment of inertia
        1,              # geometric transformation tag
        '-mass', m,     # mass per unit length → consistent mass matrix
        '-cMass'        # use consistent (not lumped) mass
    )

# =============================================================================
# 5. EIGENVALUE ANALYSIS
# =============================================================================

# Run eigensolver — returns list of eigenvalues λ = ω²
eigenvalues = ops.eigen('-fullGenLapack', n_modes)

omega = [np.sqrt(lam) for lam in eigenvalues]   # angular frequency [rad/s]
freq  = [w / (2 * np.pi) for w in omega]        # frequency [Hz]
period = [1.0 / f for f in freq]                 # period [s]

print("\n" + "=" * 55)
print("  OPENSEESPY EIGENVALUE RESULTS")
print("=" * 55)
print(f"  {'Mode':<6} {'ω [rad/s]':>12} {'f [Hz]':>12} {'T [s]':>12}  {'Δf vs theory':>14}")
print("  " + "-" * 53)
for i in range(n_modes):
    err = abs(freq[i] - f_theory[i]) / f_theory[i] * 100
    print(f"  {i+1:<6} {omega[i]:>12.4f} {freq[i]:>12.4f} {period[i]:>12.4f}  {err:>12.3f} %")

# =============================================================================
# 6. EXTRACT MODE SHAPES  (transverse displacement DOF = DOF 1 at each node)
# =============================================================================

node_heights = [i * node_spacing for i in range(n_nodes)]
mode_shapes  = []

for mode_idx in range(1, n_modes + 1):
    shape = []
    for node_id in range(1, n_nodes + 1):
        ux = ops.nodeEigenvector(node_id, mode_idx, 1)   # DOF 1 = horizontal (ux)
        shape.append(ux)
    shape = np.array(shape)
    # Normalise so the free tip has magnitude 1
    if abs(shape[-1]) > 1e-12:
        shape /= abs(shape[-1])
    mode_shapes.append(shape)

# =============================================================================
# 7. PLOT
# =============================================================================

fig, axes = plt.subplots(1, 2, figsize=(12, 7))
fig.suptitle(
    f"Cantilever Steel Pipe — Eigenvalue Analysis\n"
    f"L = {L} m | D_out = {D_out*1000:.0f} mm | D_in = {D_in*1000:.0f} mm | "
    f"E = {E/1e9:.0f} GPa | {n_elem} elements",
    fontsize=11, fontweight='bold'
)

colors = ['#e63946', '#2a9d8f', '#e9c46a']

# ── Left plot: mode shapes ────────────────────────────────────────────────────
ax1 = axes[0]
ax1.set_title("Mode Shapes (normalised to tip displacement = 1)", fontsize=10)

for i, shape in enumerate(mode_shapes):
    label = (f"Mode {i+1}: f = {freq[i]:.4f} Hz, "
             f"T = {period[i]:.4f} s")
    ax1.plot(shape, node_heights, color=colors[i], linewidth=2.5,
             marker='o', markersize=3.5, label=label)

ax1.axvline(0, color='black', linewidth=0.8, linestyle='--')
ax1.axhline(0, color='black', linewidth=2.0)       # base
ax1.set_xlabel("Normalised Transverse Displacement  [-]", fontsize=9)
ax1.set_ylabel("Height along cantilever  [m]", fontsize=9)
ax1.legend(fontsize=8.5, loc='upper left')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.set_xlim(-1.3, 1.3)

# ── Right plot: frequency comparison ─────────────────────────────────────────
ax2 = axes[1]
ax2.set_title("OpenSeesPy vs. Analytical Frequency Comparison", fontsize=10)

x = np.arange(1, n_modes + 1)
w = 0.35
bars1 = ax2.bar(x - w/2, f_theory[:n_modes], w,
                label='Analytical (Euler-Bernoulli)',
                color='#457b9d', alpha=0.85)
bars2 = ax2.bar(x + w/2, freq, w,
                label='OpenSeesPy',
                color='#e63946', alpha=0.85)

# Annotate bars with values
for bar in bars1:
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
             f"{bar.get_height():.4f}", ha='center', va='bottom', fontsize=7.5)
for bar in bars2:
    ax2.text(bar.get_x() + bar.get_width()/2, bar.get_height() + 0.005,
             f"{bar.get_height():.4f}", ha='center', va='bottom', fontsize=7.5)

ax2.set_xlabel("Mode Number", fontsize=9)
ax2.set_ylabel("Frequency  [Hz]", fontsize=9)
ax2.set_xticks(x)
ax2.set_xticklabels([f"Mode {i}" for i in x])
ax2.legend(fontsize=9)
ax2.grid(True, axis='y', linestyle=':', alpha=0.6)

plt.tight_layout()
# plt.savefig("/mnt/user-data/outputs/cantilever_eigenvalue.png", dpi=150, bbox_inches='tight')
plt.show()
# print("\n  Plot saved → cantilever_eigenvalue.png")
# print("=" * 55)