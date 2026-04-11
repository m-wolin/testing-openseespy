import openseespy.opensees as ops
import matplotlib.pyplot as plt
import matplotlib.gridspec as gridspec
import numpy as np


def MomentCurvature(secTag, axialLoad, maxK, numIncr=100):

    # Define two nodes at (0,0)
    ops.node(1, 0.0, 0.0)
    ops.node(2, 0.0, 0.0)

    # Fix all degrees of freedom except axial and bending
    ops.fix(1, 1, 1, 1)
    ops.fix(2, 0, 1, 0)

    # Define element
    #                             tag ndI ndJ  secTag
    ops.element('zeroLengthSection', 1, 1, 2, secTag)

    # Define constant axial load
    ops.timeSeries('Constant', 1)
    ops.pattern('Plain', 1, 1)
    ops.load(2, axialLoad, 0.0, 0.0)

    # Define analysis parameters
    ops.integrator('LoadControl', 0.0)
    ops.system('SparseGeneral', '-piv')
    ops.test('NormUnbalance', 1e-9, 10)
    ops.numberer('Plain')
    ops.constraints('Plain')
    ops.algorithm('Newton')
    ops.analysis('Static')

    # Do one analysis for constant axial load
    ops.analyze(1)

    # Define reference moment
    ops.timeSeries('Linear', 2)
    ops.pattern('Plain', 2, 2)
    ops.load(2, 0.0, 0.0, 1.0)

    # Compute curvature increment
    dK = maxK / numIncr

    # Use displacement control at node 2 for section analysis
    ops.integrator('DisplacementControl', 2, 3, dK, 1, dK, dK)

    # Do the section analysis — collect curvature and moment at each step
    curvature = [0.0]
    moment    = [0.0]

    for _ in range(numIncr):
        ops.analyze(1)
        curvature.append(ops.nodeDisp(2, 3))
        moment.append(ops.eleForce(1, 3))   # moment at node 1 of element = -M

    return np.array(curvature), np.array(moment)


ops.wipe()
print("Start MomentCurvature.py example")

# Define model builder
# --------------------
ops.model('basic', '-ndm', 2, '-ndf', 3)

# Define materials for nonlinear columns
# ------------------------------------------
# CONCRETE                  tag   f'c        ec0   f'cu        ecu
# Core concrete (confined)
ops.uniaxialMaterial('Concrete01', 1, -6.0, -0.004, -5.0, -0.014)

# Cover concrete (unconfined)
ops.uniaxialMaterial('Concrete01', 2, -5.0, -0.002, 0.0, -0.006)

# STEEL
# Reinforcing steel
fy = 60.0      # Yield stress
E = 30000.0    # Young's modulus

#                        tag  fy E0    b
ops.uniaxialMaterial('Steel01', 3, fy, E, 0.01)

# Define cross-section for nonlinear columns
# ------------------------------------------

# set some paramaters
colWidth = 15
colDepth = 24

cover = 1.5
As = 0.60     # area of no. 7 bars

# some variables derived from the parameters
y1 = colDepth / 2.0
z1 = colWidth / 2.0

ops.section('Fiber', 1)

# Create the concrete core fibers
ops.patch('rect', 1, 10, 1, cover - y1, cover - z1, y1 - cover, z1 - cover)

# Create the concrete cover fibers (top, bottom, left, right)
ops.patch('rect', 2, 10, 1, -y1, z1 - cover, y1, z1)
ops.patch('rect', 2, 10, 1, -y1, -z1, y1, cover - z1)
ops.patch('rect', 2, 2, 1, -y1, cover - z1, cover - y1, z1 - cover)
ops.patch('rect', 2, 2, 1, y1 - cover, cover - z1, y1, z1 - cover)

# Create the reinforcing fibers (left, middle, right)
ops.layer('straight', 3, 3, As, y1 - cover, z1 - cover, y1 - cover, cover - z1)
ops.layer('straight', 3, 2, As, 0.0,        z1 - cover, 0.0,        cover - z1)
ops.layer('straight', 3, 3, As, cover - y1, z1 - cover, cover - y1, cover - z1)

# Estimate yield curvature
# (Assuming no axial load and only top and bottom steel)
# d -- from cover to rebar
d = colDepth - cover
# steel yield strain
epsy = fy / E
Ky = epsy / (0.7 * d)

# Print estimate to standard output
print("Estimated yield curvature: ", Ky)

# Set axial load
P = -180.0

# Target ductility for analysis
mu = 15.0

# Number of analysis increments
numIncr = 100

# Call the section analysis procedure
curvature, moment = MomentCurvature(1, P, Ky * mu, numIncr)

results = open('results.out', 'a+')

u = ops.nodeDisp(2, 3)
if abs(u - 0.00190476190476190541) < 1e-12:
    results.write('PASSED : MomentCurvature.py\n')
    print("Passed!")
else:
    results.write('FAILED : MomentCurvature.py\n')
    print("Failed!")

results.close()

print("==========================")

# ===========================================================================
# POST-PROCESSING & PLOTS
# ===========================================================================

# Negate moment so positive curvature → positive moment
moment = -moment

# Identify yield point (first time moment stops increasing by much)
Ky_actual  = Ky
My_approx  = np.interp(Ky_actual, curvature, moment)

# Ultimate values
Mu = moment.max()
Ku = curvature[moment.argmax()]

# Ductility ratio
ductility = Ku / Ky_actual

fig = plt.figure(figsize=(13, 8))
fig.suptitle(
    f"Moment–Curvature Analysis\n"
    f"RC Section  {colWidth}\" × {colDepth}\"  |  "
    f"cover = {cover}\"  |  As = {As} in²  |  P = {P} kips",
    fontsize=11, fontweight='bold'
)

gs = gridspec.GridSpec(2, 2, figure=fig, hspace=0.42, wspace=0.35)

# ── 1. Moment–Curvature curve ────────────────────────────────────────────────
ax1 = fig.add_subplot(gs[:, 0])   # spans both rows on the left

ax1.plot(curvature, moment, color='#e63946', linewidth=2.2, label='M–φ curve')

# Yield point
ax1.axvline(Ky_actual, color='#2a9d8f', linewidth=1.2,
            linestyle='--', label=f'Est. yield  φ_y = {Ky_actual:.5f}')
ax1.axhline(My_approx, color='#2a9d8f', linewidth=1.2, linestyle='--')
ax1.plot(Ky_actual, My_approx, 'o', color='#2a9d8f', markersize=8, zorder=5)

# Ultimate point
ax1.plot(Ku, Mu, 's', color='#e9c46a', markersize=9, zorder=5,
         label=f'Peak  M_u = {Mu:.1f} kip·in')

# Elastic slope line
K_line = np.linspace(0, Ky_actual * 1.05, 20)
M_line = My_approx / Ky_actual * K_line
ax1.plot(K_line, M_line, color='#457b9d', linewidth=1.4,
         linestyle=':', label='Initial stiffness')

# Annotations
ax1.annotate(f'φ_y ≈ {Ky_actual:.4f}', xy=(Ky_actual, My_approx),
             xytext=(Ky_actual * 1.8, My_approx * 0.72),
             arrowprops=dict(arrowstyle='->', color='#2a9d8f'),
             fontsize=8.5, color='#2a9d8f')
ax1.annotate(f'M_u = {Mu:.1f}', xy=(Ku, Mu),
             xytext=(Ku * 0.55, Mu * 1.04),
             arrowprops=dict(arrowstyle='->', color='#e9c46a'),
             fontsize=8.5, color='#9e6f00')

ax1.set_xlabel('Curvature  φ  [1/in]', fontsize=10)
ax1.set_ylabel('Moment  M  [kip·in]', fontsize=10)
ax1.set_title('Moment–Curvature', fontsize=10, fontweight='bold')
ax1.legend(fontsize=8, loc='lower right')
ax1.grid(True, linestyle=':', alpha=0.6)
ax1.set_xlim(left=0)
ax1.set_ylim(bottom=0)

# ── 2. Normalised M/My vs φ/φy  (ductility view) ────────────────────────────
ax2 = fig.add_subplot(gs[0, 1])

phi_norm = curvature / Ky_actual
M_norm   = moment / My_approx

ax2.plot(phi_norm, M_norm, color='#457b9d', linewidth=2)
ax2.axvline(1.0, color='#2a9d8f', linewidth=1.0, linestyle='--', label='φ/φ_y = 1')
ax2.axhline(1.0, color='#2a9d8f', linewidth=1.0, linestyle='--')
ax2.fill_between(phi_norm, M_norm, alpha=0.12, color='#457b9d')

ax2.set_xlabel('Normalised curvature  φ / φ_y', fontsize=9)
ax2.set_ylabel('Normalised moment  M / M_y', fontsize=9)
ax2.set_title(f'Ductility  μ = φ_u / φ_y = {ductility:.1f}', fontsize=10,
              fontweight='bold')
ax2.legend(fontsize=8)
ax2.grid(True, linestyle=':', alpha=0.6)
ax2.set_xlim(left=0)
ax2.set_ylim(bottom=0)

# ── 3. Incremental stiffness  dM/dφ ─────────────────────────────────────────
ax3 = fig.add_subplot(gs[1, 1])

dM   = np.diff(moment)
dphi = np.diff(curvature)
stiffness = np.where(np.abs(dphi) > 1e-12, dM / dphi, np.nan)
phi_mid   = 0.5 * (curvature[:-1] + curvature[1:])

ax3.plot(phi_mid, stiffness, color='#e9c46a', linewidth=1.8)
ax3.axhline(0, color='black', linewidth=0.8)
ax3.set_xlabel('Curvature  φ  [1/in]', fontsize=9)
ax3.set_ylabel('Tangent stiffness  dM/dφ  [kip·in²]', fontsize=9)
ax3.set_title('Tangent Stiffness Degradation', fontsize=10, fontweight='bold')
ax3.grid(True, linestyle=':', alpha=0.6)
ax3.set_xlim(left=0)

# ── Section geometry sketch (inset on ax1) ───────────────────────────────────
ax_in = ax1.inset_axes([0.62, 0.04, 0.36, 0.30])
ax_in.set_xlim(-z1 - 1, z1 + 1)
ax_in.set_ylim(-y1 - 1, y1 + 1)
ax_in.set_aspect('equal')
ax_in.set_title('Section', fontsize=7)
ax_in.axis('off')

# Gross section
gross = plt.Rectangle((-z1, -y1), colWidth, colDepth,
                       linewidth=1.2, edgecolor='black',
                       facecolor='#d9d9d9')
ax_in.add_patch(gross)

# Core
core = plt.Rectangle((cover - z1, cover - y1),
                      colWidth - 2 * cover, colDepth - 2 * cover,
                      linewidth=1.0, edgecolor='#457b9d',
                      facecolor='#a8c8e8', alpha=0.7)
ax_in.add_patch(core)

# Rebar positions (approximate)
rebar_y  = [y1 - cover, y1 - cover, y1 - cover,
             0.0,        0.0,
            -(y1 - cover), -(y1 - cover), -(y1 - cover)]
rebar_z  = [z1 - cover, 0.0, -(z1 - cover),
             z1 - cover,     -(z1 - cover),
             z1 - cover, 0.0, -(z1 - cover)]
ax_in.plot(rebar_z, rebar_y, 'o', color='#e63946',
           markersize=4, zorder=5)

# plt.savefig('/mnt/user-data/outputs/moment_curvature.png',
#             dpi=150, bbox_inches='tight')
plt.show()
# print("Plot saved → moment_curvature.png")