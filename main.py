import openseespy.opensees as ops

# -------------------------------------------------------------
# 1. Initialization and Model Builder
# -------------------------------------------------------------
ops.wipe()
ops.model('basic', '-ndm', 2, '-ndf', 3)  # 2D model, 3 DOFs per node (ux, uy, rz)

# -------------------------------------------------------------
# 2. Nodes Geometry
# -------------------------------------------------------------
# Cantilever total length = 10m, split into 2x5m segments.
# We create two nodes at the 5m mark to insert the zero-length rotational spring.

ops.node(1, 0.0, 0.0)   # Fixed Base
ops.node(2, 5.0, 0.0)   # End of Segment 1
ops.node(3, 5.0, 0.0)   # Start of Segment 2 (Coincident with Node 2)
ops.node(4, 10.0, 0.0)  # Free Tip

# -------------------------------------------------------------
# 3. Boundary Conditions & Multipoint Constraints
# -------------------------------------------------------------
# Fix Node 1 (Cantilever Base)
ops.fix(1, 1, 1, 1)

# Connect Node 2 and Node 3: 
# They must share UX and UY (rigidly fixed to each other), but RZ is free to rotate independently
# except for the spring constraint we will add.
ops.equalDOF(2, 3, 1, 2)  # Constrain UX (1) and UY (2) between node 2 and 3

# -------------------------------------------------------------
# 4. Materials and Section Properties (Elastic Beam Column)
# -------------------------------------------------------------
# Cross-section properties for the beam segments
E = 2.1e11    # Elastic Modulus (Pa)
A = 116 * 1e-4      # Cross-sectional Area (m^2)
I = 48_200 * 1e-8    # Moment of Inertia (m^4)

# Define an elastic material for our rotational spring
mat_tag = 1
k_rot = 1.0e6  # Rotational spring stiffness (Nm/rad)
ops.uniaxialMaterial('Elastic', mat_tag, k_rot)

# Linear Coordinate Transformation for the beams
transf_tag = 1
ops.geomTransf('Linear', transf_tag)

# -------------------------------------------------------------
# 5. Elements
# -------------------------------------------------------------
# Beam Segment 1 (0m to 5m)
ops.element('elasticBeamColumn', 1, 1, 2, A, E, I, transf_tag)

# Rotational Hinge (Spring) at 5m connecting Node 2 and Node 3
# It acts only on the RZ DOF (direction 3 in OpenSees 2D)
ops.element('zeroLength', 2, 2, 3, '-mat', mat_tag, '-dir', 3)

# Beam Segment 2 (5m to 10m)
ops.element('elasticBeamColumn', 3, 3, 4, A, E, I, transf_tag)

# -------------------------------------------------------------
# 6. Loads & Analysis Setup
# -------------------------------------------------------------
# Create a Plain load pattern with a Linear time series
ops.timeSeries('Linear', 1)
ops.pattern('Plain', 1, 1)

# Apply a downward point load of 10 kN at the tip (Node 4, DOF 2)
ops.load(4, 0.0, -1000.0, 0.0)

# Build the analysis solver
ops.system('BandSPD')
ops.numberer('RCM')
ops.constraints('Transformation')
ops.integrator('LoadControl', 1.0)
ops.algorithm('Linear')
ops.analysis('Static')

# -------------------------------------------------------------
# 7. Run Analysis and Print Results
# -------------------------------------------------------------
analyze_status = ops.analyze(1)

if analyze_status == 0:
    print("Analysis successfully completed!\n")
    print(f"Node 1 (Base) Displacement:   {ops.nodeDisp(1)}")
    print(f"Node 2 (Hinge Left) Rotation: {ops.nodeDisp(2)[2]:.6f} rad")
    print(f"Node 3 (Hinge Right) Rotation:{ops.nodeDisp(3)[2]:.6f} rad")
    print(f"Node 4 (Tip) Deflection:      {ops.nodeDisp(4)[1]:.6f} m")
else:
    print("Analysis failed.")

ops.wipe()