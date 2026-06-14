import openseespy.opensees as ops
import matplotlib.pyplot as plt
import numpy as np

# -------------------------------------------------------------
# 1. Initialization and Model Builder
# -------------------------------------------------------------
ops.wipe()
ops.model("basic", "-ndm", 2, "-ndf", 3)

# -------------------------------------------------------------
# 2. Nodes Geometry
# -------------------------------------------------------------
ops.node(1, 0.0, 0.0)
ops.node(2, 5.0, 0.0)
ops.node(3, 5.0, 0.0)
ops.node(4, 10.0, 0.0)

# -------------------------------------------------------------
# 3. Boundary Conditions & Multipoint Constraints
# -------------------------------------------------------------
ops.fix(1, 1, 1, 1)
ops.equalDOF(2, 3, 1, 2)

# -------------------------------------------------------------
# 4. Materials and Section Properties
# -------------------------------------------------------------
E = 2.1e11
A = 116 * 1e-4
I = 4.82 * 1e-4

mat_tag = 1
k_rot = 1.0e6
ops.uniaxialMaterial("Elastic", mat_tag, k_rot)

transf_tag = 1
ops.geomTransf("Linear", transf_tag)

# -------------------------------------------------------------
# 5. Elements
# -------------------------------------------------------------
ops.element("elasticBeamColumn", 1, 1, 2, A, E, I, transf_tag)
ops.element("zeroLength", 2, 2, 3, "-mat", mat_tag, "-dir", 3)
ops.element("elasticBeamColumn", 3, 3, 4, A, E, I, transf_tag)

# -------------------------------------------------------------
# 6. Loads & Analysis Setup
# -------------------------------------------------------------
ops.timeSeries("Linear", 1)
ops.pattern("Plain", 1, 1)
ops.load(4, 0.0, -1000.0, 0.0)

ops.system("BandSPD")
ops.numberer("RCM")
ops.constraints("Transformation")
ops.integrator("LoadControl", 1.0)
ops.algorithm("Linear")
ops.analysis("Static")

# -------------------------------------------------------------
# 7. Run Analysis
# -------------------------------------------------------------
analyze_status = ops.analyze(1)

if analyze_status == 0:
    print("Analysis successfully completed!\n")

    d1 = ops.nodeDisp(1)
    d2 = ops.nodeDisp(2)
    d3 = ops.nodeDisp(3)
    d4 = ops.nodeDisp(4)

    print(f"Node 1 (Base) Displacement:   {d1}")
    print(f"Node 2 (Hinge Left) Rotation: {d2[2]:.6f} rad")
    print(f"Node 3 (Hinge Right) Rotation:{d3[2]:.6f} rad")
    print(f"Node 4 (Tip) Deflection:      {d4[1]:.6f} m")

    # -------------------------------------------------------------
    # 8. Plot
    # -------------------------------------------------------------
    # --- tweak these two values ---
    defl_scale = 10  # multiplier applied to displacements for visibility
    y_margin = 1  # fixed y-axis half-range in scaled metres (controls figure height)
    # ------------------------------

    x_orig = [0.0, 5.0, 10.0]
    y_orig = [0.0, 0.0, 0.0]

    x_mid = 5.0 + d2[0] * defl_scale
    y_mid = 0.0 + d2[1] * defl_scale

    x_def = [0.0 + d1[0] * defl_scale, 10.0 + d4[0] * defl_scale]
    y_def = [0.0 + d1[1] * defl_scale, 0.0 + d4[1] * defl_scale]

    node_x_def = [
        0.0 + d1[0] * defl_scale,
        5.0 + d2[0] * defl_scale,
        5.0 + d3[0] * defl_scale,
        10.0 + d4[0] * defl_scale,
    ]
    node_y_def = [
        0.0 + d1[1] * defl_scale,
        0.0 + d2[1] * defl_scale,
        0.0 + d3[1] * defl_scale,
        0.0 + d4[1] * defl_scale,
    ]

    fig, ax = plt.subplots(figsize=(9, 4))

    # Undeformed beam
    ax.plot(x_orig, y_orig, "k--", linewidth=1, label="Undeformed")

    # Deformed beam
    ax.plot(
        [x_def[0], x_mid, x_def[1]],
        [y_def[0], y_mid, y_def[1]],
        "b-",
        linewidth=2,
        label=f"Deformed (x{defl_scale})",
    )

    ax.plot(node_x_def, node_y_def, "bo", markersize=6)
    ax.plot(
        x_mid,
        y_mid,
        "o",
        color="orange",
        markersize=9,
        label="Spring (hinge)",
        zorder=5,
    )
    ax.plot(0.0, 0.0, "ks", markersize=8, label="Fixed support")

    # Load arrow at tip
    arrow_len = y_margin * 0.4
    ax.annotate(
        "",
        xy=(x_def[1], y_def[1]),
        xytext=(x_def[1], y_def[1] + arrow_len),
        arrowprops=dict(arrowstyle="->", color="red", lw=1.5),
    )
    ax.text(
        x_def[1] + 0.15, y_def[1] + arrow_len * 0.5, "10 kN", color="red", fontsize=9
    )

    # Tip deflection annotation
    ax.annotate(
        f"tip: {d4[1] * 1000:.2f} mm",
        xy=(x_def[1], y_def[1]),
        xytext=(x_def[1] - 2.0, y_def[1] - y_margin * 0.4),
        fontsize=8,
        color="blue",
        arrowprops=dict(arrowstyle="->", color="blue", lw=0.8),
    )

    # Fixed y limits so figure does not collapse when defl_scale is small
    ax.set_ylim(-y_margin, y_margin)

    ax.set_xlabel("x (m)")
    ax.set_ylabel(f"y (m, disp x{defl_scale})")
    ax.set_title("Cantilever beam with rotational spring - deformed shape")
    ax.legend(fontsize=8)
    ax.grid(True, linewidth=0.4, alpha=0.5)

    plt.tight_layout()
    plt.savefig("cantilever_spring.png", dpi=150)
    plt.show()
    print("Plot saved to cantilever_spring.png")

else:
    print("Analysis failed.")

ops.wipe()
