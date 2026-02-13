# optimize wing taper with asb

import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p

# Wing geometry
b = 2.5        # full span [m]
S = 0.806       # planform area [m^2]
AR = b**2 / S
v_cruise = 23.01 # m/s
CL_cruise = 0.518

# helpers
def compute_oswald_efficiency_from_solution(sol, aero, AR):
    CL = sol(aero["CL"])
    CDi = sol(aero["CD"])
    return CL**2 / (np.pi * AR * CDi)

def run_lifting_line(airplane, alpha):
    return asb.LiftingLine(
        airplane=airplane,
        op_point=asb.OperatingPoint(
            velocity=v_cruise,
            alpha=alpha
        ),
        xyz_ref=[0, 0, 0],
    ).run()

########## SHAPE OPTIMIZATION ##########

# linear taper
opti = asb.Opti()  # Initialize an optimization environment.
N = 2  # Number of chord sections to optimize
section_y = np.linspace(0, b/2, N) # y locations, half span
chords = opti.variable(init_guess=np.ones(N)*0.306) # All chords initially guessed as 0.251

wing_linTaper = asb.Wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[
                -0.25 * chords[i], # This keeps the quarter-chord-line straight.
                section_y[i], # Our (known) span locations for each section.
                0
            ],
            chord=chords[i],
        )
        for i in range(N)
    ]
)

airplane_linTaper = asb.Airplane( # Make an airplane object containing only this wing.
    wings=[
        wing_linTaper
    ]
)

opti.subject_to([  # add constraints
    chords > 0.1,  # Chords should stay positive
    wing_linTaper.area() == S,  # fixed area
    np.diff(chords) <= 0, # change in chord from one section to the next should be negative
    chords[0] <= 0.75, # constrain root chord for manufacturing - 24" longest chord
])

alpha = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)

op_point = asb.OperatingPoint(
    velocity=1, # some fixed velocity - doesn't matter since we're working nondimensionally
    alpha=alpha
)

vlm = asb.VortexLatticeMethod(
    airplane=airplane_linTaper,
    op_point=op_point,
    spanwise_resolution=16,
    chordwise_resolution=8,
)

aero = vlm.run()
opti.subject_to(
    aero["CL"] == CL_cruise
)
opti.minimize(aero["CD"]) # minimize drag!
sol = opti.solve()
vlm = sol(vlm)
#vlm.draw()
#vlm.draw(show_kwargs=dict(jupyter_backend="static"))
fig, ax = plt.subplots(figsize=(4, 2))
plt.plot(
    section_y,
    sol(chords),
    ".-",
    label="AeroSandbox VLM Result",
    zorder=4,
)
ax.set_xlabel('$y$', fontsize=8)
ax.set_ylabel('Chord [m]', fontsize=8)
ax.tick_params(labelsize=8)
plt.tight_layout()

# Extract optimized chord values
optimized_chords = sol(chords)
# Get root and tip chord
root_chord = optimized_chords[0]
tip_chord = optimized_chords[1]
# Compute taper ratio
taper_ratio = tip_chord / root_chord
# compute wing area
wing_area = b*(root_chord + tip_chord)/2
# Print values
CL = sol(aero["CL"])
CDi = sol(aero["CD"])
print('-'*60)
print('LINEAR TAPER OPTIMIZATION RESULTS')
print(f'aoa: {sol(alpha):.4f}')
print(f"wing area: {wing_area:.4f}")
print(f"CL : {CL:.4f}")
print(f"CDi : {CDi:.6f}")
print(f"Root Chord: {root_chord:.4f}")
print(f"Tip Chord: {tip_chord:.4f}")
print(f"Taper Ratio: {taper_ratio:.4f}")
print('-'*60)

########## DOUBLE TAPER (fixed middle section) SHAPE OPTIMIZATION ##########

opti = asb.Opti()

# Geometry parameters
half_span = b / 2
y_break = 0.1  # constant chord region length [m]

section_y_2 = np.array([0, y_break, half_span])

# Design variables
c_root = opti.variable(init_guess=0.3)
c_tip = opti.variable(init_guess=0.15)

# Enforce constant chord inboard
chords_2 = np.array([c_root, c_root, c_tip])

wing_doubleTaper = asb.Wing(
    symmetric=True,
    xsecs=[
        asb.WingXSec(
            xyz_le=[
                -0.25 * chords_2[i],  # straight quarter-chord line
                section_y_2[i],
                0
            ],
            chord=chords_2[i],
        )
        for i in range(3)
    ]
)

airplane_doubleTaper = asb.Airplane(
    wings=[wing_doubleTaper]
)

# Constraints
opti.subject_to([
    c_root > 0.1,
    c_tip > 0.05,
    c_tip <= c_root,                # enforce taper
    c_root <= 0.6096,               # manufacturing constraint
    wing_doubleTaper.area() == S,   # FIXED AREA (use same S!)
])

alpha_2 = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)

op_point_2 = asb.OperatingPoint(
    velocity=1,
    alpha=alpha_2
)

vlm_2 = asb.VortexLatticeMethod(
    airplane=airplane_doubleTaper,
    op_point=op_point_2,
    spanwise_resolution=16,
    chordwise_resolution=8,
)

aero_2 = vlm_2.run()

opti.subject_to(aero_2["CL"] == CL_cruise)
opti.minimize(aero_2["CD"])

sol_2 = opti.solve()
vlm_2 = sol_2(vlm_2)

# Extract results
root_chord_2 = sol_2(c_root)
tip_chord_2 = sol_2(c_tip)

print('-' * 60)
print('DOUBLE TAPER OPTIMIZATION RESULTS')
print(f'aoa: {sol_2(alpha_2):.4f}')
print(f'CL : {sol_2(aero_2["CL"]):.4f}')
print(f'CDi : {sol_2(aero_2["CD"]):.6f}')
print(f'Root Chord: {root_chord_2:.4f}')
print(f'Tip Chord: {tip_chord_2:.4f}')
print(f'Taper Ratio (outer panel): {tip_chord_2/root_chord_2:.4f}')
print('-' * 60)

# Plot chord distribution
fig, ax = plt.subplots(figsize=(4, 2))
ax.plot(section_y_2, sol_2(chords_2), ".-", label="Double Taper")
ax.set_xlabel("$y$")
ax.set_ylabel("Chord [m]")
ax.grid(True)
plt.tight_layout()

########## PURE DOUBLE TAPER SHAPE OPTIMIZATION ##########

# opti = asb.Opti()

# N = 3  # double taper (3 chord control points)
# section_y_2 = np.linspace(0, 1.25, N)
# chords_2 = opti.variable(init_guess=np.ones(N) * 0.251)

# wing_doubleTaper = asb.Wing(
#     symmetric=True,
#     xsecs=[
#         asb.WingXSec(
#             xyz_le=[
#                 -0.25 * chords_2[i],  # straight quarter-chord line
#                 section_y_2[i],
#                 0
#             ],
#             chord=chords_2[i],
#         )
#         for i in range(N)
#     ]
# )

# airplane_doubleTaper = asb.Airplane(
#     wings=[wing_doubleTaper]
# )

# # Constraints
# opti.subject_to([
#     chords_2 > 0.05,
#     wing_doubleTaper.area() == 0.98,
#     np.diff(chords_2) <= 0,          # monotonic taper
#     chords_2[0] <= 0.6096,           # manufacturing constraint
# ])

# alpha_2 = opti.variable(init_guess=5, lower_bound=0, upper_bound=30)

# op_point_2 = asb.OperatingPoint(
#     velocity=1,
#     alpha=alpha_2
# )

# vlm_2 = asb.VortexLatticeMethod(
#     airplane=airplane_doubleTaper,
#     op_point=op_point_2,
#     spanwise_resolution=16,
#     chordwise_resolution=8,
# )

# aero_2 = vlm_2.run()

# opti.subject_to(aero_2["CL"] == 0.612)
# opti.minimize(aero_2["CD"])

# sol_2 = opti.solve()
# vlm_2 = sol_2(vlm_2)

# # Extract results
# opt_chords_2 = sol_2(chords_2)

# print('-' * 60)
# print('DOUBLE TAPER OPTIMIZATION RESULTS')
# print(f'aoa: {sol_2(alpha_2):.4f}')
# print(f'CL : {sol_2(aero_2["CL"]):.4f}')
# print(f'CDi : {sol_2(aero_2["CD"]):.6f}')
# print(f'Chord distribution: {opt_chords_2}')
# print('-' * 60)

# # Plot chord distribution
# fig, ax = plt.subplots(figsize=(4, 2))
# ax.plot(section_y_2, opt_chords_2, ".-", label="Double Taper")
# ax.set_xlabel("$y$")
# ax.set_ylabel("Chord [m]")
# ax.grid(True)
# plt.tight_layout()


########## WING COMPARISON ##########
# wing_airfoil = asb.Airfoil("naca6412")

# # baseline rectangular wing
# airplane1 = asb.Airplane(
#     name="Initial Design",
#     xyz_ref=[0, 0, 0],  # CG location
#     wings=[
#         asb.Wing(
#             name="Main Wing",
#             symmetric=True,  # Should this wing be mirrored across the XZ plane?
#             xsecs=[  # The wing's cross ("X") sections
#                 asb.WingXSec(  # Root
#                     xyz_le=[0, 0, 0],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
#                     chord=0.45,
#                     twist=0,  # degrees
#                     airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
#                 ),
#                 asb.WingXSec(  # Tip
#                     xyz_le=[0, 1.5, 0],
#                     chord=0.45,
#                     twist=0,
#                     airfoil=wing_airfoil,
#                 )
#             ]
#         )
#     # TODO add in tail section??
#     ]
# )

# # optimal linear taper
# airplane2 = asb.Airplane(
#     name="Linear Taper",
#     xyz_ref=[0, 0, 0],  # CG location
#     wings=[
#         asb.Wing(
#             name="Main Wing",
#             symmetric=True,  # Should this wing be mirrored across the XZ plane?
#             xsecs=[  # The wing's cross ("X") sections
#                 asb.WingXSec(  # Root
#                     xyz_le=[0, 0, 0],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
#                     chord=root_chord,
#                     twist=0,  # degrees
#                     airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
#                 ),
#                 asb.WingXSec(  # Tip
#                     xyz_le=[0, 1.5, 0],
#                     chord=tip_chord,
#                     twist=0,
#                     airfoil=wing_airfoil,
#                 )
#             ]
#         )
#     ]
# )

# # double taper wing airplane
# airplane3 = asb.Airplane(
#     name="Double Taper",
#     xyz_ref=[0, 0, 0],
#     wings=[
#         asb.Wing(
#             name="Main Wing",
#             symmetric=True,
#             xsecs=[
#                 asb.WingXSec(
#                     xyz_le=[0, section_y_2[i], 0],
#                     chord=opt_chords_2[i],
#                     twist=0,
#                     airfoil=wing_airfoil,
#                 )
#                 for i in range(N)
#             ]
#         )
#     ]
# )


# optimize operating point, init design
# opti = asb.Opti()
# alpha = opti.variable(init_guess=5)
# vlm1 = asb.VortexLatticeMethod(
#     airplane=airplane1,
#     op_point=asb.OperatingPoint(
#         velocity=20,
#         alpha=alpha
#     ),
#     align_trailing_vortices_with_wind=False,
# )
# aero_baseline = vlm1.run()
# #L_over_D = aero["CL"] / aero["CD"]
# opti.subject_to(aero_baseline["CL"] == 0.61)
# opti.minimize(aero_baseline["CD"])
# sol = opti.solve()
# best_alpha_baseline = sol(alpha)
# print('-'*60)
# print('BASELINE OPERATING POINT OPTIMIZATION')
# print(f"Alpha for CL = 0.61, init design: {best_alpha_baseline:.3f} deg") # 1.822 deg
# CL = sol(aero_baseline["CL"])
# CDi = sol(aero_baseline["CD"])
# print(f"CL initial design: {CL:.4f}")
# print(f"CDi initial design: {CDi:.4f}")
# print('-'*60)

# e_baseline = compute_oswald_efficiency_from_solution(sol, aero_baseline, AR)

# # linear taper
# opti = asb.Opti()
# alpha = opti.variable(init_guess=5)
# vlm2 = asb.VortexLatticeMethod(
#     airplane=airplane2,
#     op_point=asb.OperatingPoint(
#         velocity=20,
#         alpha=alpha
#     ),
#     align_trailing_vortices_with_wind=False,
# )
# aero_linear = vlm2.run()
# #L_over_D = aero["CL"] / aero["CD"]
# opti.subject_to(aero_linear["CL"] == 0.61)
# opti.minimize(aero_linear["CD"])
# sol2 = opti.solve()
# best_alpha_linear = sol2(alpha)
# print('-'*60)
# print('LINEAR TAPER OPTIMIZED OPERATING POINT OPTIMIZATION')
# print(f"Alpha for CL = 0.61, lin taper: {best_alpha_linear:.3f} deg") # 1.822 deg
# CL = sol2(aero_linear["CL"])
# CDi = sol2(aero_linear["CD"])
# print(f"CL lin taper: {CL:.4f}")
# print(f"CDi lin taper: {CDi:.4f}")
# print('-'*60)

# e_linear   = compute_oswald_efficiency_from_solution(sol2, aero_linear, AR)

# # double taper operating point optimization
# opti = asb.Opti()
# alpha = opti.variable(init_guess=5)

# vlm3 = asb.VortexLatticeMethod(
#     airplane=airplane3,
#     op_point=asb.OperatingPoint(
#         velocity=20,
#         alpha=alpha
#     ),
#     align_trailing_vortices_with_wind=False,
# )

# aero_double = vlm3.run()
# opti.subject_to(aero_double["CL"] == 0.61)
# opti.minimize(aero_double["CD"])

# sol3 = opti.solve()
# best_alpha_double = sol3(alpha)

# print('-' * 60)
# print('DOUBLE TAPER OPERATING POINT OPTIMIZATION')
# print(f"Alpha for CL = 0.61: {sol3(alpha):.3f} deg")
# print(f"CL double taper: {sol3(aero_double['CL']):.4f}")
# print(f"CDi double taper: {sol3(aero_double['CD']):.6f}")
# print('-' * 60)

# e_double   = compute_oswald_efficiency_from_solution(sol3, aero_double, AR)

# print("-" * 60)
# print("OSWALD EFFICIENCY COMPARISON")
# print(f"Baseline rectangular : e = {e_baseline:.4f}")
# print(f"Linear taper         : e = {e_linear:.4f}")
# print(f"Double taper         : e = {e_double:.4f}")
# print("-" * 60)


# ########## DRAG POLARS ##########

# op_point = asb.OperatingPoint(
#     atmosphere=asb.Atmosphere(altitude=0),
#     velocity=22,  # m/s
# )
# xyz_ref = [0, 0, 0]

# ll_op_point = op_point.copy()
# ll_op_point.alpha = np.linspace(-12, 12, 25)
# ll_aeros_combined = {}

# for i in range(1, 4):
#     airplane = globals()[f"airplane{i}"]
#     ll_aeros = [
#         asb.LiftingLine(
#             airplane=airplane,
#             op_point=op,
#             xyz_ref=xyz_ref,
#         ).run()
#         for op in ll_op_point
#     ]
#     ll_aero = {}
#     for k in ll_aeros[0].keys():
#         ll_aero[k] = np.array([
#             aero[k]
#             for aero in ll_aeros
#         ])
#     ll_aero["alpha"] = ll_op_point.alpha

#     ll_aeros_combined[f"airplane{i}"] = ll_aero

# # Plot CL vs. alpha
# labels = ["Baseline", "Linear Taper", "Double Taper"]
# linestyles = ["-", "-.", "--"]

# # Plot CL vs. alpha
# fig5, ax5 = plt.subplots(figsize=(4, 3))
# for i in range(1, 4):
#     data = ll_aeros_combined[f"airplane{i}"]
#     ax5.plot(data["alpha"], data["CL"], linestyle=linestyles[i-1], label=labels[i-1])
# ax5.set_xlabel("Angle of Attack $\\alpha$ [deg]")
# ax5.set_ylabel("Lift Coefficient $C_L$")
# ax5.legend()
# ax5.grid(True)
# plt.tight_layout()

# # Plot CL vs. CD
# fig6, ax6 = plt.subplots(figsize=(4, 3))
# for i in range(1, 4):
#     data = ll_aeros_combined[f"airplane{i}"]
#     ax6.plot(data["CD"], data["CL"], linestyle=linestyles[i-1], label=labels[i-1])
# ax6.set_xlabel("Drag Coefficient $C_D$")
# ax6.set_ylabel("Lift Coefficient $C_L$")
# ax6.legend()
# ax6.grid(True)
# plt.tight_layout()

plt.show()