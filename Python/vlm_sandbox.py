# aerosandbox vortex lattice method
import aerosandbox as asb
import aerosandbox.numpy as np

wing_airfoil = asb.Airfoil("naca6412")
tail_airfoil = asb.Airfoil("naca0012")

# asb optimized - W/S = 11.52*9.81/0.628 = 132.2 N/m^2
# Geometry:
#   Aspect ratio = 9.96
#   Wing area = 0.628 m^2
#   Span = 2.50 m
#   Chord = 0.251 m (for rectangular wing)

# Flight:
#   C_L = 0.603
#   L/D = 23.35
#   Cruise speed = 22.08 m/s
#   Total Drag = 4.84 N, CD = 0.02581476955287263
#   Induced Drag = 2.29 N, CDi = 0.012228048735577483
#   Parasitic Drag = 2.547863 N, CDO = 0.01358672081729515
#   Takeoff lift = 113.04 N
#   Cruise lift = 113.04 N

# Weight:
#   Total weight = 11.52 kg (113.04 N)
#   Wing weight = 3.52 kg (34.56 N)
#   Fuselage weight = 8.00 kg (78.48 N)
#   Wing structural = 34.19 N
#   Wing surface = 0.37 N

# left W/S bound = 108.002 N/m^2 
# Geometry:
#   Aspect ratio = 6.70
#   Wing area = 0.933 m^2
#   Span = 2.50 m
#   Chord = 0.373 m (for rectangular wing)

# Flight:
#   C_L = 0.482
#   L/D = 19.64
#   Cruise speed = 19.12 m/s
#   Total Drag = 5.13 N, CD = 0.024545197394307185
#   Induced Drag = 2.43 N, CDi = 0.011626672449939448
#   Parasitic Drag = 2.699333 N, CDO = 0.01291852494436774
#   Takeoff lift = 167.98 N
#   Cruise lift = 100.75 N

# Weight:
#   Total weight = 10.27 kg (100.75 N)
#   Wing weight = 2.27 kg (22.27 N)
#   Fuselage weight = 8.00 kg (78.48 N)
#   Wing structural = 21.72 N
#   Wing surface = 0.55 N

# right W/S bound = 162.011 N/m^2
# Geometry:
#   Aspect ratio = 9.20
#   Wing area = 0.679 m^2
#   Span = 2.50 m
#   Chord = 0.272 m (for rectangular wing)

# Flight:
#   C_L = 0.577
#   L/D = 22.56
#   Cruise speed = 21.42 m/s
#   Total Drag = 4.88 N, CD = 0.025566895826090542
#   Induced Drag = 2.31 N, CDi = 0.012110634864995997
#   Parasitic Drag = 2.568040 N, CDO = 0.013456260961094544
#   Takeoff lift = 122.33 N
#   Cruise lift = 110.05 N

# Weight:
#   Total weight = 11.22 kg (110.05 N)
#   Wing weight = 3.22 kg (31.57 N)
#   Fuselage weight = 8.00 kg (78.48 N)
#   Wing structural = 31.17 N
#   Wing surface = 0.40 N

airplane = asb.Airplane(
    name="clamm",
    xyz_ref=[0, 0, 0],  # CG location
    wings=[
        asb.Wing(
            name="Main Wing",
            symmetric=True,  # Should this wing be mirrored across the XZ plane?
            xsecs=[  # The wing's cross ("X") sections
                asb.WingXSec(  # Root
                    xyz_le=[0, 0, 0],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                    chord=0.251,
                    twist=0,  # degrees
                    airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                ),
                asb.WingXSec(  # Tip
                    xyz_le=[0, 1.25, 0],
                    chord=0.251,
                    airfoil=wing_airfoil,
                ),
            ]
        )
        # asb.Wing(
        #     name="Horizontal Stabilizer",
        #     symmetric=True,
        #     xsecs=[
        #         asb.WingXSec(  # root
        #             xyz_le=[0, 0, 0],
        #             chord=0.321,
        #             twist=0,
        #             airfoil=tail_airfoil,
        #         ),
        #         asb.WingXSec(  # tip
        #             xyz_le=[0, 0.4816, 0],
        #             chord=0.321,
        #             twist=0,
        #             airfoil=tail_airfoil
        #         )
        #     ]
        # ).translate([1.44, 0, 0]),
        # asb.Wing(
        #     name="Vertical Stabilizer",
        #     symmetric=False,
        #     xsecs=[
        #         asb.WingXSec(
        #             xyz_le=[0, 0, 0],
        #             chord=0.275,
        #             twist=0,
        #             airfoil=tail_airfoil,
        #         ),
        #         asb.WingXSec(
        #             xyz_le=[0, 0, 0.488],
        #             chord=0.275,
        #             twist=0,
        #             airfoil=tail_airfoil
        #         )
        #     ]
        # ).translate([1.44, 0, 0.02])
    ],
    # fuselages=[
    #     asb.Fuselage(
    #         name="Nose",
    #         xsecs=[
    #             asb.FuselageXSec(
    #                 xyz_c=[0.8 * xi - 0.1, 0, 0.1 * xi - 0.03],
    #                 radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi)
    #             )
    #             for xi in np.cosspace(0, 1, 30)
    #         ]
    #     ),
        # asb.Fuselage(
        #     name="Box",
        #     xsecs=[
        #         asb.FuselageXSec(
        #             xyz_c=[0.8 * xi - 0.1, 0, 0.1 * xi - 0.03],
        #             radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi)
        #         )
        #         for xi in np.cosspace(0, 1, 30)
        #     ]
        # ),
        # asb.Fuselage(
        #     name="Taper",
        #     xsecs=[
        #         asb.FuselageXSec(
        #             xyz_c=[0.8 * xi - 0.1, 0, 0.1 * xi - 0.03],
        #             radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi)
        #         )
        #         for xi in np.cosspace(0, 1, 30)
        #     ]
        # ),
        # asb.Fuselage(
        #     name="TailSpar",
        #     xsecs=[
        #         asb.FuselageXSec(
        #             xyz_c=[0.8 * xi - 0.1, 0, 0.1 * xi - 0.03],
        #             radius=0.6 * asb.Airfoil("dae51").local_thickness(x_over_c=xi)
        #         )
        #         for xi in np.cosspace(0, 1, 30)
        #     ]
        # )
    #]
)

# run at set op point
vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=20.22,  # m/s
        alpha=1.107,  # degree
    )
)

L = 0; D = 0
aero = vlm.run()  # Returns a dictionary
for k, v in aero.items():
     print(f"{k.rjust(4)} : {v}")
     if k == "L":
        L = v
     if k == "D":
        D = v
LD = L/D
print(f'L/D = {LD}')
vlm.draw()

# optimize operating point
# opti = asb.Opti()

# alpha = opti.variable(init_guess=5)

# vlm = asb.VortexLatticeMethod(
#     airplane=airplane,
#     op_point=asb.OperatingPoint(
#         velocity=22.08,
#         alpha=alpha
#     ),
#     align_trailing_vortices_with_wind=False,
# )

# aero = vlm.run()

# L_over_D = aero["CL"] / aero["CD"]

# opti.maximize(L_over_D)

# opti.subject_to([
#     aero["CL"] >= 0.603, 

# ])

# sol = opti.solve()

# best_alpha = sol.value(alpha)
# print(f"Alpha for max L/D: {best_alpha:.3f} deg")