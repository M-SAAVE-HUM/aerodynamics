# aerosandbox vortex lattice method
import aerosandbox as asb
import aerosandbox.numpy as np

wing_airfoil = asb.Airfoil("naca6412")
tail_airfoil = asb.Airfoil("naca0012")

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
                    chord=0.5096, # quarter chord at 0.1274
                    twist=0,  # degrees
                    airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                ),
                asb.WingXSec(  # Tip
                    xyz_le=[0.06, 1.25, 0],
                    chord=0.2696,
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

# # run at set op point
# vlm = asb.VortexLatticeMethod(
#     airplane=airplane,
#     op_point=asb.OperatingPoint(
#         velocity=20.22,  # m/s
#         alpha=1.107,  # degree
#     )
# )

# L = 0; D = 0
# aero = vlm.run()  # Returns a dictionary
# for k, v in aero.items():
#      print(f"{k.rjust(4)} : {v}")
#      if k == "L":
#         L = v
#      if k == "D":
#         D = v
# LD = L/D
# print(f'L/D = {LD}')
# vlm.draw()

# optimize operating point
opti = asb.Opti()

alpha = opti.variable(init_guess=5)
velocity = opti.variable(init_guess=25)

vlm = asb.VortexLatticeMethod(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=velocity,
        alpha=alpha
    ),
    align_trailing_vortices_with_wind=False,
)

aero = vlm.run()

L_over_D = aero["CL"] / aero["CD"]

opti.maximize(L_over_D)

opti.subject_to([
    aero["CL"] >= 0.465, 

])

sol = opti.solve()

best_alpha = sol.value(alpha)
best_LD = sol.value(L_over_D)
best_velocity = sol.value(velocity)
print(f"Alpha for max L/D: {best_alpha:.3f} deg")
print(f"max L/D: {best_LD:.3f}")
print(f"Velocity for max L/D: {best_velocity:.3f} m/s")