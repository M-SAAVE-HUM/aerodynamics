# full aircraft analysis using AeroBuildup (built in analysis tool)

# TODOS
# fix sizing
# figure out fuselage!!
# make plots
# v stab translate in z?
# lift seems too high, drag seems too low?

import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p

wing_airfoil = asb.Airfoil("naca6412")
tail_airfoil = asb.Airfoil("naca0012")

### Define the 3D geometry you want to analyze/optimize.
# Here, all distances are in meters and all angles are in degrees.
airplane = asb.Airplane(
    name="Gonk",
    xyz_ref=[0, 0, 0],  # CG location
    wings=[
        asb.Wing(
            name="Main Wing",
            symmetric=True,  # Should this wing be mirrored across the XZ plane?
            xsecs=[  # The wing's cross ("X") sections
                asb.WingXSec(  # Root
                    xyz_le=[0, 0, 0],  # Coordinates of the XSec's leading edge, relative to the wing's leading edge.
                    chord=0.5,
                    twist=0,  # degrees
                    airfoil=wing_airfoil,  # Airfoils are blended between a given XSec and the next one.
                ),
                asb.WingXSec(  # Tip
                    xyz_le=[0, 1.5, 0],
                    chord=0.5,
                    airfoil=wing_airfoil,
                ),
            ]
        ),
        asb.Wing(
            name="Horizontal Stabilizer",
            symmetric=True,
            xsecs=[
                asb.WingXSec(  # root
                    xyz_le=[0, 0, 0],
                    chord=0.321,
                    twist=0,
                    airfoil=tail_airfoil,
                ),
                asb.WingXSec(  # tip
                    xyz_le=[0, 0.4816, 0],
                    chord=0.321,
                    twist=0,
                    airfoil=tail_airfoil
                )
            ]
        ).translate([1.44, 0, 0]),
        asb.Wing(
            name="Vertical Stabilizer",
            symmetric=False,
            xsecs=[
                asb.WingXSec(
                    xyz_le=[0, 0, 0],
                    chord=0.275,
                    twist=0,
                    airfoil=tail_airfoil,
                ),
                asb.WingXSec(
                    xyz_le=[0, 0, 0.488],
                    chord=0.275,
                    twist=0,
                    airfoil=tail_airfoil
                )
            ]
        ).translate([1.44, 0, 0.03])
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

# visualize plane
#airplane.draw_three_view()

# Define single operation point

# aero = asb.AeroBuildup(
#     airplane=airplane,
#     op_point=asb.OperatingPoint(
#         velocity=28,
#         alpha=0,
#     ),
# ).run()

# multiple alphas
alpha = np.linspace(-10, 20, 300)
aero = asb.AeroBuildup(
    airplane=airplane,
    op_point=asb.OperatingPoint(
        velocity=28,
        alpha=alpha,
        beta=0
    ),
).run()

#show numbers
# print("CL = " + str(*aero["CL"]))
# print("CD = " + str(*aero["CD"]))
# print("Cm = " + str(*aero["Cm"]))
# print("Lift = " + str(*aero["L"]) + " N")
# print("Drag = " + str(*aero["D"]) + " N")


fig, ax = plt.subplots(2, 2)

plt.sca(ax[0, 0])
plt.plot(alpha, aero["CL"])
plt.xlabel(r"$\alpha$ [deg]")
plt.ylabel(r"$C_L$")
p.set_ticks(5, 1, 0.5, 0.1)

plt.sca(ax[0, 1])
plt.plot(alpha, aero["CD"])
plt.xlabel(r"$\alpha$ [deg]")
plt.ylabel(r"$C_D$")
p.set_ticks(5, 1, 0.05, 0.01)
plt.ylim(bottom=0)

plt.sca(ax[1, 0])
plt.plot(alpha, aero["Cm"])
plt.xlabel(r"$\alpha$ [deg]")
plt.ylabel(r"$C_m$")
p.set_ticks(5, 1, 0.5, 0.1)

plt.sca(ax[1, 1])
plt.plot(aero["CD"], aero["CL"])
plt.xlabel(r"$C_D$")
plt.ylabel(r"$C_L$")
p.set_ticks(5, 1, 10, 2)

p.show_plot(
    r"M-SAAVE Aircraft Aerodynamics ($\alpha$ = [-10,20])"
)

# contour plot of alpha vs beta
# Beta, Alpha = np.meshgrid(np.linspace(-90, 90, 150), np.linspace(-90, 90, 150))
# aero = asb.AeroBuildup(
#     airplane=airplane,
#     op_point=asb.OperatingPoint(
#         velocity=10,
#         alpha=Alpha.flatten(),
#         beta=Beta.flatten()
#     ),
# ).run()


# def show():
#     p.set_ticks(15, 5, 15, 5)
#     p.equal()
#     p.show_plot(
#         "`asb.AeroBuildup` Aircraft Aerodynamics",
#         r"Sideslip angle $\beta$ [deg]",
#         r"Angle of Attack $\alpha$ [deg]",
#         set_ticks=False
#     )

# fig, ax = plt.subplots(figsize=(6, 5))
# p.contour(
#     Beta, Alpha, aero["CL"].reshape(Alpha.shape),
#     colorbar_label="Lift Coefficient $C_L$ [-]",
#     linelabels_format=lambda x: f"{x:.2f}",
#     linelabels_fontsize=7,
#     cmap="RdBu",
#     alpha=0.6
# )
# plt.clim(*np.array([-1, 1]) * np.max(np.abs(aero["CL"])))
# show()

# fig, ax = plt.subplots(figsize=(6, 5))
# p.contour(
#     Beta, Alpha, aero["CD"].reshape(Alpha.shape),
#     colorbar_label="Drag Coefficient $C_D$ [-]",
#     linelabels_format=lambda x: f"{x:.2f}",
#     linelabels_fontsize=7,
#     z_log_scale=True,
#     cmap="YlOrRd",
#     alpha=0.6
# )
# show()

# fig, ax = plt.subplots(figsize=(6, 5))
# p.contour(
#     Beta, Alpha, (aero["CL"] / aero["CD"]).reshape(Alpha.shape),
#     levels=15,
#     colorbar_label="Finesse $C_L / C_D$ [-]",
#     linelabels_format=lambda x: f"{x:.0f}",
#     linelabels_fontsize=7,
#     cmap="RdBu",
#     alpha=0.6
# )
# plt.clim(*np.array([-1, 1]) * np.max(np.abs(aero["CL"] / aero["CD"])))
# show()