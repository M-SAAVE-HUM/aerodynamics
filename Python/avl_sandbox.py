import aerosandbox as asb
import aerosandbox.numpy as np

wing_airfoil = asb.Airfoil("naca6412")
tail_airfoil = asb.Airfoil("naca0012")

### Define the 3D geometry you want to analyze/optimize.
# Here, all distances are in meters and all angles are in degrees.
gonk = asb.Airplane(
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

analysis = asb.AVL(
    airplane=gonk,
    op_point=asb.OperatingPoint(
        atmosphere=asb.Atmosphere(altitude=0),
        velocity=28,
        alpha=0,
        beta=0,
    ),
    avl_command="/Users/kabirkhwaja/Desktop/College/msaave/avl/avl/",  # If AVL is not on your PATH, then set avl_command to the filepath to your AVL executable.
)

analysis.write_avl("/Users/kabirkhwaja/Desktop/College/msaave/avl/runs/gonk_sandbox.avl")

output = analysis.run()

from pprint import pprint

pprint(output)