# example run of airfoil with NeuralFoil (built in airfoil analysis tool)

import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p

af = asb.Airfoil("naca6412")

alpha = np.linspace(-15, 15, 181)
re = np.geomspace(1e3, 1e8, 6)

Alpha, Re = np.meshgrid(alpha, re)

# plt, ax = plt.subplots(figsize=(6, 2))
# af.draw()

aero_flattened = af.get_aero_from_neuralfoil(
    alpha=Alpha.flatten(), # can also just be one alpha or a vector
    Re=Re.flatten(),
    mach=0.08,
    model_size="medium",  # Can be "xsmall", "small", "medium", "large", "xlarge", "xxlarge", or "xxxlarge"
    # control_surfaces=[
    #     asb.ControlSurface(
    #         name="aileron",
    #         deflection=10,  # Positive is trailing-edge down
    #         hinge_point=0.75,
    #     )
    # ],
)

# leave out rest to just get values

Aero = {
    key: value.reshape(Alpha.shape)
    for key, value in aero_flattened.items()
}

from matplotlib.colors import LinearSegmentedColormap
from aerosandbox.tools.string_formatting import eng_string

fig, ax = plt.subplots()
colors = LinearSegmentedColormap.from_list(
    "custom_cmap",
    colors=[
        p.adjust_lightness(c, 0.8) for c in
        ["orange", "darkseagreen", "dodgerblue"]
    ]
)(np.linspace(0, 1, len(re)))

for i in range(len(re)):
    line, = ax.plot(
        Aero["CD"][i, :],
        Aero["CL"][i, :],
        color=colors[i], alpha=0.8,
    )

    plt.annotate(
        f" $Re = \\mathrm{{{eng_string(re[i])}}}$",
        xy=(line.get_xdata()[-1], line.get_ydata()[-1]),
        color=colors[i],
        ha="left", va="center", fontsize=10
    )

afax = ax.inset_axes([0.76, 0.802, 0.23, 0.23])
afax.fill(
    af.x(), af.y(),
    facecolor=(0, 0, 0, 0.2), linewidth=1, edgecolor=(0, 0, 0, 0.7)
)
afax.annotate(
    text=f"{af.name} Airfoil\n",
    xy=(0.5, 0),
    ha="center", va="bottom", fontsize=10,
    alpha=0.7
)
afax.axis('off')
afax.axis('equal')

plt.xscale('log')
p.show_plot(
    title=f"$C_L$-$C_D$ Polar for {af.name} Airfoil",
    xlabel="Drag Coefficient $C_D$",
    ylabel="Lift Coefficient $C_L$",
)