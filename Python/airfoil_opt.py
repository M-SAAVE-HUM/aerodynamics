# optimization of airfoil 

# TODOS
# play around with initial guess
# change constraints to reflect init guess parameters (thickness)
# try changing TE angle
# try optimizing diff parameters (Cl, CD, Cl/Cd, stall angle?)

import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
import aerosandbox.tools.pretty_plots as p

CL_multipoint_targets = np.array([0.8, 1.0, 1.2, 1.4, 1.5, 1.6])
CL_multipoint_weights = np.array([5, 6, 7, 8, 9, 10])

Re = 500e3 * (CL_multipoint_targets / 1.25) ** -0.5
mach = 0.03

initial_guess_airfoil = asb.KulfanAirfoil("naca6412")
initial_guess_airfoil.name = "Initial Guess (NACA6412)"

opti = asb.Opti()

optimized_airfoil = asb.KulfanAirfoil(
    name="Optimized",
    lower_weights=opti.variable(
        init_guess=initial_guess_airfoil.lower_weights,
        lower_bound=-0.5,
        upper_bound=0.25,
    ),
    upper_weights=opti.variable(
        init_guess=initial_guess_airfoil.upper_weights,
        lower_bound=-0.25,
        upper_bound=0.5,
    ),
    leading_edge_weight=opti.variable(
        init_guess=initial_guess_airfoil.leading_edge_weight,
        lower_bound=-1,
        upper_bound=1,
    ),
    TE_thickness=0,
)

alpha = opti.variable(
    init_guess=np.degrees(CL_multipoint_targets / (2 * np.pi)),
    lower_bound=-5,
    upper_bound=18
)

aero = optimized_airfoil.get_aero_from_neuralfoil(
    alpha=alpha,
    Re=Re,
    mach=mach,
)

opti.subject_to([
    aero["CL"] == CL_multipoint_targets,
    np.diff(alpha) > 0,
    aero["CM"] >= -0.115,
    optimized_airfoil.local_thickness(x_over_c=0.301) >= 0.12,
    optimized_airfoil.local_thickness(x_over_c=0.90) >= 0.012,
    optimized_airfoil.TE_angle() >= 12, 
    optimized_airfoil.lower_weights[0] < -0.05,
    optimized_airfoil.upper_weights[0] > 0.05,
])

opti.subject_to(
    optimized_airfoil.local_thickness() > 0
)

get_wiggliness = lambda af: sum([
    np.sum(np.diff(np.diff(array)) ** 2)
    for array in [af.lower_weights, af.upper_weights]
])

opti.subject_to(
    get_wiggliness(optimized_airfoil) < 2 * get_wiggliness(initial_guess_airfoil)
)

opti.minimize(np.mean(aero["CD"] * CL_multipoint_weights))

sol = opti.solve(
    behavior_on_failure="return_last"
)

optimized_airfoil = sol(optimized_airfoil)
aero = sol(aero)

Re_plot = 500e3

fig, ax = plt.subplots(2, 2, figsize=(7, 8))

airfoils_and_colors = {
    "Initial Guess"           : (initial_guess_airfoil, "dimgray"),
    "NeuralFoil-Optimized"    : (optimized_airfoil, "blue"),
}

for i, (name, (af, color)) in enumerate(airfoils_and_colors.items()):
    color = p.adjust_lightness(color, 1)
    ax[0].fill(
        af.x(), af.y(),
        facecolor=(*color, 0.09),
        edgecolor=(*color, 0.6),
        linewidth=1,
        label=name,
        linestyle=(3 * i, (7, 2)),
        zorder=4 if "NeuralFoil" in name else 3,
    )
    
    aero = af.get_aero_from_neuralfoil(
        Re=Re_plot,
        mach=mach,
        alpha = np.linspace(0, 15, 41)
    )
    ax[1].plot(
        aero["CD"], aero["CL"],  #"--",
        color=color, alpha=0.7,
        label=name,
        zorder=4 if "NeuralFoil" in name else 3,
    )
    ax[2].plot(
        alpha, aero["CL"],  #"--",
        color=color, alpha=0.7,
        label=name,
        zorder=4 if "NeuralFoil" in name else 3,
    )
    ax[3].plot(
        alpha, aero["CL"] / aero["CD"],  #"--",
        color=color, alpha=0.7,
        label=name,
        zorder=4 if "NeuralFoil" in name else 3,
    )
ax[0].legend(fontsize=11, loc="lower right", ncol=len(airfoils_and_colors)//2)
ax[0].set_title("Airfoil Shapes")
ax[0].set_xlabel("$x/c$")
ax[0].set_ylabel("$y/c$")
ax[0].axis('equal')

ax[1].legend(fontsize=11, loc="lower right", ncol=len(airfoils_and_colors)//2)
ax[1].set_title(f"Aerodynamic Polars (analyzed with NeuralFoil, $\\mathrm{{Re}}=500\\mathrm{{k}}$)")
ax[1].set_xlabel("Drag Coefficient $C_D$")
ax[1].set_ylabel("Lift\nCoefficient\n$C_L$")
ax[1].set_xlim(0, 0.04)
ax[1].set_ylim(0, 1.8)
p.show_plot("Comparison of NeuralFoil-Optimized and Initial Guess Airfoils", legend=False)