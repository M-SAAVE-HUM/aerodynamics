# aircraft wing sizing optimization using aerosandbox

import aerosandbox as asb
import aerosandbox.numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

# ----------
# Parameters
# ----------

# environmental params
viscosity = 1.78e-5  # viscosity of air [kg/m/s]
density = 1.1952  # density of air [kg/m^3]
airfoil_thickness_fraction = 0.12  # airfoil thickness to chord ratio [-]

# aerodynamic params
stall_speed = 14 # stall speed [m/s]
airspeed_takeoff = stall_speed*1.1  # takeoff speed [m/s]
oswalds_efficiency = 0.95  # Oswald efficiency factor [-]
CL_max = 1.4  # max CL 
wetted_area_ratio = 2  # wetted area ratio [-]
form_factor = 1.2  # form factor [-]
#airspeed = 20 # m/s -- try setting this as constant or as a design variable

# structural params
ultimate_load_factor = 3.5  # ultimate load factor [-]
weight_fuselage = 10 * 9.81  # kg -> N

# mass budget estimate: 
# fuselage          = 3 kg         
# tail              = 1 kg            
# landing gear      = 0.5 kg
# prop system       = 2.5 kg
# avionics          = 1.5 kg 
# payload           = 0.5 kg
# margin            = 1 kg
# TOTAL             = 10 kg

# ----------
# Variables
# ----------
opti = asb.Opti()  # initialize an optimization environment

aspect_ratio = opti.variable(init_guess=6.55, lower_bound=4, upper_bound=20)  # aspect ratio
wing_area = opti.variable(init_guess=0.92, lower_bound=0.25, upper_bound=5)  # total wing area [m^2]
airspeed = opti.variable(init_guess=25, lower_bound=15, upper_bound=30)  # cruising speed [m/s]
weight = opti.variable(init_guess=100, lower_bound=20, upper_bound=500)  # total aircraft weight [N]
CL = opti.variable(init_guess=1, lower_bound=0.3, upper_bound=1.4)  # Lift coefficient of wing [-]
span = (wing_area * aspect_ratio) ** 0.5
chord = wing_area / span

# ----------
# Models
# ----------

# Aerodynamics model
Re = (density / viscosity) * airspeed * (wing_area / aspect_ratio) ** 0.5
Cf = 0.074 / Re ** 0.2
CD_profile = form_factor * Cf * wetted_area_ratio
CD_induced = CL ** 2 / (np.pi * aspect_ratio * oswalds_efficiency)
CD = CD_profile + CD_induced
dynamic_pressure = 0.5 * density * airspeed ** 2
drag = dynamic_pressure * wing_area * CD
lift_cruise = dynamic_pressure * wing_area * CL
lift_stall = 0.5 * density * wing_area * CL_max * stall_speed ** 2
lift_takeoff = 0.5 * density * wing_area * CL_max * airspeed_takeoff ** 2

# Wing weight structural model
W_W_coeff1 = 0.000521  # some number in units of 1/m, backed out using coeff_calc.py
weight_wing_structural = W_W_coeff1 * (
        ultimate_load_factor * aspect_ratio ** 1.5 *
        (weight_fuselage * weight * wing_area) ** 0.5
) / airfoil_thickness_fraction
W_W_coeff2 = 0.06  # monokote weight [kg/m^2]
weight_wing_surface = W_W_coeff2 * wing_area * 9.81  # convert kg to N
weight_wing = weight_wing_surface + weight_wing_structural

# ----------
# Constraints
# ----------
opti.subject_to([
    weight == lift_cruise, 
    weight <= lift_stall,
    weight == weight_fuselage + weight_wing, 
    span <= 2.5,
])

# ----------
# Objective
# ----------
opti.maximize(lift_cruise/drag)

# ----------
# Solve!
# ----------
sol = opti.solve(max_iter=500, verbose=True)

# ----------
# Post-process
# ----------
sol_span = (sol.value(wing_area) * sol.value(aspect_ratio)) ** 0.5
sol_chord = sol.value(wing_area) / sol_span
sol_weight_kg = sol.value(weight) / 9.81
sol_wing_weight_kg = sol.value(weight_wing) / 9.81
sol_AR = sol.value(aspect_ratio)
q = 0.5 * density * sol.value(airspeed)**2
sol_area = sol.value(wing_area)
sol_di = (sol.value(CL)**2 / (np.pi * sol_AR * oswalds_efficiency)) * q * sol_area
sol_d0 = sol.value(drag) - sol_di

print("\n" + "-"*60)
print("RESULTS")
print("-"*60)
print(f"\nGeometry:")
print(f"  Aspect ratio = {sol.value(aspect_ratio):.2f}")
print(f"  Wing area = {sol.value(wing_area):.3f} m^2")
print(f"  Span = {sol_span:.2f} m")
print(f"  Chord = {sol_chord:.3f} m (for rectangular wing)")
print(f"\nFlight:")
print(f"  C_L = {sol.value(CL):.3f}")
print(f"  L/D = {sol.value(lift_cruise/drag):.2f}")
print(f"  Cruise speed = {sol.value(airspeed):.2f} m/s")
print(f"  Total Drag = {sol.value(drag):.2f} N, CD = {sol.value(drag) / (q*sol_area)}")
print(f"  Induced Drag = {sol_di:.2f} N, CDi = {sol_di / (q*sol_area)}")
print(f"  Parasitic Drag = {sol_d0:2f} N, CDO = {sol_d0 / (q*sol_area)}")
print(f"  Takeoff lift = {sol.value(lift_takeoff):.2f} N")
print(f"  Cruise lift = {sol.value(lift_cruise):.2f} N")
print(f"  Stall lift = {sol.value(lift_stall):.2f} N")
print(f"\nWeight:")
print(f"  Total weight = {sol_weight_kg:.2f} kg ({sol.value(weight):.2f} N)")
print(f"  Wing weight = {sol_wing_weight_kg:.2f} kg ({sol.value(weight_wing):.2f} N)")
print(f"  Fuselage weight = {weight_fuselage/9.81:.2f} kg ({weight_fuselage:.2f} N)")
print(f"  Wing structural = {sol.value(weight_wing_structural):.2f} N")
print(f"  Wing surface = {sol.value(weight_wing_surface):.2f} N")
print("-"*60)

# --------------------
# extra plotting stuff
# --------------------
def wing_structural_weight(AR, S, W):
    return (
        W_W_coeff1 *
        ultimate_load_factor *
        AR**1.5 *
        np.sqrt(weight_fuselage * W * S)
        / airfoil_thickness_fraction
    )

def solve_weight(AR, S):
    def residual(W):
        W_wing_struct = wing_structural_weight(AR, S, W)
        W_wing_surface = W_W_coeff2 * S * 9.81
        return W - (weight_fuselage + W_wing_struct + W_wing_surface)

    W_guess = 50
    W_sol, = fsolve(residual, W_guess, maxfev=100)
    return W_sol

AR_vals = np.linspace(4, 16, 80)
S_vals = np.linspace(0.25, 2.5, 80)

AR_grid, S_grid = np.meshgrid(AR_vals, S_vals)

LD = np.full_like(AR_grid, np.nan)
W_struct = np.full_like(AR_grid, np.nan)
drag_ratio = np.full_like(AR_grid, np.nan)

# solve for W and L/D at each AR and S combo, if not valid then L/D = 0, W = nan
for i in range(AR_grid.shape[0]):
    for j in range(AR_grid.shape[1]):
        AR = AR_grid[i, j]
        S = S_grid[i, j]

        # Span constraint
        span = np.sqrt(AR * S)
        if span > 2.5:
            LD[i, j] = 0
            continue

        # Solve weight
        W = solve_weight(AR, S)

        # Lift constraint → CL
        CL = W / (q * S)

        if CL < 0.3 or CL > 2.0:
            LD[i, j] = 0
            continue

        # stall constraint
        L_stall = 0.5 * density * S * CL_max * stall_speed**2
        if W > L_stall:
            LD[i, j] = 0
            continue

        # Reynolds number
        Re = (density / viscosity) * sol.value(airspeed) * np.sqrt(S / AR)
        Cf = 0.074 / Re**0.2

        CD_profile = form_factor * Cf * wetted_area_ratio
        CD_induced = CL**2 / (np.pi * AR * oswalds_efficiency)
        CD = CD_profile + CD_induced

        D = q * S * CD

        LD[i, j] = W / D
        W_struct[i, j] = wing_structural_weight(AR, S, W)
        drag_ratio[i, j] = CD_induced / CD_profile

plt.figure(figsize=(8,6))
cs = plt.contourf(AR_grid, S_grid, LD, levels=30)
plt.colorbar(cs, label='L/D')
plt.plot(sol.value(aspect_ratio),sol.value(wing_area),'ro', label='Optimized')
plt.xlabel('Aspect Ratio')
plt.ylabel('Wing Area [m^2]')
plt.title('L/D Design Space')
plt.legend()
plt.show()

plt.figure(figsize=(8,6))
cs = plt.contourf(AR_grid, S_grid, W_struct / 9.81, levels=30)
plt.colorbar(cs, label='Wing Structural Weight [kg]')
plt.plot(sol.value(aspect_ratio), sol.value(wing_area),'ro')
plt.xlabel('Aspect Ratio')
plt.ylabel('Wing Area [m^2]')
plt.title('Structural Weight Design Space')
plt.show()