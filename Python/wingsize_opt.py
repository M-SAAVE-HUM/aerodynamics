# aircraft wing sizing optimization using aerosandbox

# TODOS
# fix constants to align with gonk
# ultimate load factor, wetted area ratio, fuselage drag area, CL_max

import aerosandbox as asb
import aerosandbox.numpy as np

### Constants
form_factor = 1.2  # form factor [-]
oswalds_efficiency = 0.95  # Oswald efficiency factor [-]
viscosity = 1.78e-5  # viscosity of air [kg/m/s]
density = 1.23  # density of air [kg/m^3]
airfoil_thickness_fraction = 0.12  # airfoil thickness to chord ratio [-]
ultimate_load_factor = 3.8  # ultimate load factor [-]
airspeed_takeoff = 14  # takeoff speed [m/s]
CL_max = 1.25  # max CL 
wetted_area_ratio = 2.05  # wetted area ratio [-]
drag_area_fuselage = 0.031  # fuselage drag area [m^2]
weight_fuselage = 200  # aircraft weight excluding wing [N]
weight_wing = 50 # N

opti = asb.Opti()  # initialize an optimization environment

### Variables
aspect_ratio = opti.variable(init_guess=10)  # aspect ratio
wing_area = opti.variable(init_guess=10)  # total wing area [m^2]
airspeed = opti.variable(init_guess=100)  # cruising speed [m/s]
weight = opti.variable(init_guess=10000)  # total aircraft weight [N]
CL = opti.variable(init_guess=1)  # Lift coefficient of wing [-]
chord = (wing_area / aspect_ratio) ** 0.5 # chord [m]
span = wing_area / chord # wingspan [m]

### Models
# Aerodynamics model
CD_fuselage = drag_area_fuselage / wing_area
Re = (density / viscosity) * airspeed * (wing_area / aspect_ratio) ** 0.5
Cf = 0.074 / Re ** 0.2
CD_profile = form_factor * Cf * wetted_area_ratio
CD_induced = CL ** 2 / (np.pi * aspect_ratio * oswalds_efficiency)
CD = CD_fuselage + CD_profile + CD_induced
dynamic_pressure = 0.5 * density * airspeed ** 2
drag = dynamic_pressure * wing_area * CD
lift_cruise = dynamic_pressure * wing_area * CL
lift_takeoff = 0.5 * density * wing_area * CL_max * airspeed_takeoff ** 2

# Wing weight model
# weight_wing_structural = W_W_coeff1 * (
#         ultimate_load_factor * aspect_ratio ** 1.5 *
#         (weight_fuselage * weight * wing_area) ** 0.5
# ) / airfoil_thickness_fraction
# weight_wing_surface = W_W_coeff2 * wing_area
# weight_wing = weight_wing_surface + weight_wing_structural

### Constraints
opti.subject_to([
    weight <= lift_cruise,
    weight <= lift_takeoff,
    weight == weight_fuselage + weight_wing,
    span <= 3
])

# Objective
opti.maximize(lift_cruise/drag)

sol = opti.solve(max_iter=100)

chord = ( sol.value(wing_area) / sol.value(aspect_ratio) ) ** 0.5
span = chord * sol.value(aspect_ratio)

print(f"Minimum drag = {sol.value(drag)} N")
print(f"Takeoff lift = {sol.value(lift_takeoff)} N")
print(f"Cruise lift = {sol.value(lift_cruise)} N")
print(f"Aspect ratio = {sol.value(aspect_ratio)}")
print(f"Wing area = {sol.value(wing_area)} m^2")
print(f"Span = " + str(span) +  " m")
print(f"Chord = " + str(chord) +  " m")
print(f"Airspeed = {sol.value(airspeed)} m/s")
print(f"Weight = {sol.value(weight)} N")
print(f"C_L = {sol.value(CL)}")
print(f"L/D = {sol.value(lift_cruise/drag)}")