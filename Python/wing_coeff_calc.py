# calculate wing weight coefficient value for asb optimization

# Known values from gonk as reference wing
weight_wing_known = 4.7  # kg
ultimate_load_factor_known = 3.5  # typical? idk check this
aspect_ratio_known = 6.0
weight_fuselage_known = 15.0  # kg
weight_total_known = 19.7  # kg (aircraft total weight)
wing_area_known = 1.5  # m^2
airfoil_thickness_fraction_known = 0.12  # 12% thick airfoil

# Back out W_W_coeff1
# Convert weights to Newtons
g = 9.81
weight_wing_N = weight_wing_known * g
weight_fuselage_N = weight_fuselage_known * g
weight_total_N = weight_total_known * g

W_W_coeff1 = (weight_wing_N * airfoil_thickness_fraction_known) / (
    ultimate_load_factor_known * aspect_ratio_known ** 1.5 *
    (weight_fuselage_N * weight_total_N * wing_area_known) ** 0.5
)

print(f"W_W_coeff1 = {W_W_coeff1:.6f} m^-1")