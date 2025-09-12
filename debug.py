import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from electromagnet import PowerSupply, Wire, Emag

# Physical constants
mu_0 = 4 * np.pi * 10**-7  # Permeability of free space

our_power = PowerSupply(name="Keysight RP7952A",
                        max_voltage=500,
                        max_current=40,
                        max_power=10000)
AWG12 = Wire(AWG=12)

helm = Emag(power_supply=our_power, wire=AWG12)

helm.add_turn(center=np.array([0, 0, 0]), radius=.125, points_per_turn=500)

'''helm.add_mirrored_turns(radius=.25,
                        separation=.25,
                        center=np.array([0, 0, 0]),
                        points_per_turn=500)'''

#helm.add_turn(center=np.array([0, 0, 0]), radius=.25, points_per_turn=500)
#helm.plot_coil_geometry()
helm.def_field(x_range=(-.25, .25), y_range=(0,0), z_range=(0, .125),
               points_per_axis=5)

helm.calc_b_field()
#print(np.linalg.norm(helm.b_field[helm.b_field.shape[0]//2]))

# Find field point closest to origin
origin = np.array([0, 0, 0])
distances = np.linalg.norm(helm.field - origin, axis=1)
origin_index = np.argmin(distances)
origin_b_field = helm.b_field[origin_index]
origin_magnitude = np.linalg.norm(origin_b_field)

print(f"Origin field magnitude: {origin_magnitude}")
print(f"Origin field vector: {origin_b_field}")

# Find field point closest to (0, 0, 0.125)
point_2 = np.array([0, 0, 0.125])
distances_2 = np.linalg.norm(helm.field - point_2, axis=1)
point_2_index = np.argmin(distances_2)
point_2_b_field = helm.b_field[point_2_index]
point_2_magnitude = np.linalg.norm(point_2_b_field)

print(f"Point (0,0,0.125) field magnitude: {point_2_magnitude}")
print(f"Point (0,0,0.125) field vector: {point_2_b_field}")

print(f"B max: {np.max(np.linalg.norm(helm.b_field, axis=1))}")

# Analytical solutions
R = 0.125  # Radius for analytical calculation
I = helm.current  # Current from the electromagnet

# Analytical solution for origin (0,0,0): B = μ₀I/(2R)
analytical_origin = mu_0 * I / (2 * R)
print(f"\nAnalytical Solutions:")
print(f"Origin (0,0,0) - B = μ₀I/(2R): {analytical_origin:.6e} T")

# Analytical solution for point (0,0,0.125): B = μ₀I*cos²(π/4)/(2R)
# Note: cos²(π/4) = 0.5
analytical_point_2 = mu_0 * I / (2**(5/2) * R)
print(f"Point (0,0,0.125) - B = μ₀I*cos²(π/4)/(2R): {analytical_point_2:.6e} T")

# Comparison with numerical results
print(f"\nComparison:")
print(f"Origin - Numerical: {origin_magnitude:.6e} T, Analytical: {analytical_origin:.6e} T")
print(f"Point (0,0,0.125) - Numerical: {point_2_magnitude:.6e} T, Analytical: {analytical_point_2:.6e} T")

# Calculate percent errors
percent_error_origin = abs(origin_magnitude - analytical_origin) / analytical_origin * 100
percent_error_point_2 = abs(point_2_magnitude - analytical_point_2) / analytical_point_2 * 100

print(f"\nPercent Errors:")
print(f"Origin: {percent_error_origin:.2f}%")
print(f"Point (0,0,0.125): {percent_error_point_2:.2f}%")

# Calculate using trapz integration method
print(f"\n--- Results using np.trapz integration ---")
b_field_trapz = helm.calc_b_field_trapz()

# Find field point closest to origin using trapz results
origin_b_field_trapz = b_field_trapz[origin_index]
origin_magnitude_trapz = np.linalg.norm(origin_b_field_trapz)

# Find field point closest to (0, 0, 0.125) using trapz results
point_2_b_field_trapz = b_field_trapz[point_2_index]
point_2_magnitude_trapz = np.linalg.norm(point_2_b_field_trapz)

print(f"Trapz Origin field magnitude: {origin_magnitude_trapz:.6e} T")
print(f"Trapz Point (0,0,0.125) field magnitude: {point_2_magnitude_trapz:.6e} T")

# Calculate percent errors for trapz method
percent_error_origin_trapz = abs(origin_magnitude_trapz - analytical_origin) / analytical_origin * 100
percent_error_point_2_trapz = abs(point_2_magnitude_trapz - analytical_point_2) / analytical_point_2 * 100

print(f"\nTrapz Percent Errors:")
print(f"Origin: {percent_error_origin_trapz:.2f}%")
print(f"Point (0,0,0.125): {percent_error_point_2_trapz:.2f}%")

# Compare sum vs trapz methods
print(f"\n--- Comparison: np.sum vs np.trapz ---")
print(f"Origin - Sum: {origin_magnitude:.6e} T, Trapz: {origin_magnitude_trapz:.6e} T")
print(f"Point (0,0,0.125) - Sum: {point_2_magnitude:.6e} T, Trapz: {point_2_magnitude_trapz:.6e} T")

percent_diff_origin = abs(origin_magnitude - origin_magnitude_trapz) / origin_magnitude * 100
percent_diff_point_2 = abs(point_2_magnitude - point_2_magnitude_trapz) / point_2_magnitude * 100

print(f"Percent difference between methods:")
print(f"Origin: {percent_diff_origin:.2f}%")
print(f"Point (0,0,0.125): {percent_diff_point_2:.2f}%")

helm.print_parameters()

# Plot with analytical solutions
analytical_points = [origin, point_2]
analytical_values = [analytical_origin, analytical_point_2]
point_labels = ['Origin (0,0,0)', 'Point (0,0,0.125)']

helm.plot_b_field_with_analytical(analytical_points=analytical_points,
                                analytical_values=analytical_values,
                                point_labels=point_labels,
                                title_suffix=" - Numerical vs Analytical")

