import numpy as np
import matplotlib.pyplot as plt
import electromagnet as emag

inside_dia = .25 # meters
mu_0 = 4 * np.pi * 10**-7  # Permeability of free space
R = inside_dia/2 # radius

keysight = emag.PowerSupply(name="Keysight RP7952A", max_voltage=500, max_current=40, max_power=10000)
AWG12 = emag.Wire(AWG=12)

helm = emag.Emag(power_supply=keysight, wire=AWG12)

helm.def_field((0, 0), (0, 0), (0, 0))

helm.add_mirrored_turns(radius=R, separation=R, points_per_turn=1000)
helm.calc_b_field()

analytical = 8*mu_0*helm.current*helm.helm_turns/((inside_dia/2) * np.sqrt(125))


print(f"geometry 1 length: {helm.geometry_1.shape[0]}")
print(f"geometry 2 length: {helm.geometry_2.shape[0]}")
print(f"field length: {helm.field.shape[0]}")
print(f"b-field length: {helm.b_field.shape[0]}")
print(f"Analytical: {analytical:.6e} T")
print(f"B-field: {np.linalg.norm(helm.b_field[0]):.6e} T")