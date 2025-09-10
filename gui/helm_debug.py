import Helmholtz_again.gui.helmholtz_utils as helm
import numpy as np

height = .3
width = .3
X, Z, grid_positions = helm.grid_positions_xz(width, height)
I = 10
radius = 0.2
wire_diameter = 0.01
turns_per_layer = 1
layers = 1
dt = 100
z_offset = 0.1
print(grid_positions.dtype)
b_vec = helm.Helm(radius, wire_diameter, turns_per_layer, layers, dt, z_offset, 
                  grid_positions, I)
print(b_vec)
