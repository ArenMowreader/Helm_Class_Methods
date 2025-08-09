import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
from electromagnet import PowerSupply, Wire, Magnet


our_power = PowerSupply(name="Keysight RP7952A",
                        max_voltage=500,
                        max_current=40,
                        max_power=10000)
AWG12 = Wire(AWG=12)

helm = Magnet(power_supply=our_power, wire=AWG12)

helm.add_turn(center=np.array([0, 0, 0]),
            radius=.25,
            geometry_index=1,
            points_per_turn=500)

helm.add_turn(center=np.array([0, 0, 0.1]),
            radius=.25,
            geometry_index=2,  # Add to geometry_2
            points_per_turn=750)
helm.add_turn(center=np.array([0, 0, 0.2]),
            radius=.25,
            geometry_index=1,
            points_per_turn=500)

helm.add_turn(center=np.array([0, 0, 0.3]),
            radius=.25,
            geometry_index=2,  # Add to geometry_2
            points_per_turn=750)

helm.print_parameters()
helm.get_geometry_info()
print(np.diff(helm.geometry_1, axis=0))
helm.plot_coil_geometry()

helm.def_field(x_range=(-1, 1), y_range=(-1, 1), z_range=(-1, 1))