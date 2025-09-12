import numpy as np
import matplotlib.pyplot as plt
from FirstOrderOptimizer import optimize_helmholtz
from GeometryOptimizer import optimize_electromagnet_geometry
from electromagnet import PowerSupply, Wire
from GeometryOptimizer import GeometryOptimizer, optimize_electromagnet_geometry


#first_order = FirstOrder(radius=0.125, weight_limit=70, max_power=10000)

#result = optimize_helmholtz(radius_m=0.125, weight_limit_lbs=70, max_power_w=10000, fixed_awg=12,
#                            fixed_current=40, fixed_voltage=500)


'''geom_result = optimize_electromagnet_geometry(power_supply_specs={
                'max_voltage': 500, 'max_current': 40, 'max_power': 10000},
                wire_awg=12,
                base_radius=0.125,
                max_turns=100,
                weight_limit=5)'''

power_supply = PowerSupply("Test", 500, 40, 10000)
wire = Wire(12)

optimizer = GeometryOptimizer(
    power_supply=power_supply,
    wire=wire,
    base_radius=0.125,
    max_turns=2000,  # Smaller for faster scan
    weight_limit=70
)

scan_results = optimizer.scan_design_space(
    radial_range=(1, 2000),
    axial_range=(1, 2000),
    plot=True
)