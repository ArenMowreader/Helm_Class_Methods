from electromagnet import PowerSupply, Wire, Emag
import numpy as np
import matplotlib.pyplot as plt

# Physical constants
mu_0 = 4 * np.pi * 10**-7  # Permeability of free space

keysight = PowerSupply(name="Keysight RP7952A",
                        max_voltage=500,
                        max_current=40,
                        max_power=10000)
AWG12 = Wire(AWG=12)

scale = Emag(power_supply=keysight, wire=AWG12)

inside_dia = .25 # meters
R = inside_dia/2 # radius
field_range = (inside_dia/2) * .9 # meters
#inside_box = ((2**.5)/4) * inside_dia # 1/2 of the side length of a square that fits inside coil
depth = 2 * .0254 # 2 inches in meters
turn_per_layer = 10
layers = 120
weight_limit = 140 # lbs
print(f"Turns per layer: {turn_per_layer}")
print(f"Wire diameter: {AWG12.diameter_nom_m:.6f} m")
print(f"Layer depth: {depth:.6f} m")

#scale.def_field(x_range=(-inside_box, inside_box), y_range=(-inside_box, inside_box), z_range=(inside_dia/2, inside_dia/2), points_per_axis=5)
scale.def_field((0, 0), (0, 0), (0, 0))
current = np.zeros(turn_per_layer * layers)
resistance = np.zeros(turn_per_layer * layers)
voltage = np.zeros(turn_per_layer * layers)
power = np.zeros(turn_per_layer * layers)
length = np.zeros(turn_per_layer * layers)
weight = np.zeros(turn_per_layer * layers)
b_field = np.zeros((turn_per_layer * layers, scale.field.shape[0], 3))
b_field_avg = np.zeros(turn_per_layer * layers)
b_field_origin = np.zeros(turn_per_layer * layers)
b_mag = np.zeros(turn_per_layer * layers)
turns = np.zeros(turn_per_layer * layers)
# Fill the arrays layer by layer
k = 0
for j in range(layers):  # for each layer
    for i in range(turn_per_layer ):  # For each turn in the layer
        # Check if adding this turn would exceed weight limit
        if scale.weight >= weight_limit:
            break
            
        scale.add_mirrored_turns(radius= R + (j * AWG12.diameter_nom_m),
                                separation= R + (2 * i * AWG12.diameter_nom_m), #2* because of mirrored turns
                                points_per_turn=500)
        scale.calc_b_field()
        
        current[k] = scale.current # the ith turn in the jth layer
        weight[k] = scale.weight # the weight of the ith turn in the jth layer
        voltage[k] = scale.voltage # the voltage of the ith turn in the jth layer
        power[k] = scale.power # the power of the ith turn in the jth layer
        resistance[k] = scale.resistance # the resistance of the ith turn in the jth layer
        length[k] = scale.coil_length # the length of wire at the ith turn in the jth layer
        b_field[k] = scale.b_field # the b-field at the ith turn in the jth layer
        #b_field_avg[k] = np.mean(np.linalg.norm(scale.b_field, axis=1)) # the b-field magnitude at the ith turn in the jth layer
        #b_field_origin[k] = np.linalg.norm(scale.b_field[scale.b_field.shape[0]//2]) # the b-field magnitude at the origin
        b_mag[k] = np.linalg.norm(scale.b_field[0]) # the b-field magnitude at the ith turn in the jth layer
        turns[k] = scale.helm_turns # the number of turns at the ith turn in the jth layer
        k += 1
    print("% complete: ", (k / (turn_per_layer * layers)) * 100)
    # Break when weight limit is exceeded
    if scale.weight >= weight_limit:
        break
    
#final field analysis at weight limit

scale.print_parameters()
print("b-field origin: ", scale.b_field[0])
print("b origin magnitude: ", np.linalg.norm(scale.b_field[0]))

#redefine field
scale.def_field ((0, 0), (0, 0), (0, 0), points_per_axis=5)
scale.calc_b_field()
scale.plot_b_field()
print("b-field max: ", np.max(scale.b_field))
print("b-field avg: ", np.mean(np.linalg.norm(scale.b_field, axis=1)))

analytical = 8*mu_0*current*turns/( R * np.sqrt(125))

# Create figure with subplots for each parameter
fig, axes = plt.subplots(3, 3, figsize=(15, 12))
fig.suptitle('Prototype Projections', fontsize=16)

# Remove any zero or invalid values
valid_mask = turns > 0

# Plot 1: Coil Length
axes[0, 0].plot(turns[valid_mask], length[valid_mask], 'b-', linewidth=2)
axes[0, 0].set_xlabel('Number of Turns')
axes[0, 0].set_ylabel('Coil Length (m)')
axes[0, 0].set_title('Coil Length vs Turns')
axes[0, 0].grid(True)

# Plot 2: Weight
axes[0, 1].plot(turns[valid_mask], weight[valid_mask], 'r-', linewidth=2)
axes[0, 1].set_xlabel('Number of Turns')
axes[0, 1].set_ylabel('Weight (lbs)')
axes[0, 1].set_title('Weight vs Turns')
axes[0, 1].grid(True)

# Plot 3: Resistance
axes[0, 2].plot(turns[valid_mask], resistance[valid_mask], 'g-', linewidth=2)
axes[0, 2].set_xlabel('Number of Turns')
axes[0, 2].set_ylabel('Resistance (Ω)')
axes[0, 2].set_title('Resistance vs Turns')
axes[0, 2].grid(True)

# Plot 4: Power
axes[1, 0].plot(turns[valid_mask], power[valid_mask], 'm-', linewidth=2)
axes[1, 0].set_xlabel('Number of Turns')
axes[1, 0].set_ylabel('Power (W)')
axes[1, 0].set_title('Power vs Turns')
axes[1, 0].grid(True)

# Plot 5: Current
axes[1, 1].plot(turns[valid_mask], current[valid_mask], 'c-', linewidth=2)
axes[1, 1].set_xlabel('Number of Turns')
axes[1, 1].set_ylabel('Current (A)')
axes[1, 1].set_title('Current vs Turns')
axes[1, 1].grid(True)

# Plot 6: Voltage
axes[1, 2].plot(turns[valid_mask], voltage[valid_mask], 'orange', linewidth=2)
axes[1, 2].set_xlabel('Number of Turns')
axes[1, 2].set_ylabel('Voltage (V)')
axes[1, 2].set_title('Voltage vs Turns')
axes[1, 2].grid(True)

# Plot 7: B-field Magnitude at Origin
axes[2, 0].plot(turns[valid_mask], b_mag[valid_mask], 'purple', linewidth=2)
axes[2, 0].plot(turns[valid_mask], analytical[valid_mask], 'green', linewidth=2)
axes[2, 0].set_xlabel('Number of Turns')
axes[2, 0].set_ylabel('B-field Magnitude at Origin (T)')
axes[2, 0].set_title('B-field Magnitude at Origin vs Turns')
axes[2, 0].grid(True)

# Plot 8: B-field Average and Magnitude at Origin
'''axes[2, 1].plot(turns[valid_mask], b_field_avg[valid_mask], 'brown', linewidth=2)
axes[2, 1].plot(turns[valid_mask], b_field_origin[valid_mask], 'teal', linewidth=2)
axes[2, 1].set_xlabel('Number of Turns')
axes[2, 1].set_ylabel('B-field Average and Magnitude at Origin (T)')
axes[2, 1].set_title('B-field Average and Magnitude at Origin vs Turns')
axes[2, 1].grid(True)
'''


plt.tight_layout()
plt.show()



print(f"Analysis complete! Generated plots for {np.sum(valid_mask)} valid data points.")
print(f"Turn range: {np.min(turns[valid_mask]):.0f} to {np.max(turns[valid_mask]):.0f} turns")


