from electromagnet import PowerSupply, Wire, Emag
import numpy as np
import matplotlib.pyplot as plt

keysight = PowerSupply(name="Keysight RP7952A",
                        max_voltage=500,
                        max_current=40,
                        max_power=10000)
AWG12 = Wire(AWG=12)

proto = Emag(power_supply=keysight, wire=AWG12)

inside_dia = 3 * .0254 # 1.5 inches in meters
inside_box = ((2**.5)/4) * inside_dia # 1/2 of the side length of a square that fits inside coil
depth = 2 * .0254 # 2 inches in meters
turn_per_layer = int(np.floor(depth / AWG12.diameter_nom_m))
weight_limit = 40
print(f"Turns per layer: {turn_per_layer}")
print(f"Wire diameter: {AWG12.diameter_nom_m:.6f} m")
print(f"Layer depth: {depth:.6f} m")

proto.def_field(x_range=(-inside_box, inside_box), y_range=(-inside_box, inside_box), z_range=(inside_dia/2, inside_dia/2), points_per_axis=5)
current = np.zeros(turn_per_layer * turn_per_layer)
resistance = np.zeros(turn_per_layer * turn_per_layer)
voltage = np.zeros(turn_per_layer * turn_per_layer)
power = np.zeros(turn_per_layer * turn_per_layer)
length = np.zeros(turn_per_layer * turn_per_layer)
weight = np.zeros(turn_per_layer * turn_per_layer)
b_field = np.zeros((turn_per_layer * turn_per_layer, proto.field.shape[0], 3))
b_avg = np.zeros(turn_per_layer * turn_per_layer)
b_max = np.zeros(turn_per_layer * turn_per_layer)
turns = np.zeros(turn_per_layer * turn_per_layer)

# Fill the arrays layer by layer
k = 0
for j in range(2):  # For each layer
    for i in range(turn_per_layer):  # For each turn in the layer
        # Check if adding this turn would exceed weight limit
        if proto.weight >= weight_limit:
            break
            
        proto.add_turn(center=np.array([0, 0, 0 - i * AWG12.diameter_nom_m]),
                       radius= (inside_dia / 2) + (j * AWG12.diameter_nom_m),
                       points_per_turn=500)
        proto.calc_b_field()
        
        current[k] = proto.current # the ith turn in the jth layer
        weight[k] = proto.weight # the weight of the ith turn in the jth layer
        voltage[k] = proto.voltage # the voltage of the ith turn in the jth layer
        power[k] = proto.power # the power of the ith turn in the jth layer
        resistance[k] = proto.resistance # the resistance of the ith turn in the jth layer
        length[k] = proto.coil_length # the length of wire at the ith turn in the jth layer
        b_field[k] = proto.b_field # the b-field at the ith turn in the jth layer
        b_avg[k] = np.mean(np.linalg.norm(proto.b_field, axis=1)) # the average b-field at the ith turn in the jth layer
        b_max[k] = np.max(np.linalg.norm(proto.b_field, axis=1)) # the maximum b-field at the ith turn in the jth layer
        turns[k] = proto.net_turns # the number of turns at the ith turn in the jth layer
        k += 1
    
    # If we've exceeded weight limit, break out of outer loop too
    if proto.weight >= weight_limit:
        break

proto.plot_coil_geometry()

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

# Plot 7: B-field Average
axes[2, 0].plot(turns[valid_mask], b_avg[valid_mask], 'purple', linewidth=2)
axes[2, 0].set_xlabel('Number of Turns')
axes[2, 0].set_ylabel('B-field Average (T)')
axes[2, 0].set_title('B-field Average vs Turns')
axes[2, 0].grid(True)

# Plot 8: B-field Maximum
axes[2, 1].plot(turns[valid_mask], b_max[valid_mask], 'brown', linewidth=2)
axes[2, 1].set_xlabel('Number of Turns')
axes[2, 1].set_ylabel('B-field Maximum (T)')
axes[2, 1].set_title('B-field Maximum vs Turns')
axes[2, 1].grid(True)

# Plot 9: B-field Ratio (Max/Avg)
b_ratio = b_max[valid_mask] / b_avg[valid_mask]
axes[2, 2].plot(turns[valid_mask], b_ratio, 'teal', linewidth=2)
axes[2, 2].set_xlabel('Number of Turns')
axes[2, 2].set_ylabel('B-field Ratio (Max/Avg)')
axes[2, 2].set_title('B-field Uniformity vs Turns')
axes[2, 2].grid(True)

plt.tight_layout()
plt.show()

# Also create a summary plot showing all parameters on one graph (normalized)
fig2, ax2 = plt.subplots(figsize=(12, 8))

# Normalize all parameters to their maximum values for comparison
coil_length_norm = length[valid_mask] / np.max(length[valid_mask])
weight_norm = weight[valid_mask] / np.max(weight[valid_mask])
resistance_norm = resistance[valid_mask] / np.max(resistance[valid_mask])
power_norm = power[valid_mask] / np.max(power[valid_mask])
current_norm = current[valid_mask] / np.max(current[valid_mask])
voltage_norm = voltage[valid_mask] / np.max(voltage[valid_mask])
b_avg_norm = b_avg[valid_mask] / np.max(b_avg[valid_mask])
b_max_norm = b_max[valid_mask] / np.max(b_max[valid_mask])

ax2.plot(turns[valid_mask], coil_length_norm, 'b-', linewidth=2, label='Coil Length')
ax2.plot(turns[valid_mask], weight_norm, 'r-', linewidth=2, label='Weight')
ax2.plot(turns[valid_mask], resistance_norm, 'g-', linewidth=2, label='Resistance')
ax2.plot(turns[valid_mask], power_norm, 'm-', linewidth=2, label='Power')
ax2.plot(turns[valid_mask], current_norm, 'c-', linewidth=2, label='Current')
ax2.plot(turns[valid_mask], voltage_norm, 'orange', linewidth=2, label='Voltage')
ax2.plot(turns[valid_mask], b_avg_norm, 'purple', linewidth=2, label='B-field Average')
ax2.plot(turns[valid_mask], b_max_norm, 'brown', linewidth=2, label='B-field Maximum')

ax2.set_xlabel('Number of Turns')
ax2.set_ylabel('Normalized Value')
ax2.set_title('All Parameters vs Turns (Normalized)')
ax2.legend()
ax2.grid(True)
plt.show()

print(f"Analysis complete! Generated plots for {np.sum(valid_mask)} valid data points.")
print(f"Turn range: {np.min(turns[valid_mask]):.0f} to {np.max(turns[valid_mask]):.0f} turns")


