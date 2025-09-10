import numpy as np
from electromagnet import PowerSupply, Wire, Emag

# Create a simple electromagnet setup
our_power = PowerSupply(name="Test Supply", max_voltage=100, max_current=10, max_power=1000)
AWG12 = Wire(AWG=12)
helm = Emag(power_supply=our_power, wire=AWG12)

# Add a simple coil
helm.add_turn(center=np.array([0, 0, 0]), radius=0.1, points_per_turn=100)

print("Testing different field cases:")
print("=" * 50)

# Test 1: Point case - analyze a single point
print("\n1. Point case - Single point at origin:")
helm.def_field(x_range=(0, 0), y_range=(0, 0), z_range=(0, 0))
print(f"Field shape: {helm.field.shape}")
print(f"Field points: {helm.field}")

# Test 2: 1D case - Line along x-axis
print("\n2. 1D case - Line along x-axis from -0.1 to 0.1:")
helm.def_field(x_range=(-0.1, 0.1), y_range=(0, 0), z_range=(0, 0), points_per_axis=5)
print(f"Field shape: {helm.field.shape}")
print(f"Field points (first 3): {helm.field[:3]}")

# Test 3: 1D case - Line along z-axis
print("\n3. 1D case - Line along z-axis from -0.1 to 0.1:")
helm.def_field(x_range=(0, 0), y_range=(0, 0), z_range=(-0.1, 0.1), points_per_axis=5)
print(f"Field shape: {helm.field.shape}")
print(f"Field points (first 3): {helm.field[:3]}")

# Test 4: 2D case - Plane in x-y
print("\n4. 2D case - Plane in x-y (z=0):")
helm.def_field(x_range=(-0.1, 0.1), y_range=(-0.1, 0.1), z_range=(0, 0), points_per_axis=3)
print(f"Field shape: {helm.field.shape}")
print(f"Field points (first 5): {helm.field[:5]}")

# Test 5: 3D case - Volume
print("\n5. 3D case - Volume:")
helm.def_field(x_range=(-0.05, 0.05), y_range=(-0.05, 0.05), z_range=(-0.05, 0.05), points_per_axis=2)
print(f"Field shape: {helm.field.shape}")
print(f"Field points (first 5): {helm.field[:5]}")

print("\n" + "=" * 50)
print("All field cases tested successfully!")
