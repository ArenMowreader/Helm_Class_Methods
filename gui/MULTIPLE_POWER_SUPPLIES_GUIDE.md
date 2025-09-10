# Multiple Power Supplies Guide

This guide explains how to expand the Helmholtz GUI to support multiple power supplies, following the same pattern used in `pmax.py`.

## Current Pattern

The GUI now follows the `pmax.py` pattern for power supply handling:

```python
# Create power supply DataFrame (matching pmax.py pattern)
power_supplies = pd.DataFrame({
    'name': ['Custom'],
    'Vmax': [params['max_voltage']],
    'Imax': [params['max_current']],
    'Pmax': [params['max_power']]
})

# Select the power supply (matching pmax.py pattern)
power_supply = power_supplies.iloc[0]
```

## Adding Multiple Power Supplies

### Method 1: Add to Existing DataFrame

```python
# Add more power supplies to the existing DataFrame
power_supplies = pd.concat([power_supplies, pd.DataFrame({
    'name': ['Mean Well LRS-350-12'],
    'Vmax': [12],
    'Imax': [29.2],
    'Pmax': [350]
})], ignore_index=True)

# Add another one
power_supplies = pd.concat([power_supplies, pd.DataFrame({
    'name': ['Mean Well LRS-600-24'],
    'Vmax': [24],
    'Imax': [25],
    'Pmax': [600]
})], ignore_index=True)
```

### Method 2: Create Complete Database

```python
def create_power_supply_database():
    """Create a database of power supplies"""
    return pd.DataFrame({
        'name': [
            'Custom Power Supply',
            'Mean Well LRS-350-12', 
            'Mean Well LRS-350-24',
            'Mean Well LRS-350-48',
            'Mean Well LRS-600-12',
            'Mean Well LRS-600-24',
            'Mean Well LRS-600-48'
        ],
        'Vmax': [60, 12, 24, 48, 12, 24, 48],
        'Imax': [10, 29.2, 14.6, 7.3, 50, 25, 12.5],
        'Pmax': [600, 350, 350, 350, 600, 600, 600]
    })
```

## Analyzing Multiple Power Supplies

### Single Power Supply (Current GUI)

```python
# Select one power supply
power_supply = power_supplies.iloc[0]

# Run analysis
results_df = helm.analyze_power_supply(power_supply, ohms_per_meter[0], turns_range, params['radius'])
```

### Multiple Power Supplies

```python
# Analyze all power supplies
all_results = {}
for i, power_supply in power_supplies.iterrows():
    results_df = helm.analyze_power_supply(power_supply, ohms_per_meter[0], turns_range, params['radius'])
    all_results[power_supply['name']] = {
        'turns': turns_range,
        'current': results_df['Current (A)'].values,
        'voltage': results_df['Voltage (V)'].values,
        'power': results_df['Power (W)'].values,
        'b_field': results_df['B-Field (T)'].values,
        'power_supply': power_supply,
        'max_b_field': np.max(results_df['B-Field (T)'].values)
    }
```

## GUI Integration Options

### Option 1: Dropdown Selection

Add a dropdown to select which power supply to analyze:

```python
# In create_power_inputs method
self.power_supply_var = tk.StringVar(value="Custom")
power_supply_frame = ttk.LabelFrame(self.power_frame, text="Power Supply Selection")
power_supply_frame.pack(fill="x", padx=5, pady=5)

power_supply_label = ttk.Label(power_supply_frame, text="Power Supply:")
power_supply_label.pack(side="left", padx=5)

power_supply_combo = ttk.Combobox(power_supply_frame, textvariable=self.power_supply_var, 
                                 values=["Custom", "Mean Well LRS-350-12", "Mean Well LRS-600-24"])
power_supply_combo.pack(side="left", padx=5)
```

### Option 2: Compare All Power Supplies

Add a button to compare all power supplies:

```python
# Add compare button
compare_button = ttk.Button(self.control_frame, text="Compare All Power Supplies", 
                           command=self.compare_all_power_supplies)
compare_button.pack(side="left", padx=5)

def compare_all_power_supplies(self):
    """Compare all power supplies and show results"""
    params = self.get_simulation_parameters()
    all_results = self.analyze_multiple_power_supplies(params)
    self.display_comparison_results(all_results)
```

### Option 3: Multi-Plot Display

Show results for multiple power supplies on the same plot:

```python
def update_multi_power_supply_plots(self, all_results):
    """Update plots to show multiple power supplies"""
    self.fig.clear()
    
    # Create subplots
    gs = self.fig.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
    
    # B-Field comparison
    ax1 = self.fig.add_subplot(gs[0, 0])
    for name, data in all_results.items():
        ax1.plot(data['turns'], data['b_field'], label=name, linewidth=2)
    ax1.set_xlabel('Number of Turns')
    ax1.set_ylabel('B-Field (T)')
    ax1.set_title('B-Field vs Turns')
    ax1.legend()
    ax1.grid(True)
    
    # ... similar for other plots
```

## Example Implementation

See `multiple_power_supplies_example.py` for a complete example of how to:

1. Create a power supply database
2. Analyze multiple power supplies
3. Compare results
4. Find the best performing power supply

## Benefits of This Pattern

1. **Easy Expansion**: Just add rows to the DataFrame
2. **Consistent Interface**: Same functions work with single or multiple power supplies
3. **Flexible Analysis**: Can analyze one, some, or all power supplies
4. **Easy Comparison**: Built-in methods for comparing performance
5. **GUI Ready**: Results can be easily displayed in the GUI

## Next Steps

To fully implement multiple power supplies in the GUI:

1. Add a power supply selection dropdown
2. Create a "Compare All" button
3. Update plotting functions to handle multiple datasets
4. Add power supply database management
5. Implement result comparison views

This pattern makes it easy to expand from single power supply analysis to comprehensive multi-supply comparison and optimization. 