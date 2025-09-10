# Helmholtz Coil Simulation GUI

A comprehensive GUI application for simulating Helmholtz coil configurations with two different calculation methods.

## Features

### Two Simulation Modes:
1. **First Order Approximation** - Fast computation, less accurate
2. **Line Integral** - Slower computation, more accurate

### Input Parameters:
- **Coil Parameters**: Radius, wire gauge, number of turns
- **Power Supply**: Max power, current, and voltage
- **Constraints**: Weight limit

### Outputs:
- **B-Field Curves**: Magnetic field strength vs turns
- **Power Supply Response**: Current, voltage, and power vs turns
- **Data Export**: Save results to CSV format
- **Real-time Plots**: Interactive matplotlib visualizations

## Installation

1. **Install Python Dependencies:**
   ```bash
   pip install -r requirements.txt
   ```

2. **Verify Installation:**
   ```bash
   python -c "import tkinter, matplotlib, numpy, pandas; print('All dependencies installed!')"
   ```

## Usage

### Running the GUI:

**Option 1: Direct execution**
```bash
python helmholtz_gui.py
```

**Option 2: Using the launcher**
```bash
python run_gui.py
```

### GUI Layout:

```
┌─────────────────────────────────────────────────────────────┐
│                    Simulation Mode                          │
│ ○ First Order Approximation (Fast, Less Accurate)          │
│ ○ Line Integral (Slow, More Accurate)                      │
├─────────────────────┬───────────────────────────────────────┤
│   Coil Parameters   │                                       │
│ Radius (m): [0.1]  │                                       │
│ Wire Gauge: [18]   │           Results Area                │
│ Number of Turns:   │         (Plots & Data)                │
│ [100]              │                                       │
├─────────────────────┤                                       │
│   Power Supply     │                                       │
│ Max Power: [1000]  │                                       │
│ Max Current: [10]  │                                       │
│ Max Voltage: [100] │                                       │
├─────────────────────┤                                       │
│   Constraints      │                                       │
│ Weight Limit: [5]  │                                       │
├─────────────────────┤                                       │
│ [Run] [Stop] [Save] [Clear]                                │
└─────────────────────────────────────────────────────────────┘
```

### Step-by-Step Usage:

1. **Select Simulation Mode:**
   - Choose between "First Order Approximation" or "Line Integral"
   - Input fields will automatically enable/disable based on mode

2. **Enter Parameters:**
   - **Coil Parameters**: Set radius, wire gauge, and turns configuration
   - **Power Supply**: Define maximum power, current, and voltage limits
   - **Constraints**: Set weight limit for the coil

3. **Run Simulation:**
   - Click "Run Simulation" to start computation
   - Progress bar shows simulation status
   - Results appear in the right panel

4. **View Results:**
   - **Plots Tab**: Interactive matplotlib plots
   - **Data Tab**: Numerical results and summary

5. **Save Results:**
   - Click "Save Results" to export data to CSV
   - Choose save location and filename

## Input Validation

The GUI includes comprehensive input validation:
- All numeric values must be positive
- Wire gauge must be between 0-50 AWG
- Automatic error messages for invalid inputs

## Technical Details

### Architecture:
- **Frontend**: Tkinter with ttk widgets
- **Plotting**: Matplotlib with TkAgg backend
- **Data Handling**: NumPy and Pandas
- **Threading**: Background simulation execution

### File Structure:
```
gui/
├── helmholtz_gui.py      # Main GUI application
├── run_gui.py           # Launcher script
├── requirements.txt     # Python dependencies
└── README.md           # This file
```

### Integration:
The GUI is designed to integrate with your existing `helmholtz_utils` module. Replace the placeholder simulation functions with actual calls to your calculation modules.

## Customization

### Adding New Parameters:
1. Add new variables in `create_*_inputs()` methods
2. Update `get_simulation_parameters()` method
3. Add validation in `validate_inputs()` method

### Adding New Plots:
1. Modify `update_plots()` method
2. Add new subplot to the gridspec
3. Update results data structure

### Styling:
- Modify colors, fonts, and layout in the respective creation methods
- Use ttk themes for consistent appearance

## Troubleshooting

### Common Issues:

1. **Import Errors:**
   - Ensure all dependencies are installed
   - Check Python path includes parent directory

2. **GUI Not Starting:**
   - Verify tkinter is available: `python -c "import tkinter"`
   - Check for display issues on headless systems

3. **Plot Display Issues:**
   - Ensure matplotlib backend is properly configured
   - Check for display environment variables

4. **Performance Issues:**
   - Use "First Order Approximation" for quick testing
   - Reduce number of data points for faster plotting

## Future Enhancements

- [ ] Add 3D visualization capabilities
- [ ] Implement parameter presets/saved configurations
- [ ] Add batch processing for multiple simulations
- [ ] Include more detailed error analysis
- [ ] Add export to different file formats
- [ ] Implement real-time parameter adjustment

## Contributing

When modifying the GUI:
1. Maintain the existing structure and naming conventions
2. Add proper error handling for new features
3. Update documentation for new functionality
4. Test with both simulation modes

## License

This GUI is part of the Helmholtz Coil Design-Build Thesis project. 