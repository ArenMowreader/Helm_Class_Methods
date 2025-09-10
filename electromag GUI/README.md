# Electromagnet Simulation GUI

A graphical user interface for simulating electromagnet performance using the Biot-Savart law.

## Features

### Input Parameters
- **Coil Parameters:**
  - Inside coil radius (m)
  - Coil height (m)
  - Radial depth (m)
  - Wire gauge (AWG)

- **Power Supply:**
  - Power supply name (optional)
  - Maximum power (W)
  - Maximum current (A)
  - Maximum voltage (V)

- **Constraints:**
  - Weight limit (lbs)

- **Field Specifications:**
  - X range (min, max)
  - Y range (min, max)
  - Z range (min, max)
  - Points per axis

**Note:** Entering the same value for x, y, or z dimensions will make the analysis on that plane.

### Outputs
- **Coil Geometry Plot:** 3D visualization of the coil structure
- **B-Field Vector Plot:** Magnetic field vectors in 3D space
- **Parameter Analysis Plots:**
  - Coil Length vs Turns
  - Weight vs Turns
  - Resistance vs Turns
  - Power vs Turns
  - Current vs Turns
  - Voltage vs Turns
  - B-field Average vs Turns
  - B-field Maximum vs Turns
  - B-field Uniformity vs Turns
- **Normalized Summary Plot:** All parameters normalized for comparison

## Installation

1. Ensure you have Python 3.6 or higher installed
2. Install required packages:
   ```bash
   pip install -r requirements.txt
   ```

## Usage

1. Run the GUI:
   ```bash
   python launch_gui.py
   ```

2. Enter your parameters in the input fields:
   - Set coil dimensions (radius, height, depth)
   - Choose wire gauge
   - Configure power supply specifications
   - Set weight limit
   - Define field analysis region

3. Click "Run Simulation" to start the analysis

4. View results in the tabbed interface:
   - **Coil Geometry:** 3D coil visualization
   - **B-Field:** Magnetic field vectors
   - **Parameter Analysis:** Individual parameter plots
   - **Summary:** Normalized comparison plot
   - **Console:** Simulation progress and messages

5. Export results to CSV using the "Export Results" button

## File Structure

```
electromag GUI/
├── electromagnet_gui.py    # Main GUI application
├── launch_gui.py          # GUI launcher script
├── requirements.txt       # Python dependencies
└── README.md             # This file
```

## Dependencies

- **numpy:** Numerical computations
- **matplotlib:** Plotting and visualization
- **pandas:** Data handling and CSV export
- **tkinter:** GUI framework (usually included with Python)

## Notes

- The simulation uses the Biot-Savart law for magnetic field calculation
- Coil geometry is built layer by layer until weight limit is reached
- Field analysis can be performed on planes or in 3D space
- All calculations are performed in SI units
- Results can be exported for further analysis

## Troubleshooting

- **Import errors:** Ensure all required packages are installed
- **File not found errors:** Check that electromagnet.py and MWS_wire_data.csv are in the parent directory
- **Simulation errors:** Verify input parameters are within reasonable ranges
- **Memory issues:** Reduce points per axis for large field regions
