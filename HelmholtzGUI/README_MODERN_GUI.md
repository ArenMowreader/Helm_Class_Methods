# Modern Helmholtz Coil Optimization GUI

A comprehensive, modern GUI for Helmholtz coil design and optimization with 5 integrated stages.

## Features

### 1. First Order Optimization
- **Required Parameters:**
  - Weight Limit (lbs): Maximum weight per coil
  - Radius (m): Coil radius
  - Max Power (W): Maximum power budget

- **Optional Parameters (leave blank to optimize):**
  - Max Current (A): Maximum current constraint
  - Max Voltage (V): Maximum voltage constraint  
  - Wire Gauge (AWG): Wire size constraint

- **Functionality:**
  - Uses `optimize_helmholtz` from `FirstOrderOptimizer.py`
  - Hybrid global+local optimization
  - Finds optimal power supply and wire specifications
  - Detailed results display with operating point

### 2. First Order Visualization
- **Performance Curves:**
  - B-field vs Turns
  - Current vs Turns
  - Voltage vs Turns
  - Power vs Turns

- **Features:**
  - Interactive matplotlib plots with zoom/pan
  - Optimal operating point markers
  - Save data as CSV
  - Save plots as PDF
  - Navigation toolbar

### 3. Geometric Optimization
- **Parameters:**
  - Automatically uses power supply and wire specs from Stage 1
  - Max Helmholtz Turns: Maximum total turns constraint
  - Weight Limit (lbs): Maximum weight constraint
  - Field Point Z (m): Target field point for optimization

- **Functionality:**
  - Uses `optimize_electromagnet_geometry` from `GeometryOptimizer.py`
  - Optimizes radial and axial layer configuration
  - Hybrid optimization with global search + local refinement
  - Detailed performance metrics

### 4. Field Visualization
- **Visualization Types:**
  - Contour plots of field magnitude
  - Vector field plots
  - Combined contour + vector plots

- **Parameters:**
  - Grid Resolution: Field calculation resolution
  - Field Range (m): Spatial extent of visualization
  - Plot Type: Choose visualization style

- **Features:**
  - Uses `electromagnet.py` methods for field calculation
  - Coil geometry overlay
  - Interactive plots with zoom/pan
  - Save plots as PDF/PNG

### 5. Design Space Visualization
- **Scan Parameters:**
  - Radial Range: Range of radial layers to scan
  - Axial Range: Range of axial layers to scan

- **Visualizations:**
  - B-field vs geometry parameters
  - Weight vs geometry parameters
  - Power vs geometry parameters
  - Valid design regions

- **Features:**
  - Uses `scan_design_space` from `GeometryOptimizer.py`
  - Constraint boundary visualization
  - Optimal point highlighting
  - Trade-off analysis

## Usage

### Quick Start
```bash
cd HelmholtzGUI
python run_modern_gui.py
```

### Workflow
1. **Start with First Order Optimization:**
   - Enter required parameters (weight, radius, max power)
   - Optionally specify constraints (current, voltage, wire gauge)
   - Click "Run Optimization"
   - Review optimal power supply and wire specifications

2. **View First Order Results:**
   - Click "Generate Plots" to see performance curves
   - Use toolbar to zoom, pan, and save plots
   - Save data as CSV for further analysis

3. **Optimize Geometry:**
   - Parameters are automatically set from Stage 1
   - Adjust additional parameters if needed
   - Click "Run Geometric Optimization"
   - Review optimal layer configuration

4. **Visualize Magnetic Field:**
   - Adjust visualization parameters
   - Choose plot type (contour, vector, or both)
   - Click "Generate Field Plot"
   - Use toolbar for interaction and saving

5. **Explore Design Space:**
   - Set scan ranges for radial and axial layers
   - Click "Scan Design Space"
   - Analyze trade-offs and constraints
   - Save comprehensive plots

## Dependencies

- Python 3.7+
- tkinter (usually included with Python)
- matplotlib
- numpy
- pandas
- scipy

## File Structure

```
HelmholtzGUI/
├── modern_helmholtz_gui.py    # Main GUI application
├── run_modern_gui.py          # Launcher script
├── README_MODERN_GUI.md       # This file
└── MWS_wire_data.csv         # Wire specifications data
```

## Integration

The GUI integrates with the existing codebase:
- `FirstOrderOptimizer.py` - First order optimization
- `GeometryOptimizer.py` - Geometric optimization and design space scanning
- `electromagnet.py` - Field calculation and visualization
- `helmholtz_utils.py` - Utility functions

## Notes

- All stages are initially disabled except the first
- Each stage enables the next when completed
- Results are passed between stages automatically
- All calculations run in background threads to keep GUI responsive
- Comprehensive error handling and user feedback
- Modern styling with clear directions for each stage

## Troubleshooting

- **Import Errors:** Ensure all required modules are in the parent directory
- **File Not Found:** Check that `MWS_wire_data.csv` is in the HelmholtzGUI folder
- **Optimization Failures:** Try adjusting parameters or constraints
- **Plot Issues:** Ensure matplotlib backend is properly configured

## Author

Based on original work by Aren Mowreader, modernized with comprehensive 5-stage workflow.
