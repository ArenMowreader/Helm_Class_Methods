# Helmholtz Coil Simulation GUI - Portable Version

This is a portable version of the Helmholtz Coil Simulation GUI that can be easily shared and run on any computer with Python installed.

## Quick Start

1. **Install Python** (if not already installed)
   - Download from [python.org](https://python.org)
   - Make sure to check "Add Python to PATH" during installation

2. **Install Required Packages**
   ```bash
   pip install -r requirements.txt
   ```

3. **Run the GUI**
   ```bash
   python launch_gui.py
   ```

## Files Included

- `launch_gui.py` - Main launcher script (run this!)
- `helmholtz_gui.py` - The GUI application
- `helmholtz_utils.py` - Physics calculations library
- `MWS_wire_data.csv` - Wire specifications database
- `requirements.txt` - Python package dependencies
- `README_PORTABLE.md` - This file

## Features

### Two Simulation Modes:

1. **First Order Approximation** (Fast, Less Accurate)
   - Quick calculations for initial design
   - Good for parameter exploration
   - Shows B-field vs turns curves

2. **Line Integral** (Slow, More Accurate)
   - Full Biot-Savart law calculations
   - 2D field visualization
   - Topological maps and vector fields
   - Single point analysis at center

### Input Parameters:

**Coil Parameters:**
- Radius (m)
- Wire Gauge (AWG)
- Number of Turns (First Order mode)
- Turns per Layer (Line Integral mode)
- Layers (Line Integral mode)
- Coil Spacing (Line Integral mode)

**Power Supply:**
- Max Power (W)
- Max Current (A)
- Max Voltage (V)
- Power Supply Name (optional)

**Constraints:**
- Weight Limit (lbs)

### Outputs:

**First Order Mode:**
- B-Field vs Turns curve
- Current vs Turns curve
- Voltage vs Turns curve
- Power vs Turns curve
- Critical values at weight limit

**Line Integral Mode:**
- B-Field magnitude contour plot
- B-Field vector direction plot
- Single point analysis at center
- Electrical parameters summary

## Troubleshooting

### "Module not found" errors:
- Make sure you're running `launch_gui.py` from the same directory as the other files
- Check that all files are present in the same folder

### "Package not found" errors:
- Run: `pip install -r requirements.txt`
- Make sure you have Python 3.6 or higher

### GUI doesn't start:
- Check that tkinter is available (usually comes with Python)
- Try running: `python -c "import tkinter; print('tkinter OK')"`

## System Requirements

- **Python**: 3.6 or higher
- **Operating System**: Windows, macOS, or Linux
- **Memory**: At least 4GB RAM recommended
- **Storage**: ~50MB free space

## Package Dependencies

The following packages are required (automatically installed via requirements.txt):
- numpy
- pandas
- matplotlib
- tkinter (usually included with Python)

## Support

If you encounter issues:
1. Check that all files are in the same directory
2. Verify Python version: `python --version`
3. Try reinstalling packages: `pip install -r requirements.txt --force-reinstall`
4. Check the console output for error messages

## Version

This is version 1.0 of the portable Helmholtz Coil Simulation GUI. 