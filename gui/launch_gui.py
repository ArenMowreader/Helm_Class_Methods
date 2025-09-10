#!/usr/bin/env python3
"""
Helmholtz Coil Simulation GUI Launcher
=====================================

This script launches the Helmholtz Coil Simulation GUI.
Simply run this file to start the application.

Requirements:
- Python 3.6 or higher
- Required packages listed in requirements.txt

To install requirements:
    pip install -r requirements.txt

To run:
    python launch_gui.py
"""

import sys
import os

# Add current directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    import tkinter as tk
    from helmholtz_gui import HelmholtzGUI, main
    
    print("Starting Helmholtz Coil Simulation GUI...")
    print("Loading...")
    
    # Launch the GUI
    main()
    
except ImportError as e:
    print(f"Error: Missing required package - {e}")
    print("\nPlease install the required packages:")
    print("pip install -r requirements.txt")
    input("\nPress Enter to exit...")
    
except Exception as e:
    print(f"Error starting GUI: {e}")
    print("\nPlease check that all files are present:")
    print("- helmholtz_gui.py")
    print("- helmholtz_utils.py") 
    print("- MWS_wire_data.csv")
    print("- requirements.txt")
    input("\nPress Enter to exit...") 