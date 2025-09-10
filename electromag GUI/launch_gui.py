#!/usr/bin/env python3
"""
Electromagnet Simulation GUI Launcher
====================================

This script launches the Electromagnet Simulation GUI.
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

# Add parent directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

try:
    import tkinter as tk
    from electromagnet_gui import ElectromagnetGUI, main
    
    print("Starting Electromagnet Simulation GUI...")
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
    print("- electromagnet_gui.py")
    print("- electromagnet.py (in parent directory)")
    print("- MWS_wire_data.csv (in parent directory)")
    print("- requirements.txt")
    input("\nPress Enter to exit...")
