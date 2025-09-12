#!/usr/bin/env python3
"""
Launcher script for the Helmholtz Coil Optimization GUI

This script launches the GUI with all 5 stages:
1. First Order Optimization
2. First Order Visualization  
3. Geometric Optimization
4. Field Visualization
5. Design Space Visualization

Author: Based on original work by Aren Mowreader
Date: 12/19/2024
"""

import sys
import os

# Add parent directory to path to import modules
current_dir = os.path.dirname(os.path.abspath(__file__))
parent_dir = os.path.dirname(current_dir)
sys.path.insert(0, parent_dir)

try:
    from helm_opt_GUI import main
    
    if __name__ == "__main__":
        print("Starting Helmholtz Coil Optimization GUI...")
        print("=" * 50)
        print("Features:")
        print("1. First Order Optimization - Find optimal power supply and wire specs")
        print("2. First Order Visualization - View performance curves")
        print("3. Geometric Optimization - Optimize coil geometry")
        print("4. Field Visualization - Plot magnetic field")
        print("5. Design Space Visualization - Explore design trade-offs")
        print("=" * 50)
        main()
        
except ImportError as e:
    print(f"Error importing modules: {e}")
    print("Make sure all required files are in the correct locations:")
    print("- helm_opt_GUI.py")
    print("- FirstOrderOptimizer.py")
    print("- GeometryOptimizer.py") 
    print("- electromagnet.py")
    print("- helmholtz_utils.py")
    print("- MWS_wire_data.csv")
    sys.exit(1)
except Exception as e:
    print(f"Error starting GUI: {e}")
    sys.exit(1)
