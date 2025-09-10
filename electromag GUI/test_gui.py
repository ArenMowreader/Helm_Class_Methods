#!/usr/bin/env python3
"""
Test script for Electromagnet GUI
================================

This script tests basic functionality of the electromagnet GUI.
"""

import sys
import os

# Add parent directory to Python path
current_dir = os.path.dirname(os.path.abspath(__file__))
if current_dir not in sys.path:
    sys.path.insert(0, current_dir)

def test_imports():
    """Test that all required modules can be imported"""
    try:
        import tkinter as tk
        print("✓ tkinter imported successfully")
        
        import numpy as np
        print("✓ numpy imported successfully")
        
        import matplotlib.pyplot as plt
        print("✓ matplotlib imported successfully")
        
        import pandas as pd
        print("✓ pandas imported successfully")
        
        # Test electromagnet classes
        sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
        from electromagnet import PowerSupply, Wire, Emag
        print("✓ electromagnet classes imported successfully")
        
        # Test GUI class
        from electromagnet_gui import ElectromagnetGUI
        print("✓ ElectromagnetGUI imported successfully")
        
        return True
        
    except ImportError as e:
        print(f"✗ Import error: {e}")
        return False
    except Exception as e:
        print(f"✗ Unexpected error: {e}")
        return False

def test_electromagnet_creation():
    """Test creating electromagnet objects"""
    try:
        from electromagnet import PowerSupply, Wire, Emag
        
        # Create power supply
        power_supply = PowerSupply(
            name="Test Supply",
            max_voltage=500,
            max_current=40,
            max_power=10000
        )
        print("✓ PowerSupply created successfully")
        
        # Create wire
        wire = Wire(AWG=12)
        print("✓ Wire created successfully")
        
        # Create magnet
        magnet = Emag(power_supply=power_supply, wire=wire)
        print("✓ Magnet created successfully")
        
        return True
        
    except Exception as e:
        print(f"✗ Error creating electromagnet objects: {e}")
        return False

def test_gui_creation():
    """Test creating GUI object"""
    try:
        import tkinter as tk
        from electromagnet_gui import ElectromagnetGUI
        
        # Create root window
        root = tk.Tk()
        root.withdraw()  # Hide the window
        
        # Create GUI
        gui = ElectromagnetGUI(root)
        print("✓ GUI created successfully")
        
        # Clean up
        root.destroy()
        
        return True
        
    except Exception as e:
        print(f"✗ Error creating GUI: {e}")
        return False

def main():
    """Run all tests"""
    print("Testing Electromagnet GUI...")
    print("=" * 40)
    
    tests = [
        ("Import Test", test_imports),
        ("Electromagnet Creation Test", test_electromagnet_creation),
        ("GUI Creation Test", test_gui_creation)
    ]
    
    passed = 0
    total = len(tests)
    
    for test_name, test_func in tests:
        print(f"\nRunning {test_name}...")
        if test_func():
            passed += 1
            print(f"✓ {test_name} PASSED")
        else:
            print(f"✗ {test_name} FAILED")
    
    print("\n" + "=" * 40)
    print(f"Tests completed: {passed}/{total} passed")
    
    if passed == total:
        print("✓ All tests passed! GUI should work correctly.")
        return True
    else:
        print("✗ Some tests failed. Please check the errors above.")
        return False

if __name__ == "__main__":
    success = main()
    if not success:
        sys.exit(1)
