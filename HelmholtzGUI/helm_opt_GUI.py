'''
Helmholtz Coil Optimization GUI

A comprehensive GUI for Helmholtz coil design with 4 main stages:
1. First Order Optimization
2. Geometric Optimization
3. Field Visualization
4. Design Space Visualization

Author: Based on original work by Aren Mowreader
Date: 12/19/2024
'''

import tkinter as tk
from tkinter import ttk, messagebox, filedialog
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg, NavigationToolbar2Tk
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
import threading
import os
import sys

# Add parent directory to path to import modules
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

from FirstOrderOptimizer import optimize_helmholtz
from GeometryOptimizer import optimize_electromagnet_geometry, GeometryOptimizer
from electromagnet import PowerSupply, Wire, Emag
import helmholtz_utils as helm

class HelmholtzGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Helmholtz Coil Optimization Suite")
        self.root.geometry("1400x900")
        
        # Style configuration
        self.setup_styles()
        
        # Data storage
        self.first_order_results = None
        self.geometric_results = None
        self.field_results = None
        self.design_space_results = None
        
        # Initialize all variables for both tabs
        self.initialize_variables()
        
        # Create main interface
        self.create_main_interface()
    
    def initialize_variables(self):
        """Initialize all tkinter variables for both tabs"""
        # First Order Optimization variables
        self.weight_var = tk.DoubleVar(value=75.0)
        self.radius_var = tk.DoubleVar(value=0.125)
        self.max_power_var = tk.DoubleVar(value=10000.0)
        self.max_current_var = tk.StringVar(value="")
        self.max_voltage_var = tk.StringVar(value="")
        self.wire_gauge_var = tk.StringVar(value="")
        
        # Geometric Optimization variables
        self.geo_weight_var = tk.DoubleVar(value=75.0)
        self.geo_radius_var = tk.DoubleVar(value=0.125)
        self.geo_max_power_var = tk.DoubleVar(value=10000.0)
        self.geo_max_current_var = tk.DoubleVar(value=100.0)
        self.geo_max_voltage_var = tk.DoubleVar(value=100.0)
        self.geo_wire_gauge_var = tk.DoubleVar(value=14.0)
        
    def setup_styles(self):
        """Configure classic styling"""
        style = ttk.Style()
        style.theme_use('default')
        
        # Configure colors - classic theme
        style.configure('Title.TLabel', font=('Arial', 14, 'bold'))
        style.configure('Stage.TLabel', font=('Arial', 11, 'bold'))
        style.configure('Info.TLabel', font=('Arial', 9))
        style.configure('Success.TLabel', font=('Arial', 9, 'bold'))
        style.configure('Warning.TLabel', font=('Arial', 9, 'bold'))
        
        # Configure buttons - classic theme
        style.configure('Action.TButton', font=('Arial', 9), padding=(8, 4))
        style.configure('Primary.TButton', font=('Arial', 9, 'bold'), padding=(8, 4))
        
    def create_main_interface(self):
        """Create the main interface with 5 stages"""
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        # Title
        title_label = ttk.Label(main_frame, text="Helmholtz Coil Optimization Suite", 
                               style='Title.TLabel')
        title_label.pack(pady=(0, 15))
        
        # Create notebook for stages
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill="both", expand=True)
        
        # Create each stage
        self.create_first_order_optimization_tab()
        self.create_geometric_optimization_tab()
        
    def create_first_order_optimization_tab(self):
        """Create First Order Optimization tab"""
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text="1. First Order Optimization")
        
        # Main content frame
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Left panel - Input parameters
        left_panel = ttk.LabelFrame(content_frame, text="Optimization Parameters", padding="10")
        left_panel.pack(side="left", fill="y", padx=(0, 10))
        
        # Required parameters
        required_frame = ttk.LabelFrame(left_panel, text="Required Parameters", padding="10")
        required_frame.pack(fill="x", pady=(0, 15))
        
        # Weight limit
        ttk.Label(required_frame, text="Weight Limit (lbs):").grid(row=0, column=0, sticky="w", pady=5)
        weight_entry = ttk.Entry(required_frame, textvariable=self.weight_var, width=15)
        weight_entry.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Radius
        ttk.Label(required_frame, text="Radius (m):").grid(row=1, column=0, sticky="w", pady=5)
        radius_entry = ttk.Entry(required_frame, textvariable=self.radius_var, width=15)
        radius_entry.grid(row=1, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Max Power
        ttk.Label(required_frame, text="Max Power (W):").grid(row=2, column=0, sticky="w", pady=5)
        max_power_entry = ttk.Entry(required_frame, textvariable=self.max_power_var, width=15)
        max_power_entry.grid(row=2, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Optional parameters
        optional_frame = ttk.LabelFrame(left_panel, text="Optional Parameters (Leave blank to optimize)", padding="10")
        optional_frame.pack(fill="x", pady=(0, 10))
        
        # Max Current
        ttk.Label(optional_frame, text="Max Current (A):").grid(row=0, column=0, sticky="w", pady=5)
        current_entry = ttk.Entry(optional_frame, textvariable=self.max_current_var, width=15)
        current_entry.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Max Voltage
        ttk.Label(optional_frame, text="Max Voltage (V):").grid(row=1, column=0, sticky="w", pady=5)
        voltage_entry = ttk.Entry(optional_frame, textvariable=self.max_voltage_var, width=15)
        voltage_entry.grid(row=1, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Wire Gauge
        ttk.Label(optional_frame, text="Wire Gauge (AWG):").grid(row=2, column=0, sticky="w", pady=5)
        wire_gauge_entry = ttk.Entry(optional_frame, textvariable=self.wire_gauge_var, width=15)
        wire_gauge_entry.grid(row=2, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Optimization note
        note_text = "Note: Leave optional parameters blank to optimize over them.\nFill them in to fix those parameters during optimization."
        note_label = ttk.Label(optional_frame, text=note_text, style='Info.TLabel', wraplength=200)
        note_label.grid(row=3, column=0, columnspan=2, sticky="w", pady=(5, 0))
        
        # Control buttons
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill="x", pady=(10, 0))
        
        self.optimize_button = ttk.Button(button_frame, text="Run Optimization", 
                                        command=self.run_first_order_optimization,
                                        style='Primary.TButton')
        self.optimize_button.pack(side="left", padx=(0, 10))
        
        self.clear_button = ttk.Button(button_frame, text="Clear Results", 
                                     command=self.clear_first_order_results,
                                     style='Action.TButton')
        self.clear_button.pack(side="left", padx=(0, 10))
        
        self.sync_to_geo_button = ttk.Button(button_frame, text="Sync to Geometric", 
                                           command=self.sync_first_order_to_geometric,
                                           style='Action.TButton')
        self.sync_to_geo_button.pack(side="left")
        
        # Right panel - Results and directions
        right_panel = ttk.Frame(content_frame)
        right_panel.pack(side="right", fill="both", expand=True)
        
        # Directions
        directions_frame = ttk.LabelFrame(right_panel, text="Directions", padding="10")
        directions_frame.pack(fill="x", pady=(0, 10))
        
        directions_text = """First Order Optimization Stage:

1. Enter required parameters:
   • Weight Limit: Maximum weight per coil in pounds
   • Radius: Coil radius in meters  
   • Max Power: Maximum power budget in watts

2. Optionally specify constraints:
   • Max Current: Leave blank to optimize
   • Max Voltage: Leave blank to optimize
   • Wire Gauge: Leave blank to optimize
   """
        
        directions_label = ttk.Label(directions_frame, text=directions_text, style='Info.TLabel', 
                                   justify="left", wraplength=400)
        directions_label.pack(anchor="w")
        
        # Results area
        results_frame = ttk.LabelFrame(right_panel, text="Optimization Results", padding="10")
        results_frame.pack(fill="both", expand=True)
        
        # Results text area
        self.results_text = tk.Text(results_frame, height=20, width=60, wrap="word", 
                                   font=("Courier", 9))
        results_scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.results_text.yview)
        self.results_text.configure(yscrollcommand=results_scrollbar.set)
        
        self.results_text.pack(side="left", fill="both", expand=True)
        results_scrollbar.pack(side="right", fill="y")
        
    def run_first_order_optimization(self):
        """Run first order optimization"""
        try:
            # Validate inputs
            if not self.validate_first_order_inputs():
                return
                
            # Get parameters
            weight = self.weight_var.get()
            radius = self.radius_var.get()
            max_power = self.max_power_var.get()
            
            # Parse optional parameters
            max_current = None
            max_voltage = None
            wire_gauge = None
            
            if self.max_current_var.get().strip():
                max_current = float(self.max_current_var.get())
            if self.max_voltage_var.get().strip():
                max_voltage = float(self.max_voltage_var.get())
            if self.wire_gauge_var.get().strip():
                wire_gauge = float(self.wire_gauge_var.get())
            
            # Update status
            self.results_text.delete(1.0, tk.END)
            self.results_text.insert(tk.END, "Running optimization...\n")
            self.results_text.insert(tk.END, f"Weight limit: {weight} lbs\n")
            self.results_text.insert(tk.END, f"Radius: {radius} m\n")
            self.results_text.insert(tk.END, f"Max power: {max_power} W\n")
            self.results_text.insert(tk.END, "\nOptimizing...\n")
            self.root.update()
            
            # Run optimization in separate thread
            thread = threading.Thread(target=self._run_optimization_thread, 
                                    args=(weight, radius, max_power, max_current, max_voltage, wire_gauge))
            thread.daemon = True
            thread.start()
            
        except Exception as e:
            messagebox.showerror("Error", f"Error starting optimization: {str(e)}")
    
    def _run_optimization_thread(self, weight, radius, max_power, max_current, max_voltage, wire_gauge):
        """Run optimization in background thread"""
        try:
            # Run optimization
            result = optimize_helmholtz(
                radius_m=radius,
                weight_limit_lbs=weight,
                max_power_w=max_power,
                fixed_current=max_current,
                fixed_voltage=max_voltage,
                fixed_awg=wire_gauge
            )
            
            # Update GUI with results
            self.root.after(0, self._display_optimization_results, result)
            
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, lambda: messagebox.showerror("Optimization Error", error_msg))
    
    def _display_optimization_results(self, result):
        """Display optimization results"""
        self.results_text.delete(1.0, tk.END)
        
        if result['success']:
            self.first_order_results = result
            
            # Display results
            self.results_text.insert(tk.END, "OPTIMIZATION SUCCESSFUL\n")
            self.results_text.insert(tk.END, "=" * 50 + "\n\n")
            
            # Wire specs
            self.results_text.insert(tk.END, "WIRE SPECIFICATION:\n")
            self.results_text.insert(tk.END, f"  AWG Size: {result['optimal_awg']:.1f} (use AWG{round(result['optimal_awg'])})\n\n")
            
            # Operating point
            op = result['operating_point']
            self.results_text.insert(tk.END, "OPTIMAL OPERATING POINT:\n")
            self.results_text.insert(tk.END, f"  Turns per Coil: {op['turns_per_coil']:.0f}\n")
            self.results_text.insert(tk.END, f"  Operating Current: {op['current_A']:.1f} A\n")
            self.results_text.insert(tk.END, f"  Operating Voltage: {op['voltage_V']:.1f} V\n")
            self.results_text.insert(tk.END, f"  Operating Power: {op['power_W']:.0f} W\n")
            self.results_text.insert(tk.END, f"  Single Coil Weight: {op['weight_single_coil_lbs']:.1f} lbs\n")
            self.results_text.insert(tk.END, f"  Both Coils Weight: {op['weight_both_coils_lbs']:.1f} lbs\n")
            self.results_text.insert(tk.END, f"  Wire Length: {op['wire_length_m']:.1f} m\n\n")
            
            # Performance
            self.results_text.insert(tk.END, "PERFORMANCE:\n")
            self.results_text.insert(tk.END, f"  Maximum B-field: {result['max_b_field_tesla']:.4f} T\n")

            self.results_text.insert(tk.END, "=" * 50 + "\n\n")
            self.results_text.insert(tk.END, "Use these results to to source components.\n")
            self.results_text.insert(tk.END, "• Wire gauges close to recommended value will perform similarly.\n")
            self.results_text.insert(tk.END, "• Power supplies should be programmable to the Operating Current, Voltage, and Power.\n")
            self.results_text.insert(tk.END, "• Geometric Optimization will provide a more accurate estimate of performance.\n")
            self.results_text.insert(tk.END, "=" * 50 + "\n\n")
            
            # Results displayed successfully
            
        else:
            self.results_text.insert(tk.END, "OPTIMIZATION FAILED\n")
            self.results_text.insert(tk.END, "=" * 50 + "\n\n")
            self.results_text.insert(tk.END, f"Error: {result.get('message', 'Unknown error')}\n")
    
    def validate_first_order_inputs(self):
        """Validate first order optimization inputs"""
        try:
            weight = self.weight_var.get()
            if weight <= 0:
                messagebox.showerror("Input Error", "Weight limit must be positive")
                return False
                
            radius = self.radius_var.get()
            if radius <= 0:
                messagebox.showerror("Input Error", "Radius must be positive")
                return False
                
            max_power = self.max_power_var.get()
            if max_power <= 0:
                messagebox.showerror("Input Error", "Max power must be positive")
                return False
                
            # Validate optional parameters if provided
            if self.max_current_var.get().strip():
                current = float(self.max_current_var.get())
                if current <= 0:
                    messagebox.showerror("Input Error", "Max current must be positive")
                    return False
                    
            if self.max_voltage_var.get().strip():
                voltage = float(self.max_voltage_var.get())
                if voltage <= 0:
                    messagebox.showerror("Input Error", "Max voltage must be positive")
                    return False
                    
            if self.wire_gauge_var.get().strip():
                awg = float(self.wire_gauge_var.get())
                if awg < 6 or awg > 56:
                    messagebox.showerror("Input Error", "Wire gauge must be between 6 and 56")
                    return False
                    
            return True
            
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numeric values")
            return False
    
    def clear_first_order_results(self):
        """Clear first order optimization results"""
        self.results_text.delete(1.0, tk.END)
        self.first_order_results = None
        # Results cleared
    
    def create_geometric_optimization_tab(self):
        """Create Geometric Optimization tab"""
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text="2. Geometric Optimization")
        
        # Main content frame
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Left panel - Input parameters
        left_panel = ttk.LabelFrame(content_frame, text="Geometric Optimization Parameters", padding="10")
        left_panel.pack(side="left", fill="y", padx=(0, 10))
        
        # All parameters in one section
        params_frame = ttk.Frame(left_panel)
        params_frame.pack(fill="x", pady=(0, 10))
        
        # Weight limit
        ttk.Label(params_frame, text="Weight Limit (lbs):").grid(row=0, column=0, sticky="w", pady=5)
        geo_weight_entry = ttk.Entry(params_frame, textvariable=self.geo_weight_var, width=15)
        geo_weight_entry.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Radius
        ttk.Label(params_frame, text="Radius (m):").grid(row=1, column=0, sticky="w", pady=5)
        geo_radius_entry = ttk.Entry(params_frame, textvariable=self.geo_radius_var, width=15)
        geo_radius_entry.grid(row=1, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Max Power
        ttk.Label(params_frame, text="Max Power (W):").grid(row=2, column=0, sticky="w", pady=5)
        geo_max_power_entry = ttk.Entry(params_frame, textvariable=self.geo_max_power_var, width=15)
        geo_max_power_entry.grid(row=2, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Max Current
        ttk.Label(params_frame, text="Max Current (A):").grid(row=3, column=0, sticky="w", pady=5)
        geo_current_entry = ttk.Entry(params_frame, textvariable=self.geo_max_current_var, width=15)
        geo_current_entry.grid(row=3, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Max Voltage
        ttk.Label(params_frame, text="Max Voltage (V):").grid(row=4, column=0, sticky="w", pady=5)
        geo_voltage_entry = ttk.Entry(params_frame, textvariable=self.geo_max_voltage_var, width=15)
        geo_voltage_entry.grid(row=4, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Wire Gauge
        ttk.Label(params_frame, text="Wire Gauge (AWG):").grid(row=5, column=0, sticky="w", pady=5)
        geo_wire_gauge_entry = ttk.Entry(params_frame, textvariable=self.geo_wire_gauge_var, width=15)
        geo_wire_gauge_entry.grid(row=5, column=1, sticky="w", padx=(10, 0), pady=5)
        
        
        # Control buttons
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill="x", pady=(10, 0))
        
        self.geo_optimize_button = ttk.Button(button_frame, text="Run Geometric Optimization", 
                                            command=self.run_geometric_optimization,
                                            style='Primary.TButton')
        self.geo_optimize_button.pack(side="left", padx=(0, 10))
        
        self.clear_geo_button = ttk.Button(button_frame, text="Clear Results", 
                                         command=self.clear_geometric_results,
                                         style='Action.TButton')
        self.clear_geo_button.pack(side="left", padx=(0, 10))
        
        self.sync_to_first_button = ttk.Button(button_frame, text="Sync to First Order", 
                                             command=self.sync_geometric_to_first_order,
                                             style='Action.TButton')
        self.sync_to_first_button.pack(side="left")
        
        # Right panel - Results and directions
        right_panel = ttk.Frame(content_frame)
        right_panel.pack(side="right", fill="both", expand=True)
        
        # Directions
        directions_frame = ttk.LabelFrame(right_panel, text="Directions", padding="10")
        directions_frame.pack(fill="x", pady=(0, 10))
        
        directions_text = """Geometric Optimization Stage:

1. Enter all required parameters:
   • Weight Limit: Maximum weight for a SINGLE COIL in pounds
   • Radius: Coil radius in meters  
   • Max Power: Maximum power budget in watts
   • Max Current: Maximum current capability
   • Max Voltage: Maximum voltage capability
   • Wire Gauge: AWG wire size

2. Max turns will be calculated automatically based on weight limit

3. Click 'Run Geometric Optimization' to find optimal coil geometry

4. Results will show optimal layer configuration and performance metrics"""
        
        directions_label = ttk.Label(directions_frame, text=directions_text, style='Info.TLabel', 
                                   justify="left", wraplength=400)
        directions_label.pack(anchor="w")
        
        # Results area
        results_frame = ttk.LabelFrame(right_panel, text="Geometric Optimization Results", padding="10")
        results_frame.pack(fill="both", expand=True)
        
        # Results text area
        self.geo_results_text = tk.Text(results_frame, height=20, width=60, wrap="word", 
                                       font=("Courier", 9))
        geo_scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.geo_results_text.yview)
        self.geo_results_text.configure(yscrollcommand=geo_scrollbar.set)
        
        self.geo_results_text.pack(side="left", fill="both", expand=True)
        geo_scrollbar.pack(side="right", fill="y")
    
    def run_geometric_optimization(self):
        """Run geometric optimization"""
        try:
            # Validate inputs
            if not self.validate_geometric_inputs():
                return
                
            # Get parameters from input fields
            weight_limit = self.geo_weight_var.get()
            radius = self.geo_radius_var.get()
            max_power = self.geo_max_power_var.get()
            max_current = self.geo_max_current_var.get()
            max_voltage = self.geo_max_voltage_var.get()
            wire_awg = self.geo_wire_gauge_var.get()
            field_z = 0.0  # Field point is always [0,0,0]
            
            # Create power supply specs
            ps_specs = {
                'max_power': max_power,
                'max_current': max_current,
                'max_voltage': max_voltage
            }
            
            # Update status
            self.geo_results_text.delete(1.0, tk.END)
            self.geo_results_text.insert(tk.END, "Running geometric optimization...\n")
            self.geo_results_text.insert(tk.END, f"Power Supply: {ps_specs['max_power']:.0f}W/{ps_specs['max_current']:.0f}A/{ps_specs['max_voltage']:.0f}V\n")
            self.geo_results_text.insert(tk.END, f"Wire: AWG {wire_awg:.0f}\n")
            self.geo_results_text.insert(tk.END, f"Radius: {radius} m\n")
            self.geo_results_text.insert(tk.END, f"Weight limit (single coil): {weight_limit} lbs\n")
            self.geo_results_text.insert(tk.END, f"Field point: [0, 0, 0]\n\n")
            self.geo_results_text.insert(tk.END, "Optimizing...\n")
            self.root.update()
            
            # Run optimization in separate thread
            thread = threading.Thread(target=self._run_geometric_optimization_thread, 
                                    args=(ps_specs, int(wire_awg), radius, weight_limit, 0.0))
            thread.daemon = True
            thread.start()
            
        except Exception as e:
            messagebox.showerror("Error", f"Error starting geometric optimization: {str(e)}")
    
    def validate_geometric_inputs(self):
        """Validate geometric optimization inputs"""
        try:
            weight = self.geo_weight_var.get()
            if weight <= 0:
                messagebox.showerror("Input Error", "Weight limit must be positive")
                return False
                
            radius = self.geo_radius_var.get()
            if radius <= 0:
                messagebox.showerror("Input Error", "Radius must be positive")
                return False
                
            max_power = self.geo_max_power_var.get()
            if max_power <= 0:
                messagebox.showerror("Input Error", "Max power must be positive")
                return False
                
            max_current = self.geo_max_current_var.get()
            if max_current <= 0:
                messagebox.showerror("Input Error", "Max current must be positive")
                return False
                
            max_voltage = self.geo_max_voltage_var.get()
            if max_voltage <= 0:
                messagebox.showerror("Input Error", "Max voltage must be positive")
                return False
                
            wire_awg = self.geo_wire_gauge_var.get()
            if wire_awg < 6 or wire_awg > 56:
                messagebox.showerror("Input Error", "Wire gauge must be between 6 and 56")
                return False
                
            return True
            
        except ValueError:
            messagebox.showerror("Input Error", "Please enter valid numeric values")
            return False
    
    def sync_first_order_to_geometric(self):
        """Sync values from First Order Optimization to Geometric Optimization"""
        try:
            # Sync basic parameters
            self.geo_weight_var.set(self.weight_var.get())
            self.geo_radius_var.set(self.radius_var.get())
            self.geo_max_power_var.set(self.max_power_var.get())
            
            # Sync optional parameters if they have values
            if self.max_current_var.get().strip():
                self.geo_max_current_var.set(float(self.max_current_var.get()))
            if self.max_voltage_var.get().strip():
                self.geo_max_voltage_var.set(float(self.max_voltage_var.get()))
            if self.wire_gauge_var.get().strip():
                self.geo_wire_gauge_var.set(float(self.wire_gauge_var.get()))
                
        except Exception as e:
            messagebox.showerror("Sync Error", f"Error syncing values: {str(e)}")
    
    def sync_geometric_to_first_order(self):
        """Sync values from Geometric Optimization to First Order Optimization"""
        try:
            # Sync basic parameters
            self.weight_var.set(self.geo_weight_var.get())
            self.radius_var.set(self.geo_radius_var.get())
            self.max_power_var.set(self.geo_max_power_var.get())
            
            # Sync geometric-specific parameters to optional fields
            self.max_current_var.set(str(self.geo_max_current_var.get()))
            self.max_voltage_var.set(str(self.geo_max_voltage_var.get()))
            self.wire_gauge_var.set(str(self.geo_wire_gauge_var.get()))
            
        except Exception as e:
            messagebox.showerror("Sync Error", f"Error syncing values: {str(e)}")
    
    def _run_geometric_optimization_thread(self, ps_specs, wire_awg, radius, weight_limit, field_z):
        """Run geometric optimization in background thread"""
        try:
            # Run optimization
            result = optimize_electromagnet_geometry(
                power_supply_specs=ps_specs,
                wire_awg=wire_awg,
                base_radius=radius,
                weight_limit=weight_limit,
                field_point=[0, 0, field_z],
                max_iterations=100,
                verbose=True
            )
            
            # Update GUI with results
            self.root.after(0, self._display_geometric_results, result, radius, wire_awg)
            
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, lambda: messagebox.showerror("Optimization Error", error_msg))
    
    def _display_geometric_results(self, result, radius, wire_awg):
        """Display geometric optimization results"""
        self.geo_results_text.delete(1.0, tk.END)
        
        if result['success']:
            self.geometric_results = result
            
            # Display results
            self.geo_results_text.insert(tk.END, "GEOMETRIC OPTIMIZATION SUCCESSFUL\n")
            self.geo_results_text.insert(tk.END, "=" * 50 + "\n\n")
            
            # Optimal geometry
            geom = result['optimal_geometry']
            self.geo_results_text.insert(tk.END, "OPTIMAL GEOMETRY:\n")
            self.geo_results_text.insert(tk.END, f"  Radial layers: {geom['radial_layers']}\n")
            
            # Calculate geometry dimensions
            from electromagnet import Wire
            wire = Wire(int(wire_awg))
            radial_depth = (geom['radial_layers'] - 1) * wire.diameter_nom_m
            outside_radius = radius + radial_depth
            axial_height = geom['axial_layers'] * wire.diameter_nom_m
            
            self.geo_results_text.insert(tk.END, f"  Radial depth: {radial_depth*1000:.1f} mm\n")
            self.geo_results_text.insert(tk.END, f"  Outside radius: {outside_radius*1000:.1f} mm\n")
            self.geo_results_text.insert(tk.END, f"  Axial layers: {geom['axial_layers']}\n")
            self.geo_results_text.insert(tk.END, f"  Axial height: {axial_height*1000:.1f} mm\n")
            self.geo_results_text.insert(tk.END, f"  Total single coil turns: {geom['radial_layers'] * geom['axial_layers']}\n\n")
            
            # Performance
            perf = result['performance']
            self.geo_results_text.insert(tk.END, "PERFORMANCE:\n")
            self.geo_results_text.insert(tk.END, f"  B-field: {perf['b_field_magnitude']:.4f} T\n")
            self.geo_results_text.insert(tk.END, f"  Weight (both coils): {perf['weight']:.1f} lbs\n")
            self.geo_results_text.insert(tk.END, f"  Weight (single coil): {perf['weight']/2:.1f} lbs\n")
            self.geo_results_text.insert(tk.END, f"  Power: {perf['power']:.0f} W\n")
            self.geo_results_text.insert(tk.END, f"  Current: {perf['current']:.1f} A\n")
            self.geo_results_text.insert(tk.END, f"  Voltage: {perf['voltage']:.1f} V\n")
            self.geo_results_text.insert(tk.END, f"  Wire length: {perf['wire_length']:.1f} m\n")
            self.geo_results_text.insert(tk.END, f"  Total turns: {perf['turns_total']}\n")
            self.geo_results_text.insert(tk.END, f"  Single coil turns: {perf['helm_turns']}\n\n")
            
            # Optimization info
            info = result['optimization_info']
            self.geo_results_text.insert(tk.END, "OPTIMIZATION INFO:\n")
            self.geo_results_text.insert(tk.END, f"  Global iterations: {info['global_iterations']}\n")
            self.geo_results_text.insert(tk.END, f"  Local iterations: {info['local_iterations']}\n")
            self.geo_results_text.insert(tk.END, f"  Total evaluations: {info['total_evaluations']}\n")
            
            # Results displayed successfully
            
        else:
            self.geo_results_text.insert(tk.END, "GEOMETRIC OPTIMIZATION FAILED\n")
            self.geo_results_text.insert(tk.END, "=" * 50 + "\n\n")
            self.geo_results_text.insert(tk.END, f"Error: {result.get('message', 'Unknown error')}\n")
    
    def clear_geometric_results(self):
        """Clear geometric optimization results"""
        self.geo_results_text.delete(1.0, tk.END)
        self.geometric_results = None
        # Results cleared
    

def main():
    """Main function to run the GUI"""
    root = tk.Tk()
    app = HelmholtzGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
