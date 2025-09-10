'''
Author: Aren Mowreader
Date: 7/2/2025
Version: 1.0
Description:

Script to create a GUI for the Helmholtz simulation.
    Two different modes:
        -First Order Approximation. Runs Faster but less accurate.
        -Line Integral. Runs slower but more accurate.
    The GUI will have the following features:
        Inputs:
            First Order Approximation:
                -Radius of the spiral
                -Wire Gauge
                -Number of turns
                -Power Supply
                    -Max Power
                    -Max Current
                    -Max Voltage
                -Weight Limit
            Line Integral:
                -Radius of the spiral
                -Wire Gauge
                -Turns per layer
                -Layers
                -Power Supply
                    -Max Power
                    -Max Current
                    -Max Voltage
                -Weight Limit
        Outputs:
            First Order Approximation:
                -B-Field Curve
                -Power Supply Response Curves
                    -Current vs Turns
                    -Voltage vs Turns
                    -Power vs Turns
                Chart of critical values at weight limit
            Line Integral:
                -B-Field Curve
                -Power Supply Response Curves
                    -Current vs Turns
                    -Voltage vs Turns
                    -Power vs Turns
                Chart of critical values at weight limit
                xz cross section of B-Field Magnitude at weight limit
                xz cross section of B-Field vector at weight limit
                Uniform volume projection


    '''

import tkinter as tk
from tkinter import ttk
from tkinter import messagebox
from tkinter import filedialog
from tkinter import scrolledtext
from tkinter import simpledialog
from tkinter import colorchooser
from tkinter import messagebox
import matplotlib.pyplot as plt
from matplotlib.backends.backend_tkagg import FigureCanvasTkAgg
from matplotlib.figure import Figure
import numpy as np
import pandas as pd
import helmholtz_utils as helm
import threading
import time
import os


class HelmholtzGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Helmholtz Coil Simulation GUI")
        self.root.geometry("1200x800")
        
        # Variables
        self.results_data = {}
        self.is_simulating = False
        self.first_order_completed = False
        self.last_first_order_params = None
        
        # Create main containers
        self.create_widgets()
        self.setup_layout()
        
        # Set up parameter change tracking
        self.setup_parameter_tracking()
        
    def create_widgets(self):
        """Create all GUI widgets"""
        # Mode selection
        self.create_mode_selection()
        
        # Input frames
        self.create_input_frames()
        
        # Control buttons
        self.create_control_buttons()
        
        # Results area
        self.create_results_area()
        
        # Status bar
        self.create_status_bar()
        
    def create_mode_selection(self):
        """Create mode selection frame"""
        mode_frame = ttk.LabelFrame(self.root, text="Simulation Workflow", padding="10")
        mode_frame.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        
        # Instructions
        ttk.Label(mode_frame, text="1. Run First Order Approximation to explore parameters", 
                 font=("Arial", 10, "bold")).pack(anchor="w", pady=2)
        ttk.Label(mode_frame, text="2. Run Field Analysis for detailed 2D field visualization", 
                 font=("Arial", 10, "bold")).pack(anchor="w", pady=2)
        
    def create_input_frames(self):
        """Create input parameter frames"""
        # Coil parameters frame
        self.coil_frame = ttk.LabelFrame(self.root, text="Coil Parameters", padding="10")
        self.coil_frame.grid(row=1, column=0, sticky="nsew", padx=10, pady=5)
        
        # Power supply frame
        self.power_frame = ttk.LabelFrame(self.root, text="Power Supply", padding="10")
        self.power_frame.grid(row=2, column=0, sticky="nsew", padx=10, pady=5)
        
        # Constraints frame
        self.constraints_frame = ttk.LabelFrame(self.root, text="Constraints", padding="10")
        self.constraints_frame.grid(row=3, column=0, sticky="nsew", padx=10, pady=5)
        
        # Create input widgets
        self.create_coil_inputs()
        self.create_power_inputs()
        self.create_constraint_inputs()
        
    def create_coil_inputs(self):
        """Create coil parameter input widgets"""
        # Radius
        ttk.Label(self.coil_frame, text="Radius (m):").grid(row=0, column=0, sticky="w", pady=2)
        self.radius_var = tk.DoubleVar(value=0.253)
        self.radius_entry = ttk.Entry(self.coil_frame, textvariable=self.radius_var, width=15)
        self.radius_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        
        # Wire gauge
        ttk.Label(self.coil_frame, text="Wire Gauge (AWG):").grid(row=1, column=0, sticky="w", pady=2)
        self.wire_gauge_var = tk.IntVar(value=12)
        self.wire_gauge_entry = ttk.Entry(self.coil_frame, textvariable=self.wire_gauge_var, width=15)
        self.wire_gauge_entry.grid(row=1, column=1, sticky="w", padx=5, pady=2)
        
        # Number of turns (First Order)
        ttk.Label(self.coil_frame, text="Number of Turns:").grid(row=2, column=0, sticky="w", pady=2)
        self.turns_var = tk.IntVar(value=2000)
        self.turns_entry = ttk.Entry(self.coil_frame, textvariable=self.turns_var, width=15)
        self.turns_entry.grid(row=2, column=1, sticky="w", padx=5, pady=2)
        
        # Coil Spacing (Line Integral)
        ttk.Label(self.coil_frame, text="Coil Spacing (m):").grid(row=3, column=0, sticky="w", pady=2)
        self.coil_spacing_var = tk.DoubleVar(value=0.253)  # Default to radius value
        self.coil_spacing_entry = ttk.Entry(self.coil_frame, textvariable=self.coil_spacing_var, width=15)
        self.coil_spacing_entry.grid(row=3, column=1, sticky="w", padx=5, pady=2)
        
    def create_power_inputs(self):
        """Create power supply input widgets"""
        # Power Supply Name
        ttk.Label(self.power_frame, text="Power Supply Name: (optional)").grid(row=0, column=0, sticky="w", pady=2)
        self.power_supply_name_var = tk.StringVar(value="")
        self.power_supply_name_entry = ttk.Entry(self.power_frame, textvariable=self.power_supply_name_var, width=15)
        self.power_supply_name_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        
        # Max Power
        ttk.Label(self.power_frame, text="Max Power (W):").grid(row=1, column=0, sticky="w", pady=2)
        self.max_power_var = tk.DoubleVar(value=10000.0)
        self.max_power_entry = ttk.Entry(self.power_frame, textvariable=self.max_power_var, width=15)
        self.max_power_entry.grid(row=1, column=1, sticky="w", padx=5, pady=2)
        
        # Max Current
        ttk.Label(self.power_frame, text="Max Current (A):").grid(row=2, column=0, sticky="w", pady=2)
        self.max_current_var = tk.DoubleVar(value=40.0)
        self.max_current_entry = ttk.Entry(self.power_frame, textvariable=self.max_current_var, width=15)
        self.max_current_entry.grid(row=2, column=1, sticky="w", padx=5, pady=2)
        
        # Max Voltage
        ttk.Label(self.power_frame, text="Max Voltage (V):").grid(row=3, column=0, sticky="w", pady=2)
        self.max_voltage_var = tk.DoubleVar(value=500.0)
        self.max_voltage_entry = ttk.Entry(self.power_frame, textvariable=self.max_voltage_var, width=15)
        self.max_voltage_entry.grid(row=3, column=1, sticky="w", padx=5, pady=2)
        
    def create_constraint_inputs(self):
        """Create constraint input widgets"""
        # Weight limit
        ttk.Label(self.constraints_frame, text="Weight Limit (lbs):").grid(row=0, column=0, sticky="w", pady=2)
        self.weight_limit_var = tk.DoubleVar(value=75.0)
        self.weight_limit_entry = ttk.Entry(self.constraints_frame, textvariable=self.weight_limit_var, width=15)
        self.weight_limit_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        
    def create_control_buttons(self):
        """Create control button frame"""
        button_frame = ttk.Frame(self.root)
        button_frame.grid(row=4, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        
        # First Order button
        self.first_order_button = ttk.Button(button_frame, text="Run First Order Simulation", 
                                            command=self.run_first_order)
        self.first_order_button.pack(side="left", padx=5)
        
        # Field Analysis button (initially disabled)
        self.field_analysis_button = ttk.Button(button_frame, text="Run Field Analysis", 
                                               command=self.run_field_analysis, state="disabled")
        self.field_analysis_button.pack(side="left", padx=5)
        
        # Stop button
        self.stop_button = ttk.Button(button_frame, text="Stop", command=self.stop_simulation, state="disabled")
        self.stop_button.pack(side="left", padx=5)
        
        # Save results button
        self.save_button = ttk.Button(button_frame, text="Save Results", command=self.save_results)
        self.save_button.pack(side="left", padx=5)
        
        # Save plots button
        self.save_plots_button = ttk.Button(button_frame, text="Save Plots", command=self.save_plots)
        self.save_plots_button.pack(side="left", padx=5)
        
        # Clear button
        self.clear_button = ttk.Button(button_frame, text="Clear Results", command=self.clear_results)
        self.clear_button.pack(side="left", padx=5)
        
    def create_results_area(self):
        """Create results display area"""
        results_frame = ttk.LabelFrame(self.root, text="Results", padding="10")
        results_frame.grid(row=1, column=1, rowspan=3, sticky="nsew", padx=10, pady=5)
        
        # Notebook for tabbed results
        self.notebook = ttk.Notebook(results_frame)
        self.notebook.pack(fill="both", expand=True)
        
        # First Order Plots tab
        self.first_order_plots_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.first_order_plots_frame, text="First Order Plots")
        
        # Field Analysis tab (Line Integral plots)
        self.field_analysis_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.field_analysis_frame, text="Field Analysis")
        
        # Data tab
        self.data_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.data_frame, text="Data")
        
        # Create matplotlib figures for each tab
        self.fig_first_order = Figure(figsize=(8, 6))
        self.canvas_first_order = FigureCanvasTkAgg(self.fig_first_order, self.first_order_plots_frame)
        self.canvas_first_order.get_tk_widget().pack(fill="both", expand=True)
        
        self.fig_field_analysis = Figure(figsize=(8, 6))
        self.canvas_field_analysis = FigureCanvasTkAgg(self.fig_field_analysis, self.field_analysis_frame)
        self.canvas_field_analysis.get_tk_widget().pack(fill="both", expand=True)
        
        # Data text area
        self.data_text = scrolledtext.ScrolledText(self.data_frame, height=20, width=60)
        self.data_text.pack(fill="both", expand=True)
        
    def create_status_bar(self):
        """Create status bar"""
        self.status_var = tk.StringVar(value="Ready")
        self.status_bar = ttk.Label(self.root, textvariable=self.status_var, relief="sunken")
        self.status_bar.grid(row=5, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        
        # Progress bar
        self.progress = ttk.Progressbar(self.root, mode='indeterminate')
        self.progress.grid(row=6, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        
    def setup_layout(self):
        """Configure grid weights for responsive layout"""
        self.root.grid_columnconfigure(1, weight=1)
        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_rowconfigure(2, weight=1)
        self.root.grid_rowconfigure(3, weight=1)
        
    def setup_parameter_tracking(self):
        """Set up tracking for parameter changes"""
        # Enable coil spacing by default (for Line Integral mode)
        self.coil_spacing_entry.config(state="normal")
        
    def check_parameter_changes(self):
        """Check if parameters have changed since last First Order run"""
        if not self.first_order_completed:
            return True
            
        current_params = self.get_simulation_parameters()
        if self.last_first_order_params != current_params:
            self.first_order_completed = False
            self.field_analysis_button.config(state="disabled")
            return True
        return False
        
    def run_first_order(self):
        """Run First Order Approximation simulation"""
        if not self.validate_inputs():
            return
            
        if self.is_simulating:
            messagebox.showwarning("Simulation Running", "A simulation is already running!")
            return
            
        self.is_simulating = True
        self.first_order_button.config(state="disabled")
        self.field_analysis_button.config(state="disabled")
        self.stop_button.config(state="normal")
        self.status_var.set("Running First Order simulation...")
        self.progress.start()
        
        # Start simulation in separate thread
        thread = threading.Thread(target=self._run_first_order_thread)
        thread.daemon = True
        thread.start()
        
    def run_field_analysis(self):
        """Run Field Analysis (Line Integral) simulation"""
        if not self.first_order_completed:
            messagebox.showwarning("First Order Required", "Please run First Order simulation first!")
            return
            
        if not self.validate_inputs():
            return
            
        if self.is_simulating:
            messagebox.showwarning("Simulation Running", "A simulation is already running!")
            return
            
        self.is_simulating = True
        self.first_order_button.config(state="disabled")
        self.field_analysis_button.config(state="disabled")
        self.stop_button.config(state="normal")
        self.status_var.set("Running Field Analysis...")
        self.progress.start()
        
        # Start simulation in separate thread
        thread = threading.Thread(target=self._run_field_analysis_thread)
        thread.daemon = True
        thread.start()
            
    def validate_inputs(self):
        """Validate all input parameters"""
        try:
            # Validate numeric inputs
            radius = self.radius_var.get()
            if radius <= 0:
                raise ValueError("Radius must be positive")
                
            wire_gauge = self.wire_gauge_var.get()
            if wire_gauge < 6 or wire_gauge > 56:
                raise ValueError("Wire gauge must be between 6 and 56")
                
            turns = self.turns_var.get()
            if turns <= 0:
                raise ValueError("Number of turns must be positive")
                
            coil_spacing = self.coil_spacing_var.get()
            if coil_spacing <= 0:
                raise ValueError("Coil spacing must be positive")
                
            max_power = self.max_power_var.get()
            if max_power <= 0:
                raise ValueError("Max power must be positive")
                
            max_current = self.max_current_var.get()
            if max_current <= 0:
                raise ValueError("Max current must be positive")
                
            max_voltage = self.max_voltage_var.get()
            if max_voltage <= 0:
                raise ValueError("Max voltage must be positive")
                
            weight_limit = self.weight_limit_var.get()
            if weight_limit <= 0:
                raise ValueError("Weight limit must be positive")
                
            return True
            
        except ValueError as e:
            messagebox.showerror("Input Error", str(e))
            return False
            
    def run_simulation(self):
        """Run the simulation in a separate thread"""
        if not self.validate_inputs():
            return
            
        if self.is_simulating:
            messagebox.showwarning("Simulation Running", "A simulation is already running!")
            return
            
        self.is_simulating = True
        self.run_button.config(state="disabled")
        self.stop_button.config(state="normal")
        self.status_var.set("Running simulation...")
        self.progress.start()
        
        # Start simulation in separate thread
        thread = threading.Thread(target=self._run_simulation_thread)
        thread.daemon = True
        thread.start()
        
    def _run_first_order_thread(self):
        """Run First Order simulation in background thread"""
        try:
            # Get parameters
            params = self.get_simulation_parameters()
            
            # Run First Order simulation
            results = self.run_first_order_simulation(params)
            
            # Mark as completed and store parameters
            self.first_order_completed = True
            self.last_first_order_params = params.copy()
            
            # Update GUI with results
            self.root.after(0, self.display_first_order_results, results)
            
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, lambda: messagebox.showerror("Simulation Error", error_msg))
        finally:
            self.root.after(0, self.simulation_finished)
            
    def _run_field_analysis_thread(self):
        """Run Field Analysis simulation in background thread"""
        try:
            # Get parameters
            params = self.get_simulation_parameters()
            
            # Run Line Integral simulation
            results = self.run_line_integral_simulation(params)
            
            # Update GUI with results
            self.root.after(0, self.display_field_analysis_results, results)
            
        except Exception as e:
            error_msg = str(e)
            self.root.after(0, lambda: messagebox.showerror("Simulation Error", error_msg))
        finally:
            self.root.after(0, self.simulation_finished)
            
    def get_simulation_parameters(self):
        """Get all simulation parameters"""
        params = {
            'radius': self.radius_var.get(),
            'wire_gauge': self.wire_gauge_var.get(),
            'turns': self.turns_var.get(),
            'coil_spacing': self.coil_spacing_var.get(),
            'power_supply_name': self.power_supply_name_var.get(),
            'max_power': self.max_power_var.get(),
            'max_current': self.max_current_var.get(),
            'max_voltage': self.max_voltage_var.get(),
            'weight_limit': self.weight_limit_var.get()
        }
        
        # Calculate turns_per_layer and layers from total turns (for Line Integral)
        total_turns = params['turns']
        # Simple calculation: assume 2 layers with equal turns per layer
        params['layers'] = 2
        params['turns_per_layer'] = total_turns // 2
        # If odd number of turns, add the remainder to the first layer
        if total_turns % 2 != 0:
            params['turns_per_layer'] += 1
        
        return params
        
    def run_first_order_simulation(self, params):
        """Run first order approximation simulation"""
        # Get wire data
        script_dir = os.path.dirname(os.path.abspath(__file__))
        csv_path = os.path.join(script_dir, 'MWS_wire_data.csv')
        wire_data = pd.read_csv(csv_path)
        ohms_per_meter, lbs_per_meter, wire_diameter = helm.get_wire_data(
            wire_data, params['wire_gauge'], params['wire_gauge'])
        
        # Create power supply DataFrame (matching pmax.py pattern)
        power_supplies = pd.DataFrame({
            'name': ['Custom'],
            'Vmax': [params['max_voltage']],
            'Imax': [params['max_current']],
            'Pmax': [params['max_power']]
        })
        
        # Select the power supply (matching pmax.py pattern)
        power_supply = power_supplies.iloc[0]
        
        
        # Create range of turns to analyze
        turns_range = np.arange(1, params['turns'] + 1) 
        
        # Run analysis
        results_df = helm.analyze_power_supply(power_supply, ohms_per_meter[0],
                                                turns_range, params['radius'])
        turns_at_weight = helm.weight_limit_n_one_coil(params['weight_limit'], 
                                           params['radius'], lbs_per_meter)
        
        # Calculate max B-field data
        max_b_field_index = np.argmax(results_df['B_power'].values)
        turns_at_max_b = turns_range[max_b_field_index] 
        weight_at_max_b = helm.weight(turns_at_max_b, params['radius'],
                                     lbs_per_meter[0])
        
        # Extract data for GUI plotting
        return {
            'mode': 'First Order Approximation',
            'turns': turns_range,
            'current': results_df['I'].values,
            'voltage': results_df['V'].values,
            'power': results_df['power'].values,
            'b_field': results_df['B_power'].values,
            'b_field_current': results_df['B_current'].values,  # Add B_current data
            'turns_at_weight': turns_at_weight,
            'turns_at_max_b': turns_at_max_b,
            'weight_at_max_b': weight_at_max_b,
            'params': params,
            'lbs_per_meter': lbs_per_meter[0]  # Add wire weight data for calculations
        }
        
    def run_line_integral_simulation(self, params):
        """Run line integral simulation"""
        # Get the stored First Order results to use turns_at_weight
        if not hasattr(self, 'results_data') or 'turns_at_weight' not in self.results_data:
            raise ValueError("First Order simulation must be run before Field Analysis")
            
        turns_at_weight = self.results_data['turns_at_weight']
        
        wire_data = pd.read_csv(os.path.join(os.path.dirname(__file__), 
                                             'MWS_wire_data.csv'))
        ohms_per_meter, lbs_per_meter, wire_diameter = helm.get_wire_data(wire_data,
                                     params['wire_gauge'], params['wire_gauge'])
        
        
        wire_length = helm.wire_length_square(params['radius'], 
                    turns_at_weight, wire_diameter)
        
        power = helm.power_at_length(params['max_power'], params['max_current'], 
                                 params['max_voltage'], wire_length,
                                 ohms_per_meter[0])
        # Convert to scalar if it's an array
        if hasattr(power, '__len__'):
            power = float(power[0])
        else:
            power = float(power)
            
        I = np.sqrt(power/(ohms_per_meter[0]*wire_length))
        # Convert to scalar if it's an array
        if hasattr(I, '__len__'):
            I = float(I[0])
        else:
            I = float(I)
            
        height = params['coil_spacing'] + .2 
        width = params['radius']*2 + .2
        X, Z, grid_positions = helm.grid_positions_xz(width, height, .05)
        b_vectors = helm.Helm(params['radius'], wire_diameter, params['turns_per_layer'],
                               params['layers'], 1000, params['coil_spacing']/2, grid_positions, I)
        print("helm good")
        b_magnitude = np.linalg.norm(b_vectors, axis=1)
        #reshape b_magnitude and b_vectors to 2d arrays using X shape
        b_magnitude = b_magnitude.reshape(X.shape)
        b_vectors = b_vectors.reshape(X.shape[0], X.shape[1], 3)
        # Extract center value for single point analysis
        center_z = X.shape[0] // 2
        center_x = X.shape[1] // 2
        b_field_center = b_magnitude[center_z, center_x]

        print("grid good")
        return {
            'mode': 'Line Integral',
            'params': params,
            # Single point analysis data
            'b_field_center': float(b_field_center),
            'li_current': I,
            'li_voltage': power / I,
            'li_power': power,
            'turns_total': params['turns_per_layer'] * params['layers'],
            'wire_length': float(wire_length),
            'coil_weight': float(wire_length * lbs_per_meter[0] / 2),  # /2 for one coil
            # 2D field data for visualization
            'b_field_2d': b_magnitude,
            'b_vectors_2d': b_vectors,
            'X_grid': X,
            'Z_grid': Z
        }
        
    def display_first_order_results(self, results):
        """Display First Order simulation results"""
        results['mode'] = 'First Order Approximation'
        self.results_data = results
        
        # Update plots
        self.update_plots(results)
        
        # Update data text
        self.update_data_text(results)
        
        # Enable Field Analysis button
        self.field_analysis_button.config(state="normal")
        
        # Switch to results tab
        self.notebook.select(0)
        
    def display_field_analysis_results(self, results):
        """Display Field Analysis simulation results"""
        results['mode'] = 'Line Integral'
        self.results_data = results
        
        # Update plots
        self.update_plots(results)
        
        # Update data text
        self.update_data_text(results)
        
        # Switch to results tab
        self.notebook.select(0)
        
    def update_plots(self, results):
        """Update matplotlib plots"""
        if results['mode'] == 'Line Integral':
            # Update Field Analysis tab with Line Integral plots
            self.fig_field_analysis.clear()
            self.update_line_integral_plots(results)
            
            # Create title with optional power supply name
            title = f'Field Analysis - {results["mode"]}'
            if results['params']['power_supply_name'].strip():
                title += f' - {results["params"]["power_supply_name"]}'
            self.fig_field_analysis.suptitle(title, fontsize=14)
            self.canvas_field_analysis.draw()
            
            # Also update First Order plots for comparison
            self.fig_first_order.clear()
            first_order_results = self.run_first_order_simulation(results['params'])
            self.update_first_order_plots(first_order_results)
            
            # Create title for First Order plots
            title = f'First Order Approximation Comparison'
            if results['params']['power_supply_name'].strip():
                title += f' - {results["params"]["power_supply_name"]}'
            self.fig_first_order.suptitle(title, fontsize=14)
            self.canvas_first_order.draw()
            
        else:
            # Update First Order plots tab
            self.fig_first_order.clear()
            self.update_first_order_plots(results)
            
            # Create title with optional power supply name
            title = f'{results["mode"]} Simulation Results'
            if results['params']['power_supply_name'].strip():
                title += f' - {results["params"]["power_supply_name"]}'
            self.fig_first_order.suptitle(title, fontsize=14)
            self.canvas_first_order.draw()
            
            # Clear Field Analysis tab
            self.fig_field_analysis.clear()
            self.canvas_field_analysis.draw()
        
    def update_first_order_plots(self, results):
        """Update plots for First Order Approximation mode"""
        # Create subplots
        gs = self.fig_first_order.add_gridspec(2, 2, hspace=0.5, wspace=0.3)
        
        # B-Field curve (both B_power and B_current)
        ax1 = self.fig_first_order.add_subplot(gs[0, 0])
        ax1.plot(results['turns'], results['b_field'], 'b-', linewidth=2, label='B-Field (Power-based)')
        ax1.plot(results['turns'], results['b_field_current'], 'r--', linewidth=2, label='B-Field (Current-based)')
        ax1.set_xlabel('Number of Turns')
        ax1.set_ylabel('B-Field (T)')
        ax1.set_title('B-Field vs Turns')
        ax1.grid(True)
        ax1.legend()
        
        # Add vertical line at weight limit (if available)
        if 'turns_at_weight' in results:
            ax1.axvline(x=results['turns_at_weight'], color='gray', linestyle=':', linewidth=2, label='Weight Limit')
            ax1.legend()
        
        # Current vs Turns
        ax2 = self.fig_first_order.add_subplot(gs[0, 1])
        ax2.plot(results['turns'], results['current'], 'r-', linewidth=2)
        ax2.set_xlabel('Number of Turns')
        ax2.set_ylabel('Current (A)')
        ax2.set_title('Current vs Turns')
        ax2.grid(True)
        
        # Add vertical line at weight limit (if available)
        if 'turns_at_weight' in results:
            ax2.axvline(x=results['turns_at_weight'], color='gray', linestyle=':', linewidth=2, label='Weight Limit')
            ax2.legend()
        
        # Voltage vs Turns
        ax3 = self.fig_first_order.add_subplot(gs[1, 0])
        ax3.plot(results['turns'], results['voltage'], 'g-', linewidth=2)
        ax3.set_xlabel('Number of Turns')
        ax3.set_ylabel('Voltage (V)')
        ax3.set_title('Voltage vs Turns')
        ax3.grid(True)
        
        # Add vertical line at weight limit (if available)
        if 'turns_at_weight' in results:
            ax3.axvline(x=results['turns_at_weight'], color='gray', linestyle=':', linewidth=2, label='Weight Limit')
            ax3.legend()
        
        # Power vs Turns
        ax4 = self.fig_first_order.add_subplot(gs[1, 1])
        ax4.plot(results['turns'], results['power'], 'm-', linewidth=2)
        ax4.set_xlabel('Number of Turns')
        ax4.set_ylabel('Power (W)')
        ax4.set_title('Power vs Turns')
        ax4.grid(True)
        
        # Add vertical line at weight limit (if available)
        if 'turns_at_weight' in results:
            ax4.axvline(x=results['turns_at_weight'], color='gray', linestyle=':', linewidth=2, label='Weight Limit')
            ax4.legend()
            
    def update_line_integral_plots(self, results):
        """Update plots for Line Integral mode with 2D field visualizations"""
        # Create subplots: 1x2 layout for just the field visualizations
        gs = self.fig_field_analysis.add_gridspec(1, 2, hspace=0.3, wspace=0.3)
        
        # Get 2D field data
        b_field_2d = results['b_field_2d']
        b_vectors_2d = results['b_vectors_2d']
        X = results['X_grid']
        Z = results['Z_grid']
        
        # Plot 1: B-Field Magnitude Contour (Topological Map)
        ax1 = self.fig_field_analysis.add_subplot(gs[0, 0])
        contour = ax1.contourf(X, Z, b_field_2d, levels=20, cmap='viridis')
        ax1.set_xlabel('X (m)')
        ax1.set_ylabel('Z (m)')
        ax1.set_title('B-Field Magnitude (T)')
        ax1.set_aspect('equal')
        ax1.grid(True, alpha=0.3)
        self.fig_field_analysis.colorbar(contour, ax=ax1, label='Field Strength (T)')
        
        # Plot 2: Vector Field
        ax2 = self.fig_field_analysis.add_subplot(gs[0, 1])
        # Sample every nth point to avoid overcrowding
        n = max(1, min(X.shape[0], X.shape[1]) // 20)
        quiver = ax2.quiver(X[::n, ::n], Z[::n, ::n], 
                           b_vectors_2d[::n, ::n, 0], b_vectors_2d[::n, ::n, 2],
                           scale=50, color='blue', alpha=0.7)
        ax2.set_xlabel('X (m)')
        ax2.set_ylabel('Z (m)')
        ax2.set_title('B-Field Vector Direction')
        ax2.set_aspect('equal')
        ax2.grid(True, alpha=0.3)
        
    def update_data_text(self, results):
        """Update data text area"""
        self.data_text.delete(1.0, tk.END)
        
        text = f"Simulation Mode: {results['mode']}\n"
        if results['params']['power_supply_name'].strip():
            text += f"Power Supply: {results['params']['power_supply_name']}\n"
        text += "=" * 50 + "\n\n"
        
        text += "Parameters:\n"
        for key, value in results['params'].items():
            text += f"  {key}: {value}\n"
        
        text += "\nResults Summary:\n"
        
        if results['mode'] == 'Line Integral':
            # First Order Approximation results (for comparison)
            text += "\nFIRST ORDER APPROXIMATION RESULTS:\n"
            text += "-" * 40 + "\n"
            
            # Run first order analysis with same parameters
            first_order_results = self.run_first_order_simulation(results['params'])
            
            text += f"  Max B-Field (Power-based): {np.max(first_order_results['b_field']):.6f} T\n"
            text += f"  Max B-Field (Current-based): {np.max(first_order_results['b_field_current']):.6f} T\n"
            if 'turns_at_max_b' in first_order_results:
                text += f"  Turns at Max B-Field: {first_order_results['turns_at_max_b']:.0f}\n"
                text += f"  Weight at Max B-Field: {first_order_results['weight_at_max_b']:.2f} lbs\n"
            text += f"  Turns at Weight Limit: {first_order_results['turns_at_weight']:.0f}\n"
            weight_limit_index = min(first_order_results['turns_at_weight']-1, len(first_order_results['b_field'])-1)
            if weight_limit_index >= 0:
                text += f"  B-Field (Power) at Weight Limit: {first_order_results['b_field'][weight_limit_index]:.6f} T\n"
                text += f"  B-Field (Current) at Weight Limit: {first_order_results['b_field_current'][weight_limit_index]:.6f} T\n"
                text += f"  Current at Weight Limit: {first_order_results['current'][weight_limit_index]:.2f} A\n"
                text += f"  Voltage at Weight Limit: {first_order_results['voltage'][weight_limit_index]:.2f} V\n"
                text += f"  Power at Weight Limit: {first_order_results['power'][weight_limit_index]:.2f} W\n"
            else:
                text += f"  Weight limit exceeds simulation range\n"
            
            text += "\nLINE INTEGRAL RESULTS:\n"
            text += "-" * 40 + "\n"
            # Single point analysis for Line Integral mode
            text += f"  B-Field at Center: {results['b_field_center']:.6f} T\n"
            text += f"  Current: {results['li_current']:.2f} A\n"
            text += f"  Voltage: {results['li_voltage']:.2f} V\n"
            text += f"  Power: {results['li_power']:.2f} W\n"
            text += f"  Total Turns: {results['turns_total']:.0f}\n"
            text += f"  Wire Length: {results['wire_length']:.2f} m\n"
            text += f"  Coil Weight: {results['coil_weight']:.2f} lbs\n"
            text += f"  Coil Spacing: {results['params']['coil_spacing']:.3f} m\n"
            
            # Add 2D field statistics
            if 'b_field_2d' in results:
                b_field_2d = results['b_field_2d']
                text += f"  Max B-Field in Grid: {np.max(b_field_2d):.6f} T\n"
                text += f"  Min B-Field in Grid: {np.min(b_field_2d):.6f} T\n"
                text += f"  Mean B-Field in Grid: {np.mean(b_field_2d):.6f} T\n"
        else:
            # Range analysis for First Order mode
            text += f"  Max B-Field (Power-based): {np.max(results['b_field']):.6f} T\n"
            text += f"  Max B-Field (Current-based): {np.max(results['b_field_current']):.6f} T\n"
            
            # Display max B-field data (pre-calculated)
            if 'turns_at_max_b' in results:
                text += f"  Turns at Max B-Field: {results['turns_at_max_b']:.0f}\n"
                text += f"  Weight at Max B-Field: {results['weight_at_max_b']:.2f} lbs\n"
            
            text += f"  Turns at Weight Limit: {results['turns_at_weight']:.0f}\n"
            weight_limit_index = min(results['turns_at_weight']-1, len(results['b_field'])-1)
            if weight_limit_index >= 0:
                text += f"  B-Field (Power) at Weight Limit: {results['b_field'][weight_limit_index]:.6f} T\n"
                text += f"  B-Field (Current) at Weight Limit: {results['b_field_current'][weight_limit_index]:.6f} T\n"
                text += f"  Current at Weight Limit: {results['current'][weight_limit_index]:.2f} A\n"
                text += f"  Voltage at Weight Limit: {results['voltage'][weight_limit_index]:.2f} V\n"
                text += f"  Power at Weight Limit: {results['power'][weight_limit_index]:.2f} W\n"
            else:
                text += f"  Weight limit exceeds simulation range\n"
        
        if results['mode'] == 'Line Integral':
            text += "\nComparison Table:\n"
            text += "Parameter\tFirst Order\tLine Integral\tDifference\n"
            text += "-" * 70 + "\n"
            
            # Get first order results for comparison
            first_order_results = self.run_first_order_simulation(results['params'])
            fo_turns = results['turns_total']
            fo_index = min(fo_turns-1, len(first_order_results['b_field'])-1)
            
            fo_b_field = first_order_results['b_field'][fo_index]
            fo_b_field_current = first_order_results['b_field_current'][fo_index]
            fo_current = first_order_results['current'][fo_index]
            fo_voltage = first_order_results['voltage'][fo_index]
            fo_power = first_order_results['power'][fo_index]
            
            text += f"Turns\t{fo_turns:.0f}\t{results['turns_total']:.0f}\t-\n"
            text += f"B-Field (Power) (T)\t{fo_b_field:.6f}\t{results['b_field_center']:.6f}\t{(results['b_field_center']-fo_b_field):.6f}\n"
            text += f"B-Field (Current) (T)\t{fo_b_field_current:.6f}\t{results['b_field_center']:.6f}\t{(results['b_field_center']-fo_b_field_current):.6f}\n"
            text += f"Current (A)\t{fo_current:.2f}\t{results['li_current']:.2f}\t{(results['li_current']-fo_current):.2f}\n"
            text += f"Voltage (V)\t{fo_voltage:.2f}\t{results['li_voltage']:.2f}\t{(results['li_voltage']-fo_voltage):.2f}\n"
            text += f"Power (W)\t{fo_power:.2f}\t{results['li_power']:.2f}\t{(results['li_power']-fo_power):.2f}\n"
            
            text += f"\nAdditional Line Integral Data:\n"
            text += f"Wire Length (m)\t{results['wire_length']:.2f}\n"
            text += f"Coil Weight (lbs)\t{results['coil_weight']:.2f}\n"
        else:
            text += "\nData Table:\n"
            text += "Turns\tB-Field(Power)(T)\tB-Field(Current)(T)\tCurrent(A)\tVoltage(V)\tPower(W)\n"
            text += "-" * 80 + "\n"
            
            for i in range(min(20, len(results['turns']))):  # Show first 20 points
                text += f"{results['turns'][i]:.0f}\t{results['b_field'][i]:.6f}\t{results['b_field_current'][i]:.6f}\t"
                text += f"{results['current'][i]:.2f}\t{results['voltage'][i]:.2f}\t"
                text += f"{results['power'][i]:.2f}\n"
            
        self.data_text.insert(1.0, text)
        
    def simulation_finished(self):
        """Handle simulation completion"""
        self.is_simulating = False
        self.first_order_button.config(state="normal")
        self.stop_button.config(state="disabled")
        self.status_var.set("Simulation completed")
        self.progress.stop()
        
    def stop_simulation(self):
        """Stop running simulation"""
        self.is_simulating = False
        self.status_var.set("Simulation stopped")
        self.progress.stop()
        self.run_button.config(state="normal")
        self.stop_button.config(state="disabled")
        
    def save_results(self):
        """Save results to file"""
        if not self.results_data:
            messagebox.showwarning("No Results", "No results to save!")
            return
            
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                # Create DataFrame
                df = pd.DataFrame({
                    'turns': self.results_data['turns'],
                    'b_field_power': self.results_data['b_field'],
                    'b_field_current': self.results_data['b_field_current'],
                    'current': self.results_data['current'],
                    'voltage': self.results_data['voltage'],
                    'power': self.results_data['power']
                })
                
                df.to_csv(filename, index=False)
                messagebox.showinfo("Success", f"Results saved to {filename}")
                
            except Exception as e:
                messagebox.showerror("Save Error", f"Error saving file: {str(e)}")
                
    def save_plots(self):
        """Save plots to PDF file"""
        if not self.results_data:
            messagebox.showwarning("No Results", "No plots to save!")
            return
            
        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                # Save both figures to PDF
                with plt.backends.backend_pdf.PdfPages(filename) as pdf:
                    # Save First Order plots
                    pdf.savefig(self.fig_first_order, bbox_inches='tight', dpi=300)
                    
                    # Save Field Analysis plots if they exist
                    if self.results_data['mode'] == 'Line Integral':
                        pdf.savefig(self.fig_field_analysis, bbox_inches='tight', dpi=300)
                
                messagebox.showinfo("Success", f"Plots saved to {filename}")
                
            except Exception as e:
                messagebox.showerror("Save Error", f"Error saving plots: {str(e)}")
                
    def clear_results(self):
        """Clear all results"""
        self.results_data = {}
        self.fig_first_order.clear()
        self.canvas_first_order.draw()
        self.fig_field_analysis.clear()
        self.canvas_field_analysis.draw()
        self.data_text.delete(1.0, tk.END)
        self.status_var.set("Results cleared")
        
        # Reset workflow state
        self.first_order_completed = False
        self.last_first_order_params = None
        self.field_analysis_button.config(state="disabled")


def main():
    """Main function to run the GUI"""
    root = tk.Tk()
    app = HelmholtzGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()