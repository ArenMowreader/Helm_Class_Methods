'''
Author: Aren Mowreader
Date: 7/2/2025
Version: 1.0
Description:

Script to create a GUI for the Electromagnet simulation.
    The GUI will have the following features:
        Inputs:
            Coil Parameters:
                -Inside coil radius (m)
                -Coil height (m) 
                -Radial depth (m)
                -Wire Gauge (AWG)
            Power Supply:
                -Power Supply Name (optional)
                -Max Power (W)
                -Max Current (A)
                -Max Voltage (V)
            Constraints:
                -Weight Limit (lbs)
            Field Specifications:
                -X Range (min, max)
                -Y Range (min, max)
                -Z Range (min, max)
                -Points per axis
                -Note: Entering the same value for x, y, or z dimensions will make the analysis on that plane
        Outputs:
            -Coil geometry plot
            -B-field vector plot
            -Parameter analysis plots:
                -Coil Length vs Turns
                -Weight vs Turns
                -Resistance vs Turns
                -Power vs Turns
                -Current vs Turns
                -Voltage vs Turns
                -B-field Average vs Turns
                -B-field Maximum vs Turns
                -B-field Uniformity vs Turns
            -Normalized summary plot
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
import threading
import time
import os
import sys

# Add parent directory to path to import electromagnet classes
sys.path.append(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from electromagnet import PowerSupply, Wire, Emag


class ElectromagnetGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Electromagnet Simulation GUI")
        self.root.geometry("1400x900")
        
        # Variables
        self.results_data = {}
        self.is_simulating = False
        self.simulation_completed = False
        
        # Create main containers
        self.create_widgets()
        self.setup_layout()
        
    def create_widgets(self):
        """Create all GUI widgets"""
        # Instructions
        self.create_instructions()
        
        # Input frames
        self.create_input_frames()
        
        # Control buttons
        self.create_control_buttons()
        
        # Results area
        self.create_results_area()
        
        # Status bar
        self.create_status_bar()
        
    def create_instructions(self):
        """Create instructions frame"""
        instructions_frame = ttk.LabelFrame(self.root, text="Instructions", padding="10")
        instructions_frame.grid(row=0, column=0, columnspan=2, sticky="ew", padx=10, pady=5)
        
        # Instructions
        ttk.Label(instructions_frame, text="1. Enter coil parameters (radius, height, depth, wire gauge)", 
                 font=("Arial", 10, "bold")).pack(anchor="w", pady=2)
        ttk.Label(instructions_frame, text="2. Set power supply specifications", 
                 font=("Arial", 10, "bold")).pack(anchor="w", pady=2)
        ttk.Label(instructions_frame, text="3. Define field specifications for B-field analysis", 
                 font=("Arial", 10, "bold")).pack(anchor="w", pady=2)
        ttk.Label(instructions_frame, text="4. Run simulation to generate plots and analysis", 
                 font=("Arial", 10, "bold")).pack(anchor="w", pady=2)
        
        # Note about plane analysis
        note_text = "Note: Entering the same value for x, y, or z dimensions will make the analysis on that plane"
        ttk.Label(instructions_frame, text=note_text, 
                 font=("Arial", 9, "italic"), foreground="blue").pack(anchor="w", pady=2)
        
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
        
        # Field specifications frame
        self.field_frame = ttk.LabelFrame(self.root, text="Field Specifications", padding="10")
        self.field_frame.grid(row=4, column=0, sticky="nsew", padx=10, pady=5)
        
        # Create input widgets
        self.create_coil_inputs()
        self.create_power_inputs()
        self.create_constraint_inputs()
        self.create_field_inputs()
        
    def create_coil_inputs(self):
        """Create coil parameter input widgets"""
        # Inside coil radius
        ttk.Label(self.coil_frame, text="Inside Coil Radius (m):").grid(row=0, column=0, sticky="w", pady=2)
        self.radius_var = tk.DoubleVar(value=0.125)
        self.radius_entry = ttk.Entry(self.coil_frame, textvariable=self.radius_var, width=15)
        self.radius_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        
        # Coil height
        ttk.Label(self.coil_frame, text="Coil Height (m):").grid(row=1, column=0, sticky="w", pady=2)
        self.height_var = tk.DoubleVar(value=0.0508)  # 2 inches
        self.height_entry = ttk.Entry(self.coil_frame, textvariable=self.height_var, width=15)
        self.height_entry.grid(row=1, column=1, sticky="w", padx=5, pady=2)
        
        # Radial depth
        ttk.Label(self.coil_frame, text="Radial Depth (m):").grid(row=2, column=0, sticky="w", pady=2)
        self.depth_var = tk.DoubleVar(value=0.0508)  # 2 inches
        self.depth_entry = ttk.Entry(self.coil_frame, textvariable=self.depth_var, width=15)
        self.depth_entry.grid(row=2, column=1, sticky="w", padx=5, pady=2)
        
        # Wire gauge
        ttk.Label(self.coil_frame, text="Wire Gauge (AWG):").grid(row=3, column=0, sticky="w", pady=2)
        self.wire_gauge_var = tk.IntVar(value=12)
        self.wire_gauge_entry = ttk.Entry(self.coil_frame, textvariable=self.wire_gauge_var, width=15)
        self.wire_gauge_entry.grid(row=3, column=1, sticky="w", padx=5, pady=2)
        
    def create_power_inputs(self):
        """Create power supply input widgets"""
        # Power Supply Name
        ttk.Label(self.power_frame, text="Power Supply Name: (optional)").grid(row=0, column=0, sticky="w", pady=2)
        self.power_supply_name_var = tk.StringVar(value="Keysight RP7952A")
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
        self.weight_limit_var = tk.DoubleVar(value=40.0)
        self.weight_limit_entry = ttk.Entry(self.constraints_frame, textvariable=self.weight_limit_var, width=15)
        self.weight_limit_entry.grid(row=0, column=1, sticky="w", padx=5, pady=2)
        
    def create_field_inputs(self):
        """Create field specification input widgets"""
        # X Range
        ttk.Label(self.field_frame, text="X Range (m):").grid(row=0, column=0, sticky="w", pady=2)
        self.x_min_var = tk.DoubleVar(value=-0.1)
        self.x_max_var = tk.DoubleVar(value=0.1)
        ttk.Entry(self.field_frame, textvariable=self.x_min_var, width=8).grid(row=0, column=1, sticky="w", padx=2, pady=2)
        ttk.Label(self.field_frame, text="to").grid(row=0, column=2, padx=2)
        ttk.Entry(self.field_frame, textvariable=self.x_max_var, width=8).grid(row=0, column=3, sticky="w", padx=2, pady=2)
        
        # Y Range
        ttk.Label(self.field_frame, text="Y Range (m):").grid(row=1, column=0, sticky="w", pady=2)
        self.y_min_var = tk.DoubleVar(value=-0.1)
        self.y_max_var = tk.DoubleVar(value=0.1)
        ttk.Entry(self.field_frame, textvariable=self.y_min_var, width=8).grid(row=1, column=1, sticky="w", padx=2, pady=2)
        ttk.Label(self.field_frame, text="to").grid(row=1, column=2, padx=2)
        ttk.Entry(self.field_frame, textvariable=self.y_max_var, width=8).grid(row=1, column=3, sticky="w", padx=2, pady=2)
        
        # Z Range
        ttk.Label(self.field_frame, text="Z Range (m):").grid(row=2, column=0, sticky="w", pady=2)
        self.z_min_var = tk.DoubleVar(value=0.125)
        self.z_max_var = tk.DoubleVar(value=0.125)
        ttk.Entry(self.field_frame, textvariable=self.z_min_var, width=8).grid(row=2, column=1, sticky="w", padx=2, pady=2)
        ttk.Label(self.field_frame, text="to").grid(row=2, column=2, padx=2)
        ttk.Entry(self.field_frame, textvariable=self.z_max_var, width=8).grid(row=2, column=3, sticky="w", padx=2, pady=2)
        
        # Points per axis
        ttk.Label(self.field_frame, text="Points per Axis:").grid(row=3, column=0, sticky="w", pady=2)
        self.points_per_axis_var = tk.IntVar(value=5)
        self.points_per_axis_entry = ttk.Entry(self.field_frame, textvariable=self.points_per_axis_var, width=15)
        self.points_per_axis_entry.grid(row=3, column=1, sticky="w", padx=5, pady=2)
        
    def create_control_buttons(self):
        """Create control button frame"""
        button_frame = ttk.Frame(self.root)
        button_frame.grid(row=5, column=0, sticky="ew", padx=10, pady=5)
        
        # Run simulation button
        self.run_button = ttk.Button(button_frame, text="Run Simulation", 
                                   command=self.run_simulation, style="Accent.TButton")
        self.run_button.pack(side="left", padx=5)
        
        # Stop simulation button
        self.stop_button = ttk.Button(button_frame, text="Stop", 
                                    command=self.stop_simulation, state="disabled")
        self.stop_button.pack(side="left", padx=5)
        
        # Clear results button
        self.clear_button = ttk.Button(button_frame, text="Clear Results", 
                                     command=self.clear_results)
        self.clear_button.pack(side="left", padx=5)
        
        # Export results button
        self.export_button = ttk.Button(button_frame, text="Export Results", 
                                      command=self.export_results, state="disabled")
        self.export_button.pack(side="left", padx=5)
        
    def create_results_area(self):
        """Create results display area"""
        results_frame = ttk.LabelFrame(self.root, text="Results", padding="10")
        results_frame.grid(row=1, column=1, rowspan=4, sticky="nsew", padx=10, pady=5)
        
        # Create notebook for tabbed results
        self.notebook = ttk.Notebook(results_frame)
        self.notebook.pack(fill="both", expand=True)
        
        # Geometry tab
        self.geometry_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.geometry_frame, text="Coil Geometry")
        
        # B-field tab
        self.bfield_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.bfield_frame, text="B-Field")
        
        # Analysis tab
        self.analysis_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.analysis_frame, text="Parameter Analysis")
        
        # Summary tab
        self.summary_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.summary_frame, text="Summary")
        
        # Console tab
        self.console_frame = ttk.Frame(self.notebook)
        self.notebook.add(self.console_frame, text="Console")
        
        # Create console
        self.console = scrolledtext.ScrolledText(self.console_frame, height=20, width=80)
        self.console.pack(fill="both", expand=True)
        
    def create_status_bar(self):
        """Create status bar"""
        self.status_var = tk.StringVar(value="Ready")
        self.status_bar = ttk.Label(self.root, textvariable=self.status_var, relief="sunken")
        self.status_bar.grid(row=6, column=0, columnspan=2, sticky="ew", padx=10, pady=2)
        
    def setup_layout(self):
        """Set up grid layout weights"""
        self.root.grid_rowconfigure(1, weight=1)
        self.root.grid_rowconfigure(2, weight=1)
        self.root.grid_rowconfigure(3, weight=1)
        self.root.grid_rowconfigure(4, weight=1)
        self.root.grid_columnconfigure(1, weight=1)
        
    def log_message(self, message):
        """Add message to console"""
        self.console.insert(tk.END, f"{message}\n")
        self.console.see(tk.END)
        self.root.update_idletasks()
        
    def run_simulation(self):
        """Run the electromagnet simulation"""
        if self.is_simulating:
            return
            
        # Validate inputs
        if not self.validate_inputs():
            return
            
        # Start simulation in separate thread
        self.is_simulating = True
        self.run_button.config(state="disabled")
        self.stop_button.config(state="normal")
        self.status_var.set("Running simulation...")
        
        # Clear previous results
        self.clear_plots()
        
        # Start simulation thread
        self.simulation_thread = threading.Thread(target=self._run_simulation_thread)
        self.simulation_thread.daemon = True
        self.simulation_thread.start()
        
    def _run_simulation_thread(self):
        """Run simulation in separate thread"""
        try:
            self.log_message("Starting electromagnet simulation...")
            
            # Get parameters
            params = self.get_parameters()
            self.log_message(f"Parameters: {params}")
            
            # Create electromagnet
            self.log_message("Creating electromagnet...")
            electromagnet = self.create_electromagnet(params)
            
            # Calculate turns per layer
            turns_per_layer = int(np.floor(params['depth'] / electromagnet.wire.diameter_nom_m))
            self.log_message(f"Turns per layer: {turns_per_layer}")
            
            # Define field
            self.log_message("Defining field...")
            electromagnet.def_field(
                x_range=(params['x_min'], params['x_max']),
                y_range=(params['y_min'], params['y_max']),
                z_range=(params['z_min'], params['z_max']),
                points_per_axis=params['points_per_axis']
            )
            
            # Initialize arrays
            max_turns = turns_per_layer * turns_per_layer
            current = np.zeros(max_turns)
            resistance = np.zeros(max_turns)
            voltage = np.zeros(max_turns)
            power = np.zeros(max_turns)
            length = np.zeros(max_turns)
            weight = np.zeros(max_turns)
            b_field = np.zeros((max_turns, electromagnet.field.shape[0], 3))
            b_avg = np.zeros(max_turns)
            b_max = np.zeros(max_turns)
            turns = np.zeros(max_turns)
            
            # Fill arrays layer by layer
            k = 0
            for j in range(turns_per_layer):  # For each layer
                for i in range(turns_per_layer):  # For each turn in the layer
                    if self.is_simulating == False:  # Check if stopped
                        break
                        
                    # Check if adding this turn would exceed weight limit
                    if electromagnet.weight >= params['weight_limit']:
                        break
                        
                    self.log_message(f"Adding turn {k+1}: layer {j+1}, turn {i+1}")
                    
                    electromagnet.add_turn(
                        center=np.array([0, 0, 0 - i * electromagnet.wire.diameter_nom_m]),
                        radius=(params['radius'] / 2) + (j * electromagnet.wire.diameter_nom_m),
                        points_per_turn=500
                    )
                    electromagnet.calc_b_field()
                    
                    current[k] = electromagnet.current
                    weight[k] = electromagnet.weight
                    voltage[k] = electromagnet.voltage
                    power[k] = electromagnet.power
                    resistance[k] = electromagnet.resistance
                    length[k] = electromagnet.coil_length
                    b_field[k] = electromagnet.b_field
                    b_avg[k] = np.mean(np.linalg.norm(electromagnet.b_field, axis=1))
                    b_max[k] = np.max(np.linalg.norm(electromagnet.b_field, axis=1))
                    turns[k] = electromagnet.net_turns
                    k += 1
                
                # If we've exceeded weight limit, break out of outer loop too
                if electromagnet.weight >= params['weight_limit']:
                    break
                    
                if self.is_simulating == False:  # Check if stopped
                    break
            
            # Store results
            self.results_data = {
                'electromagnet': electromagnet,
                'current': current[:k],
                'weight': weight[:k],
                'voltage': voltage[:k],
                'power': power[:k],
                'resistance': resistance[:k],
                'length': length[:k],
                'b_field': b_field[:k],
                'b_avg': b_avg[:k],
                'b_max': b_max[:k],
                'turns': turns[:k],
                'params': params
            }
            
            # Generate plots
            self.log_message("Generating plots...")
            self.generate_plots()
            
            self.simulation_completed = True
            self.log_message("Simulation completed successfully!")
            
        except Exception as e:
            self.log_message(f"Error during simulation: {str(e)}")
            messagebox.showerror("Simulation Error", f"An error occurred during simulation:\n{str(e)}")
        finally:
            self.is_simulating = False
            self.run_button.config(state="normal")
            self.stop_button.config(state="disabled")
            self.export_button.config(state="normal")
            self.status_var.set("Ready")
            
    def validate_inputs(self):
        """Validate user inputs"""
        try:
            # Check coil parameters
            if self.radius_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Inside coil radius must be positive")
                return False
                
            if self.height_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Coil height must be positive")
                return False
                
            if self.depth_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Radial depth must be positive")
                return False
                
            if self.wire_gauge_var.get() < 0 or self.wire_gauge_var.get() > 50:
                messagebox.showerror("Invalid Input", "Wire gauge must be between 0 and 50")
                return False
                
            # Check power supply parameters
            if self.max_power_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Max power must be positive")
                return False
                
            if self.max_current_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Max current must be positive")
                return False
                
            if self.max_voltage_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Max voltage must be positive")
                return False
                
            # Check constraints
            if self.weight_limit_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Weight limit must be positive")
                return False
                
            # Check field specifications
            if self.x_min_var.get() >= self.x_max_var.get():
                messagebox.showerror("Invalid Input", "X min must be less than X max")
                return False
                
            if self.y_min_var.get() >= self.y_max_var.get():
                messagebox.showerror("Invalid Input", "Y min must be less than Y max")
                return False
                
            if self.z_min_var.get() > self.z_max_var.get():
                messagebox.showerror("Invalid Input", "Z min must be less than or equal to Z max")
                return False
                
            if self.points_per_axis_var.get() <= 0:
                messagebox.showerror("Invalid Input", "Points per axis must be positive")
                return False
                
            return True
            
        except tk.TclError:
            messagebox.showerror("Invalid Input", "Please check that all inputs are valid numbers")
            return False
            
    def get_parameters(self):
        """Get all parameters from GUI"""
        return {
            'radius': self.radius_var.get(),
            'height': self.height_var.get(),
            'depth': self.depth_var.get(),
            'wire_gauge': self.wire_gauge_var.get(),
            'power_supply_name': self.power_supply_name_var.get(),
            'max_power': self.max_power_var.get(),
            'max_current': self.max_current_var.get(),
            'max_voltage': self.max_voltage_var.get(),
            'weight_limit': self.weight_limit_var.get(),
            'x_min': self.x_min_var.get(),
            'x_max': self.x_max_var.get(),
            'y_min': self.y_min_var.get(),
            'y_max': self.y_max_var.get(),
            'z_min': self.z_min_var.get(),
            'z_max': self.z_max_var.get(),
            'points_per_axis': self.points_per_axis_var.get()
        }
        
    def create_electromagnet(self, params):
        """Create electromagnet object from parameters"""
        power_supply = PowerSupply(
            name=params['power_supply_name'],
            max_voltage=params['max_voltage'],
            max_current=params['max_current'],
            max_power=params['max_power']
        )
        
        wire = Wire(AWG=params['wire_gauge'])
        
        electromagnet = Emag(power_supply=power_supply, wire=wire)
        
        return electromagnet
        
    def generate_plots(self):
        """Generate all plots"""
        if not self.simulation_completed:
            return
            
        # Plot coil geometry
        self.plot_coil_geometry()
        
        # Plot B-field
        self.plot_b_field()
        
        # Plot parameter analysis
        self.plot_parameter_analysis()
        
        # Plot summary
        self.plot_summary()
        
    def plot_coil_geometry(self):
        """Plot coil geometry"""
        if not self.simulation_completed:
            return
            
        # Clear previous plot
        for widget in self.geometry_frame.winfo_children():
            widget.destroy()
            
        # Create figure
        fig = Figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot coil geometry
        self.results_data['electromagnet'].plot_coil_geometry(ax=ax)
        
        # Create canvas
        canvas = FigureCanvasTkAgg(fig, self.geometry_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        
    def plot_b_field(self):
        """Plot B-field"""
        if not self.simulation_completed:
            return
            
        # Clear previous plot
        for widget in self.bfield_frame.winfo_children():
            widget.destroy()
            
        # Create figure
        fig = Figure(figsize=(8, 6))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot B-field
        self.results_data['electromagnet'].plot_b_field(ax=ax)
        
        # Create canvas
        canvas = FigureCanvasTkAgg(fig, self.bfield_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        
    def plot_parameter_analysis(self):
        """Plot parameter analysis"""
        if not self.simulation_completed:
            return
            
        # Clear previous plot
        for widget in self.analysis_frame.winfo_children():
            widget.destroy()
            
        # Get data
        turns = self.results_data['turns']
        valid_mask = turns > 0
        
        if np.sum(valid_mask) == 0:
            self.log_message("No valid data points for plotting")
            return
            
        # Create figure with subplots
        fig, axes = plt.subplots(3, 3, figsize=(15, 12))
        fig.suptitle('Electromagnet Parameters vs Number of Turns', fontsize=16)
        
        # Plot 1: Coil Length
        axes[0, 0].plot(turns[valid_mask], self.results_data['length'][valid_mask], 'b-', linewidth=2)
        axes[0, 0].set_xlabel('Number of Turns')
        axes[0, 0].set_ylabel('Coil Length (m)')
        axes[0, 0].set_title('Coil Length vs Turns')
        axes[0, 0].grid(True)
        
        # Plot 2: Weight
        axes[0, 1].plot(turns[valid_mask], self.results_data['weight'][valid_mask], 'r-', linewidth=2)
        axes[0, 1].set_xlabel('Number of Turns')
        axes[0, 1].set_ylabel('Weight (lbs)')
        axes[0, 1].set_title('Weight vs Turns')
        axes[0, 1].grid(True)
        
        # Plot 3: Resistance
        axes[0, 2].plot(turns[valid_mask], self.results_data['resistance'][valid_mask], 'g-', linewidth=2)
        axes[0, 2].set_xlabel('Number of Turns')
        axes[0, 2].set_ylabel('Resistance (Ω)')
        axes[0, 2].set_title('Resistance vs Turns')
        axes[0, 2].grid(True)
        
        # Plot 4: Power
        axes[1, 0].plot(turns[valid_mask], self.results_data['power'][valid_mask], 'm-', linewidth=2)
        axes[1, 0].set_xlabel('Number of Turns')
        axes[1, 0].set_ylabel('Power (W)')
        axes[1, 0].set_title('Power vs Turns')
        axes[1, 0].grid(True)
        
        # Plot 5: Current
        axes[1, 1].plot(turns[valid_mask], self.results_data['current'][valid_mask], 'c-', linewidth=2)
        axes[1, 1].set_xlabel('Number of Turns')
        axes[1, 1].set_ylabel('Current (A)')
        axes[1, 1].set_title('Current vs Turns')
        axes[1, 1].grid(True)
        
        # Plot 6: Voltage
        axes[1, 2].plot(turns[valid_mask], self.results_data['voltage'][valid_mask], 'orange', linewidth=2)
        axes[1, 2].set_xlabel('Number of Turns')
        axes[1, 2].set_ylabel('Voltage (V)')
        axes[1, 2].set_title('Voltage vs Turns')
        axes[1, 2].grid(True)
        
        # Plot 7: B-field Average
        axes[2, 0].plot(turns[valid_mask], self.results_data['b_avg'][valid_mask], 'purple', linewidth=2)
        axes[2, 0].set_xlabel('Number of Turns')
        axes[2, 0].set_ylabel('B-field Average (T)')
        axes[2, 0].set_title('B-field Average vs Turns')
        axes[2, 0].grid(True)
        
        # Plot 8: B-field Maximum
        axes[2, 1].plot(turns[valid_mask], self.results_data['b_max'][valid_mask], 'brown', linewidth=2)
        axes[2, 1].set_xlabel('Number of Turns')
        axes[2, 1].set_ylabel('B-field Maximum (T)')
        axes[2, 1].set_title('B-field Maximum vs Turns')
        axes[2, 1].grid(True)
        
        # Plot 9: B-field Ratio (Max/Avg)
        b_ratio = self.results_data['b_max'][valid_mask] / self.results_data['b_avg'][valid_mask]
        axes[2, 2].plot(turns[valid_mask], b_ratio, 'teal', linewidth=2)
        axes[2, 2].set_xlabel('Number of Turns')
        axes[2, 2].set_ylabel('B-field Ratio (Max/Avg)')
        axes[2, 2].set_title('B-field Uniformity vs Turns')
        axes[2, 2].grid(True)
        
        plt.tight_layout()
        
        # Create canvas
        canvas = FigureCanvasTkAgg(fig, self.analysis_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        
    def plot_summary(self):
        """Plot normalized summary"""
        if not self.simulation_completed:
            return
            
        # Clear previous plot
        for widget in self.summary_frame.winfo_children():
            widget.destroy()
            
        # Get data
        turns = self.results_data['turns']
        valid_mask = turns > 0
        
        if np.sum(valid_mask) == 0:
            return
            
        # Create figure
        fig, ax = plt.subplots(figsize=(12, 8))
        
        # Normalize all parameters to their maximum values for comparison
        length_norm = self.results_data['length'][valid_mask] / np.max(self.results_data['length'][valid_mask])
        weight_norm = self.results_data['weight'][valid_mask] / np.max(self.results_data['weight'][valid_mask])
        resistance_norm = self.results_data['resistance'][valid_mask] / np.max(self.results_data['resistance'][valid_mask])
        power_norm = self.results_data['power'][valid_mask] / np.max(self.results_data['power'][valid_mask])
        current_norm = self.results_data['current'][valid_mask] / np.max(self.results_data['current'][valid_mask])
        voltage_norm = self.results_data['voltage'][valid_mask] / np.max(self.results_data['voltage'][valid_mask])
        b_avg_norm = self.results_data['b_avg'][valid_mask] / np.max(self.results_data['b_avg'][valid_mask])
        b_max_norm = self.results_data['b_max'][valid_mask] / np.max(self.results_data['b_max'][valid_mask])
        
        ax.plot(turns[valid_mask], length_norm, 'b-', linewidth=2, label='Coil Length')
        ax.plot(turns[valid_mask], weight_norm, 'r-', linewidth=2, label='Weight')
        ax.plot(turns[valid_mask], resistance_norm, 'g-', linewidth=2, label='Resistance')
        ax.plot(turns[valid_mask], power_norm, 'm-', linewidth=2, label='Power')
        ax.plot(turns[valid_mask], current_norm, 'c-', linewidth=2, label='Current')
        ax.plot(turns[valid_mask], voltage_norm, 'orange', linewidth=2, label='Voltage')
        ax.plot(turns[valid_mask], b_avg_norm, 'purple', linewidth=2, label='B-field Average')
        ax.plot(turns[valid_mask], b_max_norm, 'brown', linewidth=2, label='B-field Maximum')
        
        ax.set_xlabel('Number of Turns')
        ax.set_ylabel('Normalized Value')
        ax.set_title('All Parameters vs Turns (Normalized)')
        ax.legend()
        ax.grid(True)
        
        # Create canvas
        canvas = FigureCanvasTkAgg(fig, self.summary_frame)
        canvas.draw()
        canvas.get_tk_widget().pack(fill="both", expand=True)
        
    def clear_plots(self):
        """Clear all plots"""
        for frame in [self.geometry_frame, self.bfield_frame, self.analysis_frame, self.summary_frame]:
            for widget in frame.winfo_children():
                widget.destroy()
                
    def stop_simulation(self):
        """Stop the simulation"""
        self.is_simulating = False
        self.log_message("Simulation stopped by user")
        
    def clear_results(self):
        """Clear all results"""
        self.results_data = {}
        self.simulation_completed = False
        self.clear_plots()
        self.console.delete(1.0, tk.END)
        self.export_button.config(state="disabled")
        self.log_message("Results cleared")
        
    def export_results(self):
        """Export results to file"""
        if not self.simulation_completed:
            messagebox.showwarning("No Results", "No results to export. Please run a simulation first.")
            return
            
        # Ask user for file location
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")],
            title="Export Results"
        )
        
        if filename:
            try:
                # Create DataFrame
                data = {
                    'Turns': self.results_data['turns'],
                    'Coil_Length_m': self.results_data['length'],
                    'Weight_lbs': self.results_data['weight'],
                    'Resistance_Ohm': self.results_data['resistance'],
                    'Power_W': self.results_data['power'],
                    'Current_A': self.results_data['current'],
                    'Voltage_V': self.results_data['voltage'],
                    'B_Field_Avg_T': self.results_data['b_avg'],
                    'B_Field_Max_T': self.results_data['b_max']
                }
                
                df = pd.DataFrame(data)
                
                # Export to CSV
                df.to_csv(filename, index=False)
                
                self.log_message(f"Results exported to {filename}")
                messagebox.showinfo("Export Successful", f"Results exported to {filename}")
                
            except Exception as e:
                self.log_message(f"Error exporting results: {str(e)}")
                messagebox.showerror("Export Error", f"Error exporting results:\n{str(e)}")


def main():
    """Main function to launch the GUI"""
    root = tk.Tk()
    app = ElectromagnetGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
