'''
Modern Helmholtz Coil Optimization GUI

A comprehensive GUI for Helmholtz coil design with 5 main stages:
1. First Order Optimization
2. First Order Visualization  
3. Geometric Optimization
4. Field Visualization
5. Design Space Visualization

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

class ModernHelmholtzGUI:
    def __init__(self, root):
        self.root = root
        self.root.title("Modern Helmholtz Coil Optimization Suite")
        self.root.geometry("1400x900")
        self.root.configure(bg='#f0f0f0')
        
        # Style configuration
        self.setup_styles()
        
        # Data storage
        self.first_order_results = None
        self.geometric_results = None
        self.field_results = None
        self.design_space_results = None
        
        # Create main interface
        self.create_main_interface()
        
    def setup_styles(self):
        """Configure modern styling"""
        style = ttk.Style()
        style.theme_use('clam')
        
        # Configure colors
        style.configure('Title.TLabel', font=('Arial', 16, 'bold'), foreground='#2c3e50')
        style.configure('Stage.TLabel', font=('Arial', 12, 'bold'), foreground='#34495e')
        style.configure('Info.TLabel', font=('Arial', 10), foreground='#7f8c8d')
        style.configure('Success.TLabel', font=('Arial', 10, 'bold'), foreground='#27ae60')
        style.configure('Warning.TLabel', font=('Arial', 10, 'bold'), foreground='#e67e22')
        
        # Configure buttons
        style.configure('Action.TButton', font=('Arial', 10, 'bold'), padding=(10, 5))
        style.configure('Primary.TButton', font=('Arial', 10, 'bold'), padding=(10, 5))
        
    def create_main_interface(self):
        """Create the main interface with 5 stages"""
        # Main container
        main_frame = ttk.Frame(self.root, padding="10")
        main_frame.pack(fill="both", expand=True)
        
        # Title
        title_label = ttk.Label(main_frame, text="Helmholtz Coil Optimization Suite", 
                               style='Title.TLabel')
        title_label.pack(pady=(0, 20))
        
        # Create notebook for stages
        self.notebook = ttk.Notebook(main_frame)
        self.notebook.pack(fill="both", expand=True)
        
        # Create each stage
        self.create_first_order_optimization_tab()
        self.create_first_order_visualization_tab()
        self.create_geometric_optimization_tab()
        self.create_field_visualization_tab()
        self.create_design_space_visualization_tab()
        
    def create_first_order_optimization_tab(self):
        """Create First Order Optimization tab"""
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text="1. First Order Optimization")
        
        # Main content frame
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Left panel - Input parameters
        left_panel = ttk.LabelFrame(content_frame, text="Optimization Parameters", padding="15")
        left_panel.pack(side="left", fill="y", padx=(0, 10))
        
        # Required parameters
        required_frame = ttk.LabelFrame(left_panel, text="Required Parameters", padding="10")
        required_frame.pack(fill="x", pady=(0, 15))
        
        # Weight limit
        ttk.Label(required_frame, text="Weight Limit (lbs):").grid(row=0, column=0, sticky="w", pady=5)
        self.weight_var = tk.DoubleVar(value=75.0)
        weight_entry = ttk.Entry(required_frame, textvariable=self.weight_var, width=15)
        weight_entry.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Radius
        ttk.Label(required_frame, text="Radius (m):").grid(row=1, column=0, sticky="w", pady=5)
        self.radius_var = tk.DoubleVar(value=0.125)
        radius_entry = ttk.Entry(required_frame, textvariable=self.radius_var, width=15)
        radius_entry.grid(row=1, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Max Power
        ttk.Label(required_frame, text="Max Power (W):").grid(row=2, column=0, sticky="w", pady=5)
        self.max_power_var = tk.DoubleVar(value=10000.0)
        max_power_entry = ttk.Entry(required_frame, textvariable=self.max_power_var, width=15)
        max_power_entry.grid(row=2, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Optional parameters
        optional_frame = ttk.LabelFrame(left_panel, text="Optional Parameters (Leave blank to optimize)", padding="10")
        optional_frame.pack(fill="x", pady=(0, 15))
        
        # Max Current
        ttk.Label(optional_frame, text="Max Current (A):").grid(row=0, column=0, sticky="w", pady=5)
        self.max_current_var = tk.StringVar(value="")
        current_entry = ttk.Entry(optional_frame, textvariable=self.max_current_var, width=15)
        current_entry.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Max Voltage
        ttk.Label(optional_frame, text="Max Voltage (V):").grid(row=1, column=0, sticky="w", pady=5)
        self.max_voltage_var = tk.StringVar(value="")
        voltage_entry = ttk.Entry(optional_frame, textvariable=self.max_voltage_var, width=15)
        voltage_entry.grid(row=1, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Wire Gauge
        ttk.Label(optional_frame, text="Wire Gauge (AWG):").grid(row=2, column=0, sticky="w", pady=5)
        self.wire_gauge_var = tk.StringVar(value="")
        wire_gauge_entry = ttk.Entry(optional_frame, textvariable=self.wire_gauge_var, width=15)
        wire_gauge_entry.grid(row=2, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Optimization note
        note_text = "Note: Leave optional parameters blank to optimize over them.\nFill them in to fix those parameters during optimization."
        note_label = ttk.Label(optional_frame, text=note_text, style='Info.TLabel', wraplength=200)
        note_label.grid(row=3, column=0, columnspan=2, sticky="w", pady=(10, 0))
        
        # Control buttons
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill="x", pady=(15, 0))
        
        self.optimize_button = ttk.Button(button_frame, text="Run Optimization", 
                                        command=self.run_first_order_optimization,
                                        style='Primary.TButton')
        self.optimize_button.pack(side="left", padx=(0, 10))
        
        self.clear_button = ttk.Button(button_frame, text="Clear Results", 
                                     command=self.clear_first_order_results,
                                     style='Action.TButton')
        self.clear_button.pack(side="left")
        
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

3. Click 'Run Optimization' to find optimal power supply and wire specifications

4. Results will show optimal operating point and performance metrics"""
        
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
            self.root.after(0, lambda: messagebox.showerror("Optimization Error", str(e)))
    
    def _display_optimization_results(self, result):
        """Display optimization results"""
        self.results_text.delete(1.0, tk.END)
        
        if result['success']:
            self.first_order_results = result
            
            # Display results
            self.results_text.insert(tk.END, "OPTIMIZATION SUCCESSFUL\n")
            self.results_text.insert(tk.END, "=" * 50 + "\n\n")
            
            # Power supply specs
            ps = result['optimal_power_supply']
            self.results_text.insert(tk.END, "OPTIMAL POWER SUPPLY:\n")
            self.results_text.insert(tk.END, f"  Max Power: {ps['Pmax']:.0f} W\n")
            self.results_text.insert(tk.END, f"  Max Current: {ps['Imax']:.1f} A\n")
            self.results_text.insert(tk.END, f"  Max Voltage: {ps['Vmax']:.1f} V\n\n")
            
            # Wire specs
            self.results_text.insert(tk.END, "OPTIMAL WIRE:\n")
            self.results_text.insert(tk.END, f"  AWG Size: {result['optimal_awg']:.1f} (use AWG{round(result['optimal_awg'])})\n\n")
            
            # Operating point
            op = result['operating_point']
            self.results_text.insert(tk.END, "OPTIMAL OPERATING POINT:\n")
            self.results_text.insert(tk.END, f"  Turns per Coil: {op['turns_per_coil']:.0f}\n")
            self.results_text.insert(tk.END, f"  Current: {op['current_A']:.1f} A\n")
            self.results_text.insert(tk.END, f"  Voltage: {op['voltage_V']:.1f} V\n")
            self.results_text.insert(tk.END, f"  Power: {op['power_W']:.0f} W\n")
            self.results_text.insert(tk.END, f"  Single Coil Weight: {op['weight_single_coil_lbs']:.1f} lbs\n")
            self.results_text.insert(tk.END, f"  Both Coils Weight: {op['weight_both_coils_lbs']:.1f} lbs\n")
            self.results_text.insert(tk.END, f"  Wire Length: {op['wire_length_m']:.1f} m\n\n")
            
            # Performance
            self.results_text.insert(tk.END, "PERFORMANCE:\n")
            self.results_text.insert(tk.END, f"  Maximum B-field: {result['max_b_field_tesla']:.4f} T\n")
            
            # Enable next stage
            self.notebook.tab(1, state="normal")  # Enable visualization tab
            
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
        self.notebook.tab(1, state="disabled")  # Disable visualization tab
    
    def create_first_order_visualization_tab(self):
        """Create First Order Visualization tab"""
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text="2. First Order Visualization")
        self.notebook.tab(1, state="disabled")  # Initially disabled
        
        # Main content frame
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Left panel - Controls
        left_panel = ttk.LabelFrame(content_frame, text="Visualization Controls", padding="15")
        left_panel.pack(side="left", fill="y", padx=(0, 10))
        
        # Directions
        directions_frame = ttk.LabelFrame(left_panel, text="Directions", padding="10")
        directions_frame.pack(fill="x", pady=(0, 15))
        
        directions_text = """First Order Visualization Stage:

1. Run First Order Optimization first to enable this stage

2. This stage shows detailed performance curves:
   • B-field vs Turns
   • Current vs Turns  
   • Voltage vs Turns
   • Power vs Turns

3. Use the toolbar to:
   • Zoom and pan
   • Save plots
   • Adjust view

4. Results can be saved as CSV or plots as PDF"""
        
        directions_label = ttk.Label(directions_frame, text=directions_text, style='Info.TLabel', 
                                   justify="left", wraplength=200)
        directions_label.pack(anchor="w")
        
        # Control buttons
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill="x", pady=(15, 0))
        
        self.plot_button = ttk.Button(button_frame, text="Generate Plots", 
                                    command=self.generate_first_order_plots,
                                    style='Primary.TButton')
        self.plot_button.pack(side="left", padx=(0, 10))
        
        self.save_data_button = ttk.Button(button_frame, text="Save Data", 
                                         command=self.save_first_order_data,
                                         style='Action.TButton')
        self.save_data_button.pack(side="left", padx=(0, 10))
        
        self.save_plots_button = ttk.Button(button_frame, text="Save Plots", 
                                          command=self.save_first_order_plots,
                                          style='Action.TButton')
        self.save_plots_button.pack(side="left")
        
        # Right panel - Plots
        right_panel = ttk.Frame(content_frame)
        right_panel.pack(side="right", fill="both", expand=True)
        
        # Create matplotlib figure
        self.fig_first_order = Figure(figsize=(10, 8), dpi=100)
        self.canvas_first_order = FigureCanvasTkAgg(self.fig_first_order, right_panel)
        self.canvas_first_order.get_tk_widget().pack(fill="both", expand=True)
        
        # Add navigation toolbar
        toolbar_frame = ttk.Frame(right_panel)
        toolbar_frame.pack(fill="x")
        self.toolbar_first_order = NavigationToolbar2Tk(self.canvas_first_order, toolbar_frame)
        self.toolbar_first_order.update()
    
    def generate_first_order_plots(self):
        """Generate first order visualization plots"""
        if self.first_order_results is None:
            messagebox.showwarning("No Data", "Please run First Order Optimization first")
            return
            
        try:
            # Clear previous plots
            self.fig_first_order.clear()
            
            # Get performance data
            performance = self.first_order_results['full_performance']
            
            # Create 2x2 subplot layout
            gs = self.fig_first_order.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
            
            # B-field vs Turns
            ax1 = self.fig_first_order.add_subplot(gs[0, 0])
            ax1.plot(performance['turns_single'], performance['b_field'], 'b-', linewidth=2, label='B-Field')
            ax1.set_xlabel('Turns per Coil')
            ax1.set_ylabel('B-Field (T)')
            ax1.set_title('B-Field vs Turns')
            ax1.grid(True, alpha=0.3)
            ax1.legend()
            
            # Current vs Turns
            ax2 = self.fig_first_order.add_subplot(gs[0, 1])
            ax2.plot(performance['turns_single'], performance['current'], 'r-', linewidth=2, label='Current')
            ax2.set_xlabel('Turns per Coil')
            ax2.set_ylabel('Current (A)')
            ax2.set_title('Current vs Turns')
            ax2.grid(True, alpha=0.3)
            ax2.legend()
            
            # Voltage vs Turns
            ax3 = self.fig_first_order.add_subplot(gs[1, 0])
            ax3.plot(performance['turns_single'], performance['voltage'], 'g-', linewidth=2, label='Voltage')
            ax3.set_xlabel('Turns per Coil')
            ax3.set_ylabel('Voltage (V)')
            ax3.set_title('Voltage vs Turns')
            ax3.grid(True, alpha=0.3)
            ax3.legend()
            
            # Power vs Turns
            ax4 = self.fig_first_order.add_subplot(gs[1, 1])
            ax4.plot(performance['turns_single'], performance['power'], 'm-', linewidth=2, label='Power')
            ax4.set_xlabel('Turns per Coil')
            ax4.set_ylabel('Power (W)')
            ax4.set_title('Power vs Turns')
            ax4.grid(True, alpha=0.3)
            ax4.legend()
            
            # Add optimal operating point markers
            op = self.first_order_results['operating_point']
            optimal_turns = op['turns_per_coil']
            
            # Find index for optimal turns
            turns_array = performance['turns_single']
            optimal_idx = np.argmin(np.abs(turns_array - optimal_turns))
            
            ax1.axvline(x=optimal_turns, color='red', linestyle='--', alpha=0.7, label='Optimal')
            ax2.axvline(x=optimal_turns, color='red', linestyle='--', alpha=0.7, label='Optimal')
            ax3.axvline(x=optimal_turns, color='red', linestyle='--', alpha=0.7, label='Optimal')
            ax4.axvline(x=optimal_turns, color='red', linestyle='--', alpha=0.7, label='Optimal')
            
            # Update legends
            ax1.legend()
            ax2.legend()
            ax3.legend()
            ax4.legend()
            
            # Set main title
            self.fig_first_order.suptitle('First Order Optimization Results', fontsize=14, fontweight='bold')
            
            # Refresh canvas
            self.canvas_first_order.draw()
            
        except Exception as e:
            messagebox.showerror("Plot Error", f"Error generating plots: {str(e)}")
    
    def save_first_order_data(self):
        """Save first order data to CSV"""
        if self.first_order_results is None:
            messagebox.showwarning("No Data", "No data to save")
            return
            
        filename = filedialog.asksaveasfilename(
            defaultextension=".csv",
            filetypes=[("CSV files", "*.csv"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                performance = self.first_order_results['full_performance']
                
                # Create DataFrame
                df = pd.DataFrame({
                    'turns_per_coil': performance['turns_single'],
                    'b_field_tesla': performance['b_field'],
                    'current_amps': performance['current'],
                    'voltage_volts': performance['voltage'],
                    'power_watts': performance['power'],
                    'weight_single_coil_lbs': performance['weights'],
                    'wire_length_m': performance['wire_length']
                })
                
                df.to_csv(filename, index=False)
                messagebox.showinfo("Success", f"Data saved to {filename}")
                
            except Exception as e:
                messagebox.showerror("Save Error", f"Error saving data: {str(e)}")
    
    def save_first_order_plots(self):
        """Save first order plots to PDF"""
        if self.first_order_results is None:
            messagebox.showwarning("No Data", "No plots to save")
            return
            
        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                self.fig_first_order.savefig(filename, bbox_inches='tight', dpi=300)
                messagebox.showinfo("Success", f"Plots saved to {filename}")
                
            except Exception as e:
                messagebox.showerror("Save Error", f"Error saving plots: {str(e)}")
    
    def create_geometric_optimization_tab(self):
        """Create Geometric Optimization tab"""
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text="3. Geometric Optimization")
        self.notebook.tab(2, state="disabled")  # Initially disabled
        
        # Main content frame
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Left panel - Parameters
        left_panel = ttk.LabelFrame(content_frame, text="Geometric Optimization Parameters", padding="15")
        left_panel.pack(side="left", fill="y", padx=(0, 10))
        
        # Directions
        directions_frame = ttk.LabelFrame(left_panel, text="Directions", padding="10")
        directions_frame.pack(fill="x", pady=(0, 15))
        
        directions_text = """Geometric Optimization Stage:

1. Run First Order Optimization first to get power supply and wire specs

2. This stage optimizes coil geometry:
   • Radial layers (depth)
   • Axial layers (height)
   • Uses power supply specs from Stage 1

3. Parameters are automatically set from First Order results

4. Click 'Run Geometric Optimization' to find optimal coil geometry

5. Results show optimal layer configuration and performance"""
        
        directions_label = ttk.Label(directions_frame, text=directions_text, style='Info.TLabel', 
                                   justify="left", wraplength=200)
        directions_label.pack(anchor="w")
        
        # Parameters display
        params_frame = ttk.LabelFrame(left_panel, text="Parameters from First Order", padding="10")
        params_frame.pack(fill="x", pady=(0, 15))
        
        # Power supply specs
        self.ps_label = ttk.Label(params_frame, text="Power Supply: Not available", style='Info.TLabel')
        self.ps_label.pack(anchor="w", pady=2)
        
        self.wire_label = ttk.Label(params_frame, text="Wire: Not available", style='Info.TLabel')
        self.wire_label.pack(anchor="w", pady=2)
        
        self.radius_label = ttk.Label(params_frame, text="Radius: Not available", style='Info.TLabel')
        self.radius_label.pack(anchor="w", pady=2)
        
        # Additional parameters
        additional_frame = ttk.LabelFrame(left_panel, text="Additional Parameters", padding="10")
        additional_frame.pack(fill="x", pady=(0, 15))
        
        # Max turns
        ttk.Label(additional_frame, text="Max Helmholtz Turns:").grid(row=0, column=0, sticky="w", pady=5)
        self.max_turns_var = tk.IntVar(value=100)
        max_turns_entry = ttk.Entry(additional_frame, textvariable=self.max_turns_var, width=15)
        max_turns_entry.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Weight limit
        ttk.Label(additional_frame, text="Weight Limit (lbs):").grid(row=1, column=0, sticky="w", pady=5)
        self.geo_weight_var = tk.DoubleVar(value=140.0)
        geo_weight_entry = ttk.Entry(additional_frame, textvariable=self.geo_weight_var, width=15)
        geo_weight_entry.grid(row=1, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Field point
        ttk.Label(additional_frame, text="Field Point Z (m):").grid(row=2, column=0, sticky="w", pady=5)
        self.field_z_var = tk.DoubleVar(value=0.0)
        field_z_entry = ttk.Entry(additional_frame, textvariable=self.field_z_var, width=15)
        field_z_entry.grid(row=2, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Control buttons
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill="x", pady=(15, 0))
        
        self.geo_optimize_button = ttk.Button(button_frame, text="Run Geometric Optimization", 
                                            command=self.run_geometric_optimization,
                                            style='Primary.TButton')
        self.geo_optimize_button.pack(side="left", padx=(0, 10))
        
        self.clear_geo_button = ttk.Button(button_frame, text="Clear Results", 
                                         command=self.clear_geometric_results,
                                         style='Action.TButton')
        self.clear_geo_button.pack(side="left")
        
        # Right panel - Results
        right_panel = ttk.Frame(content_frame)
        right_panel.pack(side="right", fill="both", expand=True)
        
        # Results text area
        results_frame = ttk.LabelFrame(right_panel, text="Geometric Optimization Results", padding="10")
        results_frame.pack(fill="both", expand=True)
        
        self.geo_results_text = tk.Text(results_frame, height=25, width=60, wrap="word", 
                                       font=("Courier", 9))
        geo_scrollbar = ttk.Scrollbar(results_frame, orient="vertical", command=self.geo_results_text.yview)
        self.geo_results_text.configure(yscrollcommand=geo_scrollbar.set)
        
        self.geo_results_text.pack(side="left", fill="both", expand=True)
        geo_scrollbar.pack(side="right", fill="y")
    
    def run_geometric_optimization(self):
        """Run geometric optimization"""
        if self.first_order_results is None:
            messagebox.showwarning("No Data", "Please run First Order Optimization first")
            return
            
        try:
            # Get parameters from first order results
            ps_specs = self.first_order_results['optimal_power_supply']
            wire_awg = round(self.first_order_results['optimal_awg'])
            radius = self.radius_var.get()
            
            # Get additional parameters
            max_turns = self.max_turns_var.get()
            weight_limit = self.geo_weight_var.get()
            field_z = self.field_z_var.get()
            
            # Update status
            self.geo_results_text.delete(1.0, tk.END)
            self.geo_results_text.insert(tk.END, "Running geometric optimization...\n")
            self.geo_results_text.insert(tk.END, f"Power Supply: {ps_specs['Pmax']:.0f}W/{ps_specs['Imax']:.0f}A/{ps_specs['Vmax']:.0f}V\n")
            self.geo_results_text.insert(tk.END, f"Wire: AWG {wire_awg}\n")
            self.geo_results_text.insert(tk.END, f"Radius: {radius} m\n")
            self.geo_results_text.insert(tk.END, f"Max turns: {max_turns}\n")
            self.geo_results_text.insert(tk.END, f"Weight limit: {weight_limit} lbs\n")
            self.geo_results_text.insert(tk.END, f"Field point: [0, 0, {field_z}]\n\n")
            self.geo_results_text.insert(tk.END, "Optimizing...\n")
            self.root.update()
            
            # Run optimization in separate thread
            thread = threading.Thread(target=self._run_geometric_optimization_thread, 
                                    args=(ps_specs, wire_awg, radius, max_turns, weight_limit, field_z))
            thread.daemon = True
            thread.start()
            
        except Exception as e:
            messagebox.showerror("Error", f"Error starting geometric optimization: {str(e)}")
    
    def _run_geometric_optimization_thread(self, ps_specs, wire_awg, radius, max_turns, weight_limit, field_z):
        """Run geometric optimization in background thread"""
        try:
            # Run optimization
            result = optimize_electromagnet_geometry(
                power_supply_specs=ps_specs,
                wire_awg=wire_awg,
                base_radius=radius,
                max_turns=max_turns,
                weight_limit=weight_limit,
                field_point=[0, 0, field_z],
                max_iterations=100,
                verbose=False
            )
            
            # Update GUI with results
            self.root.after(0, self._display_geometric_results, result)
            
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Optimization Error", str(e)))
    
    def _display_geometric_results(self, result):
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
            self.geo_results_text.insert(tk.END, f"  Axial layers: {geom['axial_layers']}\n")
            self.geo_results_text.insert(tk.END, f"  Total Helmholtz turns: {geom['radial_layers'] * geom['axial_layers']}\n\n")
            
            # Performance
            perf = result['performance']
            self.geo_results_text.insert(tk.END, "PERFORMANCE:\n")
            self.geo_results_text.insert(tk.END, f"  B-field: {perf['b_field_magnitude']:.4f} T\n")
            self.geo_results_text.insert(tk.END, f"  Weight: {perf['weight']:.1f} lbs\n")
            self.geo_results_text.insert(tk.END, f"  Power: {perf['power']:.0f} W\n")
            self.geo_results_text.insert(tk.END, f"  Current: {perf['current']:.1f} A\n")
            self.geo_results_text.insert(tk.END, f"  Voltage: {perf['voltage']:.1f} V\n")
            self.geo_results_text.insert(tk.END, f"  Wire length: {perf['wire_length']:.1f} m\n")
            self.geo_results_text.insert(tk.END, f"  Total turns: {perf['turns_total']}\n")
            self.geo_results_text.insert(tk.END, f"  Helmholtz turns: {perf['helm_turns']}\n\n")
            
            # Optimization info
            info = result['optimization_info']
            self.geo_results_text.insert(tk.END, "OPTIMIZATION INFO:\n")
            self.geo_results_text.insert(tk.END, f"  Global iterations: {info['global_iterations']}\n")
            self.geo_results_text.insert(tk.END, f"  Local iterations: {info['local_iterations']}\n")
            self.geo_results_text.insert(tk.END, f"  Total evaluations: {info['total_evaluations']}\n")
            
            # Enable next stage
            self.notebook.tab(3, state="normal")  # Enable field visualization tab
            
        else:
            self.geo_results_text.insert(tk.END, "GEOMETRIC OPTIMIZATION FAILED\n")
            self.geo_results_text.insert(tk.END, "=" * 50 + "\n\n")
            self.geo_results_text.insert(tk.END, f"Error: {result.get('message', 'Unknown error')}\n")
    
    def clear_geometric_results(self):
        """Clear geometric optimization results"""
        self.geo_results_text.delete(1.0, tk.END)
        self.geometric_results = None
        self.notebook.tab(3, state="disabled")  # Disable field visualization tab
    
    def create_field_visualization_tab(self):
        """Create Field Visualization tab"""
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text="4. Field Visualization")
        self.notebook.tab(3, state="disabled")  # Initially disabled
        
        # Main content frame
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Left panel - Controls
        left_panel = ttk.LabelFrame(content_frame, text="Field Visualization Controls", padding="15")
        left_panel.pack(side="left", fill="y", padx=(0, 10))
        
        # Directions
        directions_frame = ttk.LabelFrame(left_panel, text="Directions", padding="10")
        directions_frame.pack(fill="x", pady=(0, 15))
        
        directions_text = """Field Visualization Stage:

1. Run Geometric Optimization first to get optimal coil geometry

2. This stage visualizes the magnetic field:
   • 2D field magnitude contours
   • Vector field plots
   • Cross-sections through the field

3. Adjust visualization parameters:
   • Grid resolution
   • Field range
   • Plot type

4. Use toolbar to zoom, pan, and save plots"""
        
        directions_label = ttk.Label(directions_frame, text=directions_text, style='Info.TLabel', 
                                   justify="left", wraplength=200)
        directions_label.pack(anchor="w")
        
        # Visualization parameters
        params_frame = ttk.LabelFrame(left_panel, text="Visualization Parameters", padding="10")
        params_frame.pack(fill="x", pady=(0, 15))
        
        # Grid resolution
        ttk.Label(params_frame, text="Grid Resolution:").grid(row=0, column=0, sticky="w", pady=5)
        self.grid_res_var = tk.IntVar(value=20)
        grid_res_entry = ttk.Entry(params_frame, textvariable=self.grid_res_var, width=10)
        grid_res_entry.grid(row=0, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Field range
        ttk.Label(params_frame, text="Field Range (m):").grid(row=1, column=0, sticky="w", pady=5)
        self.field_range_var = tk.DoubleVar(value=0.3)
        field_range_entry = ttk.Entry(params_frame, textvariable=self.field_range_var, width=10)
        field_range_entry.grid(row=1, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Plot type
        ttk.Label(params_frame, text="Plot Type:").grid(row=2, column=0, sticky="w", pady=5)
        self.plot_type_var = tk.StringVar(value="contour")
        plot_type_combo = ttk.Combobox(params_frame, textvariable=self.plot_type_var, 
                                      values=["contour", "vector", "both"], width=10, state="readonly")
        plot_type_combo.grid(row=2, column=1, sticky="w", padx=(10, 0), pady=5)
        
        # Control buttons
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill="x", pady=(15, 0))
        
        self.field_plot_button = ttk.Button(button_frame, text="Generate Field Plot", 
                                          command=self.generate_field_plot,
                                          style='Primary.TButton')
        self.field_plot_button.pack(side="left", padx=(0, 10))
        
        self.save_field_button = ttk.Button(button_frame, text="Save Plot", 
                                          command=self.save_field_plot,
                                          style='Action.TButton')
        self.save_field_button.pack(side="left")
        
        # Right panel - Field plot
        right_panel = ttk.Frame(content_frame)
        right_panel.pack(side="right", fill="both", expand=True)
        
        # Create matplotlib figure
        self.fig_field = Figure(figsize=(10, 8), dpi=100)
        self.canvas_field = FigureCanvasTkAgg(self.fig_field, right_panel)
        self.canvas_field.get_tk_widget().pack(fill="both", expand=True)
        
        # Add navigation toolbar
        toolbar_frame = ttk.Frame(right_panel)
        toolbar_frame.pack(fill="x")
        self.toolbar_field = NavigationToolbar2Tk(self.canvas_field, toolbar_frame)
        self.toolbar_field.update()
    
    def generate_field_plot(self):
        """Generate magnetic field visualization"""
        if self.geometric_results is None:
            messagebox.showwarning("No Data", "Please run Geometric Optimization first")
            return
            
        try:
            # Clear previous plot
            self.fig_field.clear()
            
            # Get optimal geometry
            geom = self.geometric_results['optimal_geometry']
            perf = self.geometric_results['performance']
            
            # Get parameters
            grid_res = self.grid_res_var.get()
            field_range = self.field_range_var.get()
            plot_type = self.plot_type_var.get()
            
            # Create power supply and wire objects
            ps_specs = self.first_order_results['optimal_power_supply']
            wire_awg = round(self.first_order_results['optimal_awg'])
            radius = self.radius_var.get()
            
            power_supply = PowerSupply("Optimized", ps_specs['max_voltage'], 
                                     ps_specs['max_current'], ps_specs['max_power'])
            wire = Wire(wire_awg)
            
            # Create electromagnet with optimal geometry
            emag = Emag(power_supply, wire)
            
            # Build coil geometry
            for layer_idx in range(geom['radial_layers']):
                coil_radius = radius + (layer_idx * wire.diameter_nom_m)
                
                for turn_idx in range(geom['axial_layers']):
                    separation = radius + (2 * turn_idx * wire.diameter_nom_m)
                    
                    emag.add_mirrored_turns(
                        radius=coil_radius,
                        separation=separation,
                        center=np.array([0, 0, 0]),
                        points_per_turn=100
                    )
            
            emag.update_parameters()
            
            # Define field grid
            x_range = (-field_range, field_range)
            y_range = (-field_range, field_range)
            z_range = (0, 0)  # XZ plane
            
            emag.def_field(x_range, y_range, z_range, points_per_axis=grid_res)
            emag.calc_b_field()
            
            # Reshape for plotting
            X = emag.field[:, 0].reshape(grid_res, grid_res)
            Z = emag.field[:, 2].reshape(grid_res, grid_res)
            B_magnitude = np.linalg.norm(emag.b_field, axis=1).reshape(grid_res, grid_res)
            Bx = emag.b_field[:, 0].reshape(grid_res, grid_res)
            Bz = emag.b_field[:, 2].reshape(grid_res, grid_res)
            
            # Create plot
            if plot_type in ["contour", "both"]:
                ax1 = self.fig_field.add_subplot(111)
                contour = ax1.contourf(X, Z, B_magnitude, levels=20, cmap='viridis')
                ax1.set_xlabel('X (m)')
                ax1.set_ylabel('Z (m)')
                ax1.set_title('Magnetic Field Magnitude (T)')
                ax1.set_aspect('equal')
                ax1.grid(True, alpha=0.3)
                self.fig_field.colorbar(contour, ax=ax1, label='Field Strength (T)')
                
                if plot_type == "both":
                    # Add vector field overlay
                    step = max(1, grid_res // 15)
                    ax1.quiver(X[::step, ::step], Z[::step, ::step], 
                              Bx[::step, ::step], Bz[::step, ::step],
                              scale=50, color='white', alpha=0.7)
            
            elif plot_type == "vector":
                ax1 = self.fig_field.add_subplot(111)
                step = max(1, grid_res // 15)
                ax1.quiver(X[::step, ::step], Z[::step, ::step], 
                          Bx[::step, ::step], Bz[::step, ::step],
                          scale=50, color='blue', alpha=0.7)
                ax1.set_xlabel('X (m)')
                ax1.set_ylabel('Z (m)')
                ax1.set_title('Magnetic Field Vectors')
                ax1.set_aspect('equal')
                ax1.grid(True, alpha=0.3)
            
            # Add coil geometry
            # Plot coil positions (simplified)
            coil_centers = []
            for layer_idx in range(geom['radial_layers']):
                coil_radius = radius + (layer_idx * wire.diameter_nom_m)
                for turn_idx in range(geom['axial_layers']):
                    separation = radius + (2 * turn_idx * wire.diameter_nom_m)
                    coil_centers.append((0, separation/2))
                    coil_centers.append((0, -separation/2))
            
            for center in coil_centers:
                circle = plt.Circle(center, radius, fill=False, color='red', linewidth=2, alpha=0.7)
                ax1.add_patch(circle)
            
            # Refresh canvas
            self.canvas_field.draw()
            
        except Exception as e:
            messagebox.showerror("Plot Error", f"Error generating field plot: {str(e)}")
    
    def save_field_plot(self):
        """Save field visualization plot"""
        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("PNG files", "*.png"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                self.fig_field.savefig(filename, bbox_inches='tight', dpi=300)
                messagebox.showinfo("Success", f"Field plot saved to {filename}")
                
            except Exception as e:
                messagebox.showerror("Save Error", f"Error saving plot: {str(e)}")
    
    def create_design_space_visualization_tab(self):
        """Create Design Space Visualization tab"""
        tab_frame = ttk.Frame(self.notebook)
        self.notebook.add(tab_frame, text="5. Design Space Visualization")
        self.notebook.tab(4, state="disabled")  # Initially disabled
        
        # Main content frame
        content_frame = ttk.Frame(tab_frame)
        content_frame.pack(fill="both", expand=True, padx=10, pady=10)
        
        # Left panel - Controls
        left_panel = ttk.LabelFrame(content_frame, text="Design Space Controls", padding="15")
        left_panel.pack(side="left", fill="y", padx=(0, 10))
        
        # Directions
        directions_frame = ttk.LabelFrame(left_panel, text="Directions", padding="10")
        directions_frame.pack(fill="x", pady=(0, 15))
        
        directions_text = """Design Space Visualization Stage:

1. Run Geometric Optimization first to get power supply and wire specs

2. This stage scans the design space:
   • B-field vs geometry parameters
   • Weight vs geometry parameters
   • Power vs geometry parameters
   • Valid design regions

3. Adjust scan parameters:
   • Radial layer range
   • Axial layer range
   • Scan resolution

4. Visualize trade-offs and constraints"""
        
        directions_label = ttk.Label(directions_frame, text=directions_text, style='Info.TLabel', 
                                   justify="left", wraplength=200)
        directions_label.pack(anchor="w")
        
        # Scan parameters
        params_frame = ttk.LabelFrame(left_panel, text="Scan Parameters", padding="10")
        params_frame.pack(fill="x", pady=(0, 15))
        
        # Radial range
        ttk.Label(params_frame, text="Radial Range:").grid(row=0, column=0, sticky="w", pady=5)
        self.radial_min_var = tk.IntVar(value=1)
        self.radial_max_var = tk.IntVar(value=10)
        ttk.Entry(params_frame, textvariable=self.radial_min_var, width=5).grid(row=0, column=1, sticky="w", padx=(5, 0), pady=5)
        ttk.Label(params_frame, text="to").grid(row=0, column=2, sticky="w", padx=5, pady=5)
        ttk.Entry(params_frame, textvariable=self.radial_max_var, width=5).grid(row=0, column=3, sticky="w", padx=(5, 0), pady=5)
        
        # Axial range
        ttk.Label(params_frame, text="Axial Range:").grid(row=1, column=0, sticky="w", pady=5)
        self.axial_min_var = tk.IntVar(value=1)
        self.axial_max_var = tk.IntVar(value=10)
        ttk.Entry(params_frame, textvariable=self.axial_min_var, width=5).grid(row=1, column=1, sticky="w", padx=(5, 0), pady=5)
        ttk.Label(params_frame, text="to").grid(row=1, column=2, sticky="w", padx=5, pady=5)
        ttk.Entry(params_frame, textvariable=self.axial_max_var, width=5).grid(row=1, column=3, sticky="w", padx=(5, 0), pady=5)
        
        # Control buttons
        button_frame = ttk.Frame(left_panel)
        button_frame.pack(fill="x", pady=(15, 0))
        
        self.scan_button = ttk.Button(button_frame, text="Scan Design Space", 
                                    command=self.scan_design_space,
                                    style='Primary.TButton')
        self.scan_button.pack(side="left", padx=(0, 10))
        
        self.save_scan_button = ttk.Button(button_frame, text="Save Plots", 
                                         command=self.save_design_space_plots,
                                         style='Action.TButton')
        self.save_scan_button.pack(side="left")
        
        # Right panel - Design space plots
        right_panel = ttk.Frame(content_frame)
        right_panel.pack(side="right", fill="both", expand=True)
        
        # Create matplotlib figure
        self.fig_design_space = Figure(figsize=(12, 10), dpi=100)
        self.canvas_design_space = FigureCanvasTkAgg(self.fig_design_space, right_panel)
        self.canvas_design_space.get_tk_widget().pack(fill="both", expand=True)
        
        # Add navigation toolbar
        toolbar_frame = ttk.Frame(right_panel)
        toolbar_frame.pack(fill="x")
        self.toolbar_design_space = NavigationToolbar2Tk(self.canvas_design_space, toolbar_frame)
        self.toolbar_design_space.update()
    
    def scan_design_space(self):
        """Scan and visualize design space"""
        if self.first_order_results is None:
            messagebox.showwarning("No Data", "Please run First Order Optimization first")
            return
            
        try:
            # Get parameters
            ps_specs = self.first_order_results['optimal_power_supply']
            wire_awg = round(self.first_order_results['optimal_awg'])
            radius = self.radius_var.get()
            weight_limit = self.geo_weight_var.get()
            max_turns = self.max_turns_var.get()
            
            # Get scan ranges
            radial_range = (self.radial_min_var.get(), self.radial_max_var.get())
            axial_range = (self.axial_min_var.get(), self.axial_max_var.get())
            
            # Create power supply and wire objects
            power_supply = PowerSupply("Optimized", ps_specs['max_voltage'], 
                                     ps_specs['max_current'], ps_specs['max_power'])
            wire = Wire(wire_awg)
            
            # Create optimizer
            optimizer = GeometryOptimizer(
                power_supply=power_supply,
                wire=wire,
                base_radius=radius,
                max_turns=max_turns,
                weight_limit=weight_limit,
                field_point=np.array([0, 0, 0])
            )
            
            # Run scan in separate thread
            thread = threading.Thread(target=self._run_design_space_scan_thread, 
                                    args=(optimizer, radial_range, axial_range))
            thread.daemon = True
            thread.start()
            
        except Exception as e:
            messagebox.showerror("Error", f"Error starting design space scan: {str(e)}")
    
    def _run_design_space_scan_thread(self, optimizer, radial_range, axial_range):
        """Run design space scan in background thread"""
        try:
            # Run scan
            scan_result = optimizer.scan_design_space(
                radial_range=radial_range,
                axial_range=axial_range,
                plot=False
            )
            
            # Update GUI with results
            self.root.after(0, self._display_design_space_results, scan_result)
            
        except Exception as e:
            self.root.after(0, lambda: messagebox.showerror("Scan Error", str(e)))
    
    def _display_design_space_results(self, scan_result):
        """Display design space scan results"""
        try:
            # Clear previous plots
            self.fig_design_space.clear()
            
            # Create 2x2 subplot layout
            gs = self.fig_design_space.add_gridspec(2, 2, hspace=0.3, wspace=0.3)
            
            radial_vals = scan_result['radial_values']
            axial_vals = scan_result['axial_values']
            X, Y = np.meshgrid(axial_vals, radial_vals)
            
            # B-field
            ax1 = self.fig_design_space.add_subplot(gs[0, 0])
            im1 = ax1.contourf(X, Y, scan_result['b_field_grid'], levels=20, cmap='viridis')
            ax1.set_xlabel('Axial Layers')
            ax1.set_ylabel('Radial Layers')
            ax1.set_title('B-field (T)')
            self.fig_design_space.colorbar(im1, ax=ax1)
            
            # Mark optimal
            opt = scan_result['optimal_point']
            ax1.plot(opt['axial_layers'], opt['radial_layers'], 'r*', markersize=15)
            
            # Max turns constraint line
            axial_line = np.linspace(1, max(axial_vals), 100)
            radial_line = scan_result.get('max_turns', 100) / axial_line
            mask = (radial_line >= 1) & (radial_line <= max(radial_vals))
            ax1.plot(axial_line[mask], radial_line[mask], 'w--', linewidth=2)
            
            # Weight
            ax2 = self.fig_design_space.add_subplot(gs[0, 1])
            im2 = ax2.contourf(X, Y, scan_result['weight_grid'], levels=20, cmap='plasma')
            ax2.set_xlabel('Axial Layers')
            ax2.set_ylabel('Radial Layers')
            ax2.set_title('Weight (lbs)')
            self.fig_design_space.colorbar(im2, ax=ax2)
            ax2.plot(axial_line[mask], radial_line[mask], 'w--', linewidth=2)
            
            # Power
            ax3 = self.fig_design_space.add_subplot(gs[1, 0])
            im3 = ax3.contourf(X, Y, scan_result['power_grid'], levels=20, cmap='inferno')
            ax3.set_xlabel('Axial Layers')
            ax3.set_ylabel('Radial Layers')
            ax3.set_title('Power (W)')
            self.fig_design_space.colorbar(im3, ax=ax3)
            ax3.plot(axial_line[mask], radial_line[mask], 'w--', linewidth=2)
            
            # Valid designs
            ax4 = self.fig_design_space.add_subplot(gs[1, 1])
            ax4.contourf(X, Y, scan_result['valid_grid'].astype(int), 
                        levels=[0, 0.5, 1], colors=['red', 'green'], alpha=0.7)
            ax4.set_xlabel('Axial Layers')
            ax4.set_ylabel('Radial Layers')
            ax4.set_title('Valid Designs')
            ax4.plot(axial_line[mask], radial_line[mask], 'k--', linewidth=2)
            
            # Set main title
            self.fig_design_space.suptitle(f'Design Space (Max Turns = {scan_result.get("max_turns", "N/A")})', 
                                         fontsize=14, fontweight='bold')
            
            # Store results
            self.design_space_results = scan_result
            
            # Refresh canvas
            self.canvas_design_space.draw()
            
        except Exception as e:
            messagebox.showerror("Display Error", f"Error displaying design space results: {str(e)}")
    
    def save_design_space_plots(self):
        """Save design space plots"""
        if self.design_space_results is None:
            messagebox.showwarning("No Data", "No design space plots to save")
            return
            
        filename = filedialog.asksaveasfilename(
            defaultextension=".pdf",
            filetypes=[("PDF files", "*.pdf"), ("PNG files", "*.png"), ("All files", "*.*")]
        )
        
        if filename:
            try:
                self.fig_design_space.savefig(filename, bbox_inches='tight', dpi=300)
                messagebox.showinfo("Success", f"Design space plots saved to {filename}")
                
            except Exception as e:
                messagebox.showerror("Save Error", f"Error saving plots: {str(e)}")


def main():
    """Main function to run the GUI"""
    root = tk.Tk()
    app = ModernHelmholtzGUI(root)
    root.mainloop()


if __name__ == "__main__":
    main()
