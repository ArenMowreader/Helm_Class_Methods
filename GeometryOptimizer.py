'''
Electromagnet Geometry Optimizer

Optimizes coil geometry (radial layers and axial turns) for maximum magnetic field
using hybrid optimization (global search + local refinement).

Design space: All combinations where radial_layers * axial_layers <= max_turns

Author: Based on electromagnet.py by Aren Mowreader
Date: 09/11/2025
Oregon State University
'''

import numpy as np
from scipy.optimize import minimize, differential_evolution
import matplotlib.pyplot as plt
from electromagnet import PowerSupply, Wire, Emag

class GeometryOptimizer:
    """
    Optimize electromagnet geometry for maximum B-field using hybrid optimization.
    
    Parameters
    ----------
    power_supply : PowerSupply
        Power supply specifications
    wire : Wire
        Wire specifications (AWG)
    base_radius : float
        Inner coil radius (meters)
    weight_limit : float
        Maximum weight limit for a single coil (pounds)
        Note: The electromagnet.weight will return weight for both coils
    field_point : ndarray, optional
        Target field point. Default [0, 0, 0]
    """
    
    def __init__(self, power_supply, wire, base_radius, weight_limit, field_point=None):
        self.power_supply = power_supply
        self.wire = wire
        self.base_radius = base_radius
        self.weight_limit = weight_limit  # Single coil weight limit
        self.field_point = field_point if field_point is not None else np.array([0, 0, 0])
        
        # Calculate max_turns based on weight limit
        # For a single coil: weight = 2 * pi * radius * turns * lbs_per_meter
        # Solving for turns: turns = weight / (2 * pi * radius * lbs_per_meter)
        # We use a conservative estimate with average radius
        avg_radius = base_radius + (5 * wire.diameter_nom_m)  # Assume ~5 layers average
        self.max_turns = int(weight_limit / (2 * np.pi * avg_radius * wire.lbs_per_meter))
        
        # Set search bounds based on calculated max_turns
        self.max_search_dim = int(np.sqrt(self.max_turns)) + 10
        
        print(f"Optimizer initialized:")
        print(f"  Max turns (calculated): {self.max_turns}")
        print(f"  Weight limit (single coil): {weight_limit} lbs")
        print(f"  Wire: AWG {wire.AWG}")
        print(f"  Base radius: {base_radius*1000:.1f} mm")
    
    def is_valid_design(self, radial_layers, axial_layers):
        """Check max_turns constraint."""
        return radial_layers * axial_layers <= self.max_turns
    
    def create_coil_geometry(self, radial_layers, axial_layers):
        """Create coil geometry matching scale_run.py pattern."""
        if not self.is_valid_design(radial_layers, axial_layers):
            raise ValueError(f"Design violates max_turns constraint")
        
        emag = Emag(self.power_supply, self.wire)
        
        for layer_idx in range(radial_layers):
            radius = self.base_radius + (layer_idx * self.wire.diameter_nom_m)
            
            for turn_idx in range(axial_layers):
                separation = self.base_radius + (2 * turn_idx * self.wire.diameter_nom_m)
                
                emag.add_mirrored_turns(
                    radius=radius,
                    separation=separation,
                    center=np.array([0, 0, 0]),
                    points_per_turn=200
                )
        
        return emag
    
    def evaluate_design(self, params):
        """Evaluate design performance."""
        radial_layers = max(1, int(round(params[0])))
        axial_layers = max(1, int(round(params[1])))
        
        # Check constraints
        if not self.is_valid_design(radial_layers, axial_layers):
            return {'b_field_magnitude': 0, 'valid': False, 'constraint': 'max_turns'}
        
        try:
            emag = self.create_coil_geometry(radial_layers, axial_layers)
            emag.update_parameters()
            
            # emag.weight is for both coils, but weight_limit is for single coil
            if emag.weight > (2 * self.weight_limit):
                return {'b_field_magnitude': 0, 'valid': False, 'constraint': 'weight', 'weight': emag.weight}
            
            # Calculate B-field
            emag.def_field(
                x_range=(self.field_point[0], self.field_point[0]),
                y_range=(self.field_point[1], self.field_point[1]), 
                z_range=(self.field_point[2], self.field_point[2]),
                points_per_axis=1
            )
            emag.calc_b_field()
            
            return {
                'b_field_magnitude': np.linalg.norm(emag.b_field[0]),
                'weight': emag.weight,
                'power': emag.power,
                'current': emag.current,
                'voltage': emag.voltage,
                'resistance': emag.resistance,
                'turns_total': emag.net_turns,
                'helm_turns': emag.helm_turns,
                'wire_length': emag.coil_length,
                'radial_layers': radial_layers,
                'axial_layers': axial_layers,
                'valid': True,
                'emag': emag
            }
            
        except Exception as e:
            return {'b_field_magnitude': 0, 'valid': False, 'constraint': f'error: {str(e)}'}
    
    def objective_function(self, params):
        """Objective for scipy (negative B-field for minimization)."""
        result = self.evaluate_design(params)
        return 1e6 if not result['valid'] else -result['b_field_magnitude']
    
    def optimize(self, max_iterations=100, verbose=True):
        """
        Run hybrid optimization (global + local).
        
        Parameters
        ----------
        max_iterations : int
            Maximum total iterations
        verbose : bool
            Print progress
            
        Returns
        -------
        dict
            Optimization results
        """
        bounds = [(1, self.max_search_dim), (1, self.max_search_dim)]
        
        if verbose:
            print(f"\nRunning hybrid optimization...")
            print(f"Search bounds: {self.max_search_dim} x {self.max_search_dim}")
        
        # Global search
        if verbose:
            print("Phase 1: Global search...")
        
        global_result = differential_evolution(
            func=self.objective_function,
            bounds=bounds,
            maxiter=max_iterations//2,
            popsize=15,
            seed=42,
            disp=verbose
        )
        
        # Local refinement
        if global_result.success:
            if verbose:
                print("Phase 2: Local refinement...")
            
            result = minimize(
                fun=self.objective_function,
                x0=global_result.x,
                method='L-BFGS-B',
                bounds=bounds,
                options={'maxiter': max_iterations//2, 'disp': verbose}
            )
        else:
            if verbose:
                print("Global search failed, using global result")
            result = global_result
        
        # Process results
        if result.success:
            optimal_radial = max(1, int(round(result.x[0])))
            optimal_axial = max(1, int(round(result.x[1])))
            performance = self.evaluate_design([optimal_radial, optimal_axial])
            
            optimization_result = {
                'success': True,
                'optimal_geometry': {
                    'radial_layers': optimal_radial,
                    'axial_layers': optimal_axial
                },
                'performance': performance,
                'optimization_info': {
                    'global_iterations': getattr(global_result, 'nit', 'N/A'),
                    'local_iterations': getattr(result, 'nit', 'N/A'),
                    'total_evaluations': getattr(result, 'nfev', 'N/A')
                }
            }
            
            if verbose:
                self.print_results(optimization_result)
            
            return optimization_result
        else:
            return {'success': False, 'message': result.message}
    
    def scan_design_space(self, radial_range=(1, 15), axial_range=(1, 15), plot=True):
        """
        Scan and visualize design space.
        
        Parameters
        ----------
        radial_range : tuple
            Radial layer range to scan
        axial_range : tuple  
            Axial layer range to scan
        plot : bool
            Create visualization
            
        Returns
        -------
        dict
            Scan results with performance grids
        """
        radial_values = np.arange(radial_range[0], radial_range[1] + 1)
        axial_values = np.arange(axial_range[0], axial_range[1] + 1)
        
        print(f"Scanning {len(radial_values)} x {len(axial_values)} design space...")
        
        # Initialize grids
        b_field_grid = np.zeros((len(radial_values), len(axial_values)))
        weight_grid = np.zeros((len(radial_values), len(axial_values)))
        power_grid = np.zeros((len(radial_values), len(axial_values)))
        valid_grid = np.zeros((len(radial_values), len(axial_values)), dtype=bool)
        
        # Evaluate each point
        for i, radial in enumerate(radial_values):
            for j, axial in enumerate(axial_values):
                if self.is_valid_design(radial, axial):
                    result = self.evaluate_design([radial, axial])
                    if result['valid']:
                        b_field_grid[i, j] = result['b_field_magnitude']
                        weight_grid[i, j] = result['weight']
                        power_grid[i, j] = result['power']
                        valid_grid[i, j] = True
        
        # Find optimal
        valid_b_field = np.where(valid_grid, b_field_grid, 0)
        max_idx = np.unravel_index(np.argmax(valid_b_field), valid_b_field.shape)
        
        scan_result = {
            'radial_values': radial_values,
            'axial_values': axial_values,
            'b_field_grid': b_field_grid,
            'weight_grid': weight_grid,
            'power_grid': power_grid,
            'valid_grid': valid_grid,
            'optimal_point': {
                'radial_layers': radial_values[max_idx[0]],
                'axial_layers': axial_values[max_idx[1]],
                'b_field_tesla': valid_b_field[max_idx]
            }
        }
        
        if plot:
            self._plot_scan(scan_result)
        
        return scan_result
    
    def _plot_scan(self, scan_result):
        """Create scan visualization."""
        fig, axes = plt.subplots(2, 2, figsize=(12, 10))
        
        radial_vals = scan_result['radial_values']
        axial_vals = scan_result['axial_values']
        X, Y = np.meshgrid(axial_vals, radial_vals)
        
        # B-field
        im1 = axes[0,0].contourf(X, Y, scan_result['b_field_grid'], levels=20, cmap='viridis')
        axes[0,0].set_xlabel('Axial Layers')
        axes[0,0].set_ylabel('Radial Layers')
        axes[0,0].set_title('B-field (T)')
        plt.colorbar(im1, ax=axes[0,0])
        
        # Mark optimal
        opt = scan_result['optimal_point']
        axes[0,0].plot(opt['axial_layers'], opt['radial_layers'], 'r*', markersize=15)
        
        # Max turns constraint line
        axial_line = np.linspace(1, max(axial_vals), 100)
        radial_line = self.max_turns / axial_line
        mask = (radial_line >= 1) & (radial_line <= max(radial_vals))
        axes[0,0].plot(axial_line[mask], radial_line[mask], 'w--', linewidth=2)
        
        # Weight
        im2 = axes[0,1].contourf(X, Y, scan_result['weight_grid'], levels=20, cmap='plasma')
        axes[0,1].set_xlabel('Axial Layers')
        axes[0,1].set_ylabel('Radial Layers')
        axes[0,1].set_title('Weight (lbs)')
        plt.colorbar(im2, ax=axes[0,1])
        axes[0,1].plot(axial_line[mask], radial_line[mask], 'w--', linewidth=2)
        
        # Power
        im3 = axes[1,0].contourf(X, Y, scan_result['power_grid'], levels=20, cmap='inferno')
        axes[1,0].set_xlabel('Axial Layers')
        axes[1,0].set_ylabel('Radial Layers')
        axes[1,0].set_title('Power (W)')
        plt.colorbar(im3, ax=axes[1,0])
        axes[1,0].plot(axial_line[mask], radial_line[mask], 'w--', linewidth=2)
        
        # Valid designs
        axes[1,1].contourf(X, Y, scan_result['valid_grid'].astype(int), 
                          levels=[0, 0.5, 1], colors=['red', 'green'], alpha=0.7)
        axes[1,1].set_xlabel('Axial Layers')
        axes[1,1].set_ylabel('Radial Layers')
        axes[1,1].set_title('Valid Designs')
        axes[1,1].plot(axial_line[mask], radial_line[mask], 'k--', linewidth=2)
        
        plt.suptitle(f'Design Space (Max Turns = {self.max_turns})')
        plt.tight_layout()
        plt.show()
    
    def print_results(self, result):
        """Print formatted results."""
        if not result['success']:
            print(f"Optimization failed: {result['message']}")
            return
        
        print("\n" + "="*50)
        print("OPTIMAL ELECTROMAGNET GEOMETRY")
        print("="*50)
        
        geom = result['optimal_geometry']
        perf = result['performance']
        info = result['optimization_info']
        
        print(f"Geometry:")
        print(f"  Radial layers: {geom['radial_layers']}")
        
        # Calculate geometry dimensions
        radial_depth = (geom['radial_layers'] - 1) * self.wire.diameter_nom_m
        outside_radius = self.base_radius + radial_depth
        axial_height = geom['axial_layers'] * self.wire.diameter_nom_m
        
        print(f"  Radial depth: {radial_depth*1000:.1f} mm")
        print(f"  Outside radius: {outside_radius*1000:.1f} mm")
        print(f"  Axial layers: {geom['axial_layers']}")
        print(f"  Axial height: {axial_height*1000:.1f} mm")
        print(f"  Single coil turns: {perf['helm_turns']}/{self.max_turns}")
        
        print(f"\nPerformance:")
        print(f"  B-field: {perf['b_field_magnitude']:.4f} T")
        print(f"  Weight (both coils): {perf['weight']:.1f} lbs")
        print(f"  Weight (single coil): {perf['weight']/2:.1f}/{self.weight_limit} lbs")
        print(f"  Power: {perf['power']:.0f} W")
        print(f"  Current: {perf['current']:.1f} A")
        print(f"  Voltage: {perf['voltage']:.1f} V")
        print(f"  Wire length: {perf['wire_length']:.1f} m")
        
        print(f"\nOptimization:")
        print(f"  Global iterations: {info['global_iterations']}")
        print(f"  Local iterations: {info['local_iterations']}")
        print(f"  Total evaluations: {info['total_evaluations']}")
        
        print("="*50)


def optimize_electromagnet_geometry(power_supply_specs, wire_awg, base_radius, 
                                  weight_limit, field_point=None, 
                                  max_iterations=100, verbose=True):
    """
    Optimize electromagnet geometry for maximum B-field.
    
    Parameters
    ----------
    power_supply_specs : dict
        {'max_voltage': V, 'max_current': A, 'max_power': W}
    wire_awg : int
        Wire gauge
    base_radius : float
        Inner coil radius (meters)
    weight_limit : float
        Maximum weight for a single coil (pounds)
    field_point : array-like, optional
        Target point [x,y,z]. Default [0,0,0]
    max_iterations : int
        Maximum optimization iterations
    verbose : bool
        Print progress and results
        
    Returns
    -------
    dict
        Optimization results with optimal geometry and performance
        
    Examples
    --------
    # Basic optimization
    result = optimize_electromagnet_geometry(
        power_supply_specs={'max_voltage': 500, 'max_current': 40, 'max_power': 10000},
        wire_awg=12,
        base_radius=0.125,
        weight_limit=140
    )
    
    # Optimize for off-center field point
    result = optimize_electromagnet_geometry(
        power_supply_specs={'max_voltage': 300, 'max_current': 50, 'max_power': 10000},
        wire_awg=14,
        base_radius=0.1,
        weight_limit=75,
        field_point=[0, 0, 0.01]
    )
    """
    power_supply = PowerSupply(
        name="Optimizer",
        max_voltage=power_supply_specs['max_voltage'],
        max_current=power_supply_specs['max_current'],
        max_power=power_supply_specs['max_power']
    )
    
    wire = Wire(wire_awg)
    
    optimizer = GeometryOptimizer(
        power_supply=power_supply,
        wire=wire,
        base_radius=base_radius,
        weight_limit=weight_limit,
        field_point=np.array(field_point) if field_point else None
    )
    
    return optimizer.optimize(max_iterations=max_iterations, verbose=verbose)


# Example usage
if __name__ == "__main__":
    
    # Example 1: Basic optimization
    print("Example 1: Basic optimization")
    result = optimize_electromagnet_geometry(
        power_supply_specs={'max_voltage': 500, 'max_current': 40, 'max_power': 10000},
        wire_awg=12,
        base_radius=0.125,
        max_turns=100,
        weight_limit=140
    )
    
    if result['success']:
        geom = result['optimal_geometry']
        perf = result['performance']
        print(f"\nResult: {geom['radial_layers']} x {geom['axial_layers']} layers")
        print(f"B-field: {perf['b_field_magnitude']:.4f} T")
    
    # Example 2: Design space scan
    print("\n" + "="*60)
    print("Example 2: Design space visualization")
    
    power_supply = PowerSupply("Test", 500, 40, 10000)
    wire = Wire(12)
    
    optimizer = GeometryOptimizer(
        power_supply=power_supply,
        wire=wire,
        base_radius=0.125,
        max_turns=50,  # Smaller for faster scan
        weight_limit=140
    )
    
    scan_results = optimizer.scan_design_space(
        radial_range=(1, 10),
        axial_range=(1, 10),
        plot=True
    )
    
    opt = scan_results['optimal_point']
    print(f"Scan optimal: {opt['radial_layers']} x {opt['axial_layers']} = {opt['b_field_tesla']:.4f} T")