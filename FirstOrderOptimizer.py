'''
Author: Aren Mowreader
Date: 9/11/2025
Description:
    Several methods to optimize the design of a Helmholtz coil.
    -First Order Method to find optimal power supply and wire gauge.
    -Line Integral Method to find optimal coil geometry with given power supply
     and wire gauge.
    Uses scipy.optimize.minimize to find the optimal design.
    Optimization Assumes programable power supply that works in Constant Current (CC)
    and Constant Voltage (CV) mode.
'''
import numpy as np
import pandas as pd
from scipy.optimize import minimize, differential_evolution
from scipy.interpolate import interp1d
import matplotlib.pyplot as plt
from helmholtz_utils import b_helm_current, weight, analyze_power_supply

class FirstOrderOptimizer:
    """
    Optimization using realistic power supply characteristics with CC/CV/CP modes.
    Uses hybrid global+local optimization to avoid local optima.
    """
    
    def __init__(self, radius, weight_limit, max_power=50000, wire_data_file='MWS_wire_data.csv'):
        self.radius = radius
        self.weight_limit = weight_limit  # Weight limit for SINGLE coil
        self.max_power = max_power  # Maximum power budget
        
        # Load and interpolate wire data
        self.wire_data = pd.read_csv(wire_data_file)
        self._setup_wire_interpolations()
        
    def _setup_wire_interpolations(self):
        """Create interpolation functions for wire properties vs AWG."""
        awg_sizes = self.wire_data['SIZE_AWG'].values
        ohms_per_meter = self.wire_data['ohms_per_meter'].values
        lbs_per_meter = self.wire_data['lbs_per_meter'].values
        diameter_m = self.wire_data['DIAMETER_NOM'].values * 0.0254
        
        self.ohms_interp = interp1d(awg_sizes, ohms_per_meter, 
                                   kind='linear', bounds_error=False, fill_value='extrapolate')
        self.lbs_interp = interp1d(awg_sizes, lbs_per_meter, 
                                  kind='linear', bounds_error=False, fill_value='extrapolate')
        self.diameter_interp = interp1d(awg_sizes, diameter_m, 
                                       kind='linear', bounds_error=False, fill_value='extrapolate')
        
        self.awg_min = awg_sizes.min()
        self.awg_max = awg_sizes.max()
    
    def get_wire_properties(self, awg_continuous):
        """Get interpolated wire properties for continuous AWG value."""
        return {
            'ohms_per_meter': float(self.ohms_interp(awg_continuous)),
            'lbs_per_meter': float(self.lbs_interp(awg_continuous)),
            'diameter_m': float(self.diameter_interp(awg_continuous))
        }
    
    def analyze_power_supply_performance(self, power_supply_specs, awg):
        """
        Analyze power supply performance across range of turns using realistic CC/CV/CP behavior.
        """
        wire_props = self.get_wire_properties(awg)
        
        # Calculate maximum turns from weight constraint (SINGLE coil weight limit)
        circumference = 2 * np.pi * self.radius
        max_turns_weight = self.weight_limit / (wire_props['lbs_per_meter'] * circumference)
        
        max_turns = int(max_turns_weight)
        
        if max_turns < 1:
            return None
        
        # Create array of turns to test (single coil)
        n_single_coil = np.arange(1, max_turns + 1)
        
        # Use the helmholtz_utils analyze_power_supply function
        ps_analysis = analyze_power_supply(
            power_supply_specs,
            wire_props['ohms_per_meter'], 
            n_single_coil,
            self.radius
        )
        
        # Calculate weights for each configuration (SINGLE coil only) - vectorized
        weights = n_single_coil * circumference * wire_props['lbs_per_meter']
        
        # Filter to only configurations under weight limit (single coil)
        valid_mask = weights <= self.weight_limit
        
        if not np.any(valid_mask):
            return None
        
        return {
            'turns_single': n_single_coil[valid_mask],
            'b_field': ps_analysis['B_current'][valid_mask],
            'current': ps_analysis['I'][valid_mask],
            'voltage': ps_analysis['V'][valid_mask], 
            'power': ps_analysis['power'][valid_mask],
            'weights': weights[valid_mask],  # Single coil weights
            'wire_length': n_single_coil[valid_mask] * circumference * 2,  # Both coils
            'resistance': ps_analysis['V'][valid_mask] / ps_analysis['I'][valid_mask]
        }
    
    def objective_function(self, params):
        """
        Objective function that finds maximum B-field and validates constraints.
        """
        Pmax, Imax, Vmax, awg = params
        
        try:
            power_supply_specs = {'Pmax': Pmax, 'Imax': Imax, 'Vmax': Vmax}
            
            # Get power supply performance analysis (single calculation)
            performance = self.analyze_power_supply_performance(power_supply_specs, awg)
            
            if performance is None or len(performance['b_field']) == 0:
                return 1e6  # No valid configurations
            
            # Find maximum B-field point
            max_idx = np.argmax(performance['b_field'])
            max_b_field = performance['b_field'][max_idx]
            
            # Check constraints at the optimal operating point
            operating_voltage = performance['voltage'][max_idx]
            operating_current = performance['current'][max_idx]
            operating_power = performance['power'][max_idx]
            operating_weight_single_coil = performance['weights'][max_idx]  # Single coil weight
            
            # Weight constraint check (single coil)
            if operating_weight_single_coil > self.weight_limit:
                return 1e6 + (operating_weight_single_coil - self.weight_limit) * 1000
            
            # Power budget constraint check
            if operating_power > self.max_power:
                return 1e6 + (operating_power - self.max_power) * 1000
            
            # Power supply capability constraints
            voltage_violation = max(0, operating_voltage - Vmax)
            current_violation = max(0, operating_current - Imax)
            power_violation = max(0, operating_power - Pmax)
            
            total_violation = voltage_violation + current_violation + power_violation
            
            if total_violation > 0:
                return 1e6 + total_violation * 1000
            
            # All constraints satisfied, return negative B-field for minimization
            return -max_b_field
            
        except Exception as e:
            return 1e6
    
    def optimize_electromagnet(self, power_range=None, current_range=None, 
                              voltage_range=None, awg_range=None,
                              fixed_power=None, fixed_current=None, fixed_voltage=None, fixed_awg=None):
        """
        Optimize power supply specifications and wire gauge for maximum B-field.
        Uses hybrid global+local optimization to avoid local optima.
        """
        if power_range is None:
            power_range = (0, 10000)
            
        if current_range is None:
            current_range = (0, 1000)
            
        if voltage_range is None:
            voltage_range = (0, self.max_power)
            
        if awg_range is None:
            awg_range = (self.awg_min, self.awg_max)
        
        # Determine which parameters to optimize
        optimize_params = []
        param_bounds = []
        param_names = []
        
        # Build optimization parameter list based on what's fixed
        if fixed_power is None:
            optimize_params.append((power_range[0] + power_range[1]) / 2)
            param_bounds.append(power_range)
            param_names.append('power')
        
        if fixed_current is None:
            optimize_params.append((current_range[0] + current_range[1]) / 2)
            param_bounds.append(current_range)
            param_names.append('current')
            
        if fixed_voltage is None:
            optimize_params.append((voltage_range[0] + voltage_range[1]) / 2)
            param_bounds.append(voltage_range)
            param_names.append('voltage')
            
        if fixed_awg is None:
            optimize_params.append((awg_range[0] + awg_range[1]) / 2)
            param_bounds.append(awg_range)
            param_names.append('awg')
        
        if len(optimize_params) == 0:
            # All parameters fixed - just evaluate
            return self._evaluate_fixed_configuration(fixed_power, fixed_current, fixed_voltage, fixed_awg)
        
        # Create objective function wrapper for partial optimization
        def partial_objective(opt_params):
            # Reconstruct full parameter vector
            full_params = []
            opt_idx = 0
            
            if fixed_power is None:
                full_params.append(opt_params[opt_idx])
                opt_idx += 1
            else:
                full_params.append(fixed_power)
                
            if fixed_current is None:
                full_params.append(opt_params[opt_idx])
                opt_idx += 1
            else:
                full_params.append(fixed_current)
                
            if fixed_voltage is None:
                full_params.append(opt_params[opt_idx])
                opt_idx += 1
            else:
                full_params.append(fixed_voltage)
                
            if fixed_awg is None:
                full_params.append(opt_params[opt_idx])
                opt_idx += 1
            else:
                full_params.append(fixed_awg)
            
            return self.objective_function(full_params)
        
        print("Optimizing power supply specifications and wire gauge...")
        print(f"Weight limit: {self.weight_limit} lbs per coil")
        print(f"Power budget: {self.max_power} W")
        print(f"Radius: {self.radius*1000:.1f} mm")
        
        # Print what's being optimized vs fixed
        fixed_params = []
        if fixed_power is not None:
            fixed_params.append(f"Power: {fixed_power}W")
        if fixed_current is not None:
            fixed_params.append(f"Current: {fixed_current}A")
        if fixed_voltage is not None:
            fixed_params.append(f"Voltage: {fixed_voltage}V")
        if fixed_awg is not None:
            fixed_params.append(f"AWG: {fixed_awg}")
            
        if fixed_params:
            print(f"Fixed parameters: {', '.join(fixed_params)}")
        print(f"Optimizing: {', '.join(param_names)} using hybrid method")
        
        # Hybrid optimization: Global search followed by local refinement
        print("Running global optimization (differential_evolution)...")
        global_result = differential_evolution(
            func=partial_objective,
            bounds=param_bounds,
            maxiter=100,
            popsize=15,
            seed=42,
            disp=True
        )
        
        if global_result.success:
            print("Global optimization succeeded. Running local refinement...")
            # Use global result as starting point for local optimization
            result = minimize(
                fun=partial_objective,
                x0=global_result.x,
                method='L-BFGS-B',
                bounds=param_bounds,
                options={'disp': True, 'maxiter': 50}
            )
        else:
            print("Global optimization failed, using global result")
            result = global_result
        
        if result.success:
            # Reconstruct optimal parameters
            opt_idx = 0
            optimal_power = fixed_power if fixed_power is not None else result.x[opt_idx] 
            if fixed_power is None: opt_idx += 1
            
            optimal_current = fixed_current if fixed_current is not None else result.x[opt_idx]
            if fixed_current is None: opt_idx += 1
            
            optimal_voltage = fixed_voltage if fixed_voltage is not None else result.x[opt_idx]
            if fixed_voltage is None: opt_idx += 1
            
            optimal_awg = fixed_awg if fixed_awg is not None else result.x[opt_idx]
            
            optimal_b_field = -result.fun
            
            # Get detailed performance for optimal configuration
            power_supply_specs = {
                'Pmax': optimal_power, 
                'Imax': optimal_current, 
                'Vmax': optimal_voltage
            }
            
            performance = self.analyze_power_supply_performance(power_supply_specs, optimal_awg)
            
            if performance is not None:
                # Find the specific operating point that gives maximum B-field
                max_idx = np.argmax(performance['b_field'])
                
                recommendation = {
                    'success': True,
                    'optimal_power_supply': {
                        'Pmax': optimal_power,
                        'Imax': optimal_current, 
                        'Vmax': optimal_voltage
                    },
                    'optimal_awg': optimal_awg,
                    'max_b_field_tesla': optimal_b_field,
                    'operating_point': {
                        'turns_per_coil': performance['turns_single'][max_idx],
                        'current_A': performance['current'][max_idx],
                        'voltage_V': performance['voltage'][max_idx],
                        'power_W': performance['power'][max_idx],
                        'weight_single_coil_lbs': performance['weights'][max_idx],
                        'weight_both_coils_lbs': performance['weights'][max_idx] * 2,
                        'wire_length_m': performance['wire_length'][max_idx]
                    },
                    'optimization_info': {
                        'fixed_parameters': dict(zip(['power', 'current', 'voltage', 'awg'], 
                                                   [fixed_power, fixed_current, fixed_voltage, fixed_awg])),
                        'optimized_parameters': param_names,
                        'method': 'hybrid'
                    },
                    'full_performance': performance
                }
                
                self._print_recommendation(recommendation)
                return recommendation
            
        print(f"Optimization failed: {result.message}")
        return {'success': False, 'message': result.message}
    
    def _evaluate_fixed_configuration(self, fixed_power, fixed_current, fixed_voltage, fixed_awg):
        """Evaluate a completely fixed configuration."""
        print("Evaluating fixed configuration...")
        print(f"Fixed Power: {fixed_power}W, Current: {fixed_current}A, Voltage: {fixed_voltage}V, AWG: {fixed_awg}")
        
        power_supply_specs = {'Pmax': fixed_power, 'Imax': fixed_current, 'Vmax': fixed_voltage}
        performance = self.analyze_power_supply_performance(power_supply_specs, fixed_awg)
        
        if performance is not None and len(performance['b_field']) > 0:
            max_idx = np.argmax(performance['b_field'])
            max_b_field = performance['b_field'][max_idx]
            
            recommendation = {
                'success': True,
                'optimal_power_supply': power_supply_specs,
                'optimal_awg': fixed_awg,
                'max_b_field_tesla': max_b_field,
                'operating_point': {
                    'turns_per_coil': performance['turns_single'][max_idx],
                    'current_A': performance['current'][max_idx],
                    'voltage_V': performance['voltage'][max_idx],
                    'power_W': performance['power'][max_idx],
                    'weight_single_coil_lbs': performance['weights'][max_idx],
                    'weight_both_coils_lbs': performance['weights'][max_idx] * 2,
                    'wire_length_m': performance['wire_length'][max_idx]
                },
                'optimization_info': {
                    'fixed_parameters': {'power': fixed_power, 'current': fixed_current, 
                                       'voltage': fixed_voltage, 'awg': fixed_awg},
                    'optimized_parameters': [],
                    'method': 'evaluation'
                },
                'full_performance': performance
            }
            
            self._print_recommendation(recommendation)
            return recommendation
        else:
            return {'success': False, 'message': 'No valid configurations found for fixed parameters'}
    
    def _print_recommendation(self, rec):
        """Print formatted recommendation."""
        print("\n" + "="*60)
        print("OPTIMAL POWER SUPPLY AND WIRE CONFIGURATION")
        print("="*60)
        print(f"Power Supply Specifications:")
        print(f"  Maximum Power: {rec['optimal_power_supply']['Pmax']:.0f} W")
        print(f"  Maximum Current: {rec['optimal_power_supply']['Imax']:.1f} A")
        print(f"  Maximum Voltage: {rec['optimal_power_supply']['Vmax']:.1f} V")
        print(f"\nWire Specification:")
        print(f"  AWG Size: {rec['optimal_awg']:.1f} (use AWG{round(rec['optimal_awg'])})")
        print(f"\nOptimal Operating Point:")
        print(f"  Turns per Coil: {rec['operating_point']['turns_per_coil']:.0f}")
        print(f"  Operating Current: {rec['operating_point']['current_A']:.1f} A")
        print(f"  Operating Voltage: {rec['operating_point']['voltage_V']:.1f} V") 
        print(f"  Operating Power: {rec['operating_point']['power_W']:.0f} W")
        print(f"  Single Coil Weight: {rec['operating_point']['weight_single_coil_lbs']:.1f} lbs")
        print(f"  Both Coils Weight: {rec['operating_point']['weight_both_coils_lbs']:.1f} lbs")
        print(f"  Wire Length: {rec['operating_point']['wire_length_m']:.1f} m")
        print(f"\nPerformance:")
        print(f"  Maximum B-field: {rec['max_b_field_tesla']:.4f} T")
        
        # Show optimization info
        if 'optimization_info' in rec:
            print(f"  Optimization Method: {rec['optimization_info']['method']}")
            fixed = [f"{k}={v}" for k, v in rec['optimization_info']['fixed_parameters'].items() if v is not None]
            if fixed:
                print(f"  Fixed Parameters: {', '.join(fixed)}")
            if rec['optimization_info']['optimized_parameters']:
                print(f"  Optimized Parameters: {', '.join(rec['optimization_info']['optimized_parameters'])}")
        
        print("="*60)

def optimize_helmholtz(radius_m, weight_limit_lbs, max_power_w=10000,
                                   fixed_power=None, fixed_current=None, 
                                   fixed_voltage=None, fixed_awg=None):
    """
    Optimize electromagnet design using hybrid global+local optimization.
    
    Parameters
    ----------
    radius_m : float
        Coil radius in meters
    weight_limit_lbs : float
        Weight limit in pounds PER COIL
    max_power_w : float
        Maximum power budget in watts
    fixed_power : float, optional
        Fix power supply max power to this value
    fixed_current : float, optional
        Fix power supply max current to this value
    fixed_voltage : float, optional
        Fix power supply max voltage to this value
    fixed_awg : float, optional
        Fix wire AWG to this value
        
    Returns
    -------
    dict
        Optimal power supply specifications and configuration
        
    Examples
    --------
    # Optimize everything:
    result = optimize_electromagnet_realistic(0.125, 70, 20000)
    
    # Fix AWG, optimize power supply:
    result = optimize_electromagnet_realistic(0.125, 70, 20000, fixed_awg=12)
    
    # Fix power supply voltage, optimize rest:
    result = optimize_electromagnet_realistic(0.125, 70, 20000, fixed_voltage=500)
    """
    optimizer = FirstOrderOptimizer(
        radius=radius_m,
        weight_limit=weight_limit_lbs,  # Per coil
        max_power=max_power_w
    )
    
    return optimizer.optimize_electromagnet(
        current_range=(0, max_power_w),
        voltage_range=(0, max_power_w),
        fixed_power=fixed_power,
        fixed_current=fixed_current,
        fixed_voltage=fixed_voltage,
        fixed_awg=fixed_awg
    )

# Example usage
if __name__ == "__main__":
    # Test optimization
    result = optimize_helmholtz(
        radius_m=0.125,         # 125mm radius
        weight_limit_lbs=70,    # 70 lbs per coil (140 lbs total)
        max_power_w=20000       # 20kW power budget
    )
    
    if result['success']:
        ps = result['optimal_power_supply']
        op = result['operating_point']
        print(f"\nRecommendation: {ps['Pmax']:.0f}W/{ps['Imax']:.0f}A/{ps['Vmax']:.0f}V power supply")
        print(f"with AWG{round(result['optimal_awg'])} wire")
        print(f"Operating at {op['current_A']:.1f}A, {op['voltage_V']:.1f}V")
        print(f"Expected B-field: {result['max_b_field_tesla']:.3f} T")
    else:
        print("Optimization failed:", result['message'])