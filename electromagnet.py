'''
Electromagnet Field Simulation

Provides classes and methods to simulate the magnetic field of an electromagnet.

Uses Finite Element Method to create and analyze the magnetic field of three
dimensional coils using the Biot-Savart law.

Classes:
    PowerSupply: Power parameters.
    Wire: Wire parameters.
    Electromagnet: Electromagnet methods and parameters.
    Field: Field methods and parameters.


Author: Aren Mowreader
Date: 08/01/2025
Oregon State University
Strongly Coupled Systems Lab

'''
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd

class PowerSupply:
    """
    Power supply parameters and specifications.
    
    Parameters
    ----------
    name : str
        Name or model of the power supply
    max_voltage : float (V)
        Maximum voltage output in volts
    max_current : float (A)
        Maximum current output in amperes
    max_power : float (W)
        Maximum power output in watts
    constant_current_resistance : float (ohms)
        Resistance threshold at which the power supply switches from constant
        current to crossover state
    constant_voltage_resistance : float (ohms)
        Resistance threshold at which the power supply switches from crossover
        state to constant voltage
    """
    def __init__(self, name, max_voltage, max_current, max_power):
        self.name = name
        self.max_voltage = max_voltage
        self.max_current = max_current
        self.max_power = max_power

        self.constant_current_resistance = (self.max_power/(self.max_current**2))
        self.constant_voltage_resistance = ((self.max_voltage**2)/self.max_power)

class Wire:
    """
    Wire parameters loaded from standard AWG specifications.
    
    Parameters
    ----------
    AWG : int
        American Wire Gauge size (6-56)
        
    Attributes
    ----------
    AWG : int
        American Wire Gauge size
    wire_data : DataFrame
        Complete wire data loaded from MWS_wire_data.csv
    diameter_min_m : float
        Minimum diameter in meters
    diameter_nom_m : float
        Nominal diameter in meters
    diameter_max_m : float
        Maximum diameter in meters
    resistance_min_ohm_per_ft : float
        Minimum resistance per foot in ohms
    resistance_nom_ohm_per_ft : float
        Nominal resistance per foot in ohms
    resistance_max_ohm_per_ft : float
        Maximum resistance per foot in ohms
    feet_per_pound : float
        Length of wire per pound in feet
    pounds_per_1000ft : float
        Weight per 1000 feet in pounds
    circular_mils : float
        Cross-sectional area in circular mils
    ohms_per_meter : float
        Resistance per meter in ohms
    lbs_per_meter : float
        Weight per meter in pounds
        
    Notes
    -----
    Loads wire specifications from MWS_wire_data.csv file.
    Diameter values are converted from inches to meters.
    """
    def __init__(self, AWG):
        # Load wire data from CSV
        import os
        # Get the directory where this script is located
        script_dir = os.path.dirname(os.path.abspath(__file__))
        csv_path = os.path.join(script_dir, 'MWS_wire_data.csv')
        self.wire_data = pd.read_csv(csv_path)
        
        # Find the row corresponding to the specified AWG
        wire_row = self.wire_data[self.wire_data['SIZE_AWG'] == AWG]
        
        if wire_row.empty:
            raise ValueError(f"AWG {AWG} not found in wire data")
        
        # Extract wire properties from the CSV data
        self.AWG = AWG
        self.diameter_min_m = wire_row['DIAMETER_MIN'].iloc[0] * 0.0254  # Convert inches to meters
        self.diameter_nom_m = wire_row['DIAMETER_NOM'].iloc[0] * 0.0254  # Convert inches to meters
        self.diameter_max_m = wire_row['DIAMETER_MAX'].iloc[0] * 0.0254  # Convert inches to meters
        self.resistance_min_ohm_per_ft = wire_row['RESISTANCE_MIN'].iloc[0]
        self.resistance_nom_ohm_per_ft = wire_row['RESISTANCE_NOM'].iloc[0]
        self.resistance_max_ohm_per_ft = wire_row['RESISTANCE_MAX'].iloc[0]
        self.feet_per_pound = wire_row['FEET_PER_POUND'].iloc[0]
        self.pounds_per_1000ft = wire_row['POUNDS_PER_1000FT'].iloc[0]
        self.circular_mils = wire_row['CIRCULAR_MILS'].iloc[0]
        self.ohms_per_meter = wire_row['ohms_per_meter'].iloc[0]
        self.lbs_per_meter = wire_row['lbs_per_meter'].iloc[0]

class Emag:
    """
    Electromagnet with configurable geometry and power parameters.
    
    Parameters
    ----------
    power_supply : PowerSupply
        Power supply object containing voltage, current, and power specifications
    clearance_radius : float
        Available radius for the coil in meters
    clearance_height : float
        Available height for the coil in meters
    single_coil_lbs_limit : float
        Maximum weight limit for a single coil in pounds
    wire : Wire
        Wire object containing AWG specifications and properties
        
    Attributes
    ----------
    geometry_1 : ndarray
         Array of shape (n_points, 3) containing x, y, z coordinates 
         of wire path for first coil geometry.
    geometry_2 : ndarray
         Array of shape (n_points, 3) containing x, y, z coordinates 
         of wire path for second coil geometry.
    delta_geometry_1 : ndarray
         Array of shape (n_points, 3) containing differences between 
         consecutive points in each geometry.
    delta_geometry_2 : ndarray
         Array of shape (n_points, 3) containing differences between 
         consecutive points in each geometry.
    net_turns : int
        Total number of turns. IMPORTANT: This is NET turns!
    geometry_1_turns : int
        Number of turns in first coil geometry.
    geometry_2_turns : int
        Number of turns in second coil geometry.
    turns_depth : int
        Number of turns in the radial direction
    turns_height : int
        Number of turns in the axial direction
    turns_remainder : int
        Remaining turns after filling depth and height
    wire_length : float
        Total length (m) of wire 
    wire_resistance : float
        Resistance (ohms) of wire at present length
    wire_weight : float
        Weight (lbs) of wire at present length
    current : float
        Remaining turns after filling depth and height
        Current (A) at present length of wire
    voltage : float
        Voltage (V) at present length of wire
    power : float
        Power (W) at present length of wire
    field : ndarray
        Array of shape (n_points, 3) containing the field values at each point.
    """
    def __init__(self, power_supply, wire):
        
        # Initial Conditions
        self.power_supply = power_supply
        self.wire = wire

        # Parameters filled by methods
        
        # Coil Geometry - Two separate arrays for two geometries
        self.geometry_1 = np.array([])  # First geometry array
        self.geometry_2 = np.array([])  # Second geometry array
        self.delta_geometry_1 = np.array([])  # Differences for geometry 1
        self.delta_geometry_2 = np.array([])  # Differences for geometry 2
        self.geometry_1_turns = 0
        self.geometry_2_turns = 0
        self.net_turns = 0       # Total number of turns in the coil
        self.helm_turns = 0      # Number of turns in one Helmholtz coil configuration
        self.turns_depth = 0     # Depth of the coil in turns
        self.turns_height = 0   # Height of the coil in turns
        self.turns_remainder = 0 # Remaining turns. (turns - turns_depth * turns_height)

        self.coil_length = 0
        self.weight = 0
        self.resistance = 0
        self.power = 0
        self.current = 0
        self.voltage = 0

        
        # Power Supply Responses to Coil
        self.current = 0
        self.voltage = 0
        self.power = 0

        # Field Parameters
        self.field = np.array([])
        self.b_field = np.array([])

    def calc_length(self):
        """
        Calculate the total length of the coil including all geometries and connections.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.coil_length by summing:
        1. Distances between consecutive points within each geometry
        2. Distances between the last point of each geometry and the first point 
           of the next geometry (connecting wires)
        """
        if self.geometry_1.size == 0 and self.geometry_2.size == 0:
            self.coil_length = 0
            return
            
        # Two separate arrays: geometry_1 and geometry_2
        total_length = 0
        self.calc_delta()
        
        # Calculate length within geometry 1
        if self.geometry_1.size > 0:
            total_length += np.sum(np.linalg.norm(self.delta_geometry_1, axis=1))
        
        # Calculate length within geometry 2
        if self.geometry_2.size > 0:
            total_length += np.sum(np.linalg.norm(self.delta_geometry_2, axis=1))
        
        # Calculate connecting length between geometries if both exist
        if self.geometry_1.size > 0 and self.geometry_2.size > 0:
            # Distance from last point of geometry 1 to first point of geometry 2
            last_point = self.geometry_1[-1]  # Last point of geometry 1
            first_point = self.geometry_2[0]  # First point of geometry 2
            connecting_length = np.linalg.norm(first_point - last_point)
            total_length += connecting_length
        
        self.coil_length = total_length

    def calc_delta(self):
        """
        Calculate the differences between consecutive points in each geometry.
        
        Returns
        -------
        None

        Notes
        -----
        Modifies self.delta_geometry by calculating the differences between 
        consecutive points in each geometry.
        """
        if self.geometry_1.size == 0 and self.geometry_2.size == 0:
            self.delta_geometry_1 = np.array([])
            self.delta_geometry_2 = np.array([])
            return
        
        self.delta_geometry_1 = np.diff(self.geometry_1, axis=0)
        self.delta_geometry_2 = np.diff(self.geometry_2, axis=0)

    def calc_weight(self):
        """
        Calculate the weight of wire in the electromagnet.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.weight by multiplying the coil length by the wire's weight
        per meter.
        """
        self.weight = self.coil_length * self.wire.lbs_per_meter

    def calc_resistance(self):
        """
        Calculate the resistance of wire in the electromagnet.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.resistance by multiplying the coil length by the wire's
        resistance per meter.
        """
        self.resistance = self.coil_length * self.wire.ohms_per_meter

    def calc_power(self):
        """
        Calculate the power draw of the electromagnet.

        This calculation is done assuming the given power supply is programmable
        and provides constant current (CC) to a specified threshold 
        (constant_current_resistance). From here it enters a crossover state
        where max power (P = IV) is being drawn from the power supply but neither the 
        maximum current nor maximum voltage is being drawn.
        When V = I*R, reaches the maximum voltage (max_voltage) the function 
        switches to constant voltage (CV) state.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.power using step function and power supply parameters.s
        """

        if self.resistance <= 0:
            self.power = 0
        elif self.resistance < self.power_supply.constant_current_resistance:
            self.power = (self.power_supply.max_current**2) * self.resistance
        elif self.resistance < self.power_supply.constant_voltage_resistance:
            self.power = self.power_supply.max_power
        else:
            self.power = self.power_supply.max_voltage**2/self.resistance

    def calc_current(self):
        """
        Calculate the current the Power Supply is capable of providing to the 
        electromagnet.

        This calculation is done assuming the given power supply is programmable
        and provides constant current (CC) to a specified threshold 
        (constant_current_resistance). From here it enters a crossover state
        where max power (P = IV) is being drawn from the power supply but neither the 
        maximum current nor maximum voltage is being drawn.
        When V = I*R, reaches the maximum voltage (max_voltage) the function 
        switches to constant voltage (CV) state.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.current using step function.

        """

        if self.resistance <= 0:
            self.current = 0
        elif self.resistance < self.power_supply.constant_current_resistance:
            self.current = self.power_supply.max_current
        elif self.resistance < self.power_supply.constant_voltage_resistance:
            self.current = np.sqrt(self.power_supply.max_power/self.resistance)
        else:
            self.current = self.power_supply.max_voltage/self.resistance

    def calc_voltage(self):
        """
        Calculate the voltage the Power Supply is capable of providing to the 
        electromagnet.

        This calculation is done assuming the given power supply is programmable
        and provides constant current (CC) to a specified threshold 
        (constant_current_resistance). From here it enters a crossover state
        where max power (P = IV) is being drawn from the power supply but neither the 
        maximum current nor maximum voltage is being drawn.
        When V = I*R, reaches the maximum voltage (max_voltage) the function 
        switches to constant voltage (CV) state.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.voltage using step function.
        """

        if self.resistance <= 0:
            self.voltage = 0
        elif self.resistance < self.power_supply.constant_current_resistance:
            self.voltage = self.power_supply.max_current * self.resistance
        elif self.resistance < self.power_supply.constant_voltage_resistance:
            self.voltage = np.sqrt(self.power_supply.max_power * self.resistance)
        else:
            self.voltage = self.power_supply.max_voltage

    def update_parameters(self):
        """
        Update the parameters of the electromagnet.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.coil_length, self.weight, self.resistance, self.power,
        self.current, and self.voltage.

        Calls calc_length, calc_weight, calc_resistance, calc_power, calc_current,
        and calc_voltage.
        """
        self.calc_length()
        self.calc_weight()
        self.calc_resistance()
        self.calc_power()
        self.calc_current()
        self.calc_voltage()

    def print_parameters(self):
        """
        Print the parameters of the electromagnet.
        """
        print(f"Coil Length: {self.coil_length:.3f} m")
        print(f"Weight: {self.weight:.3f} lbs")
        print(f"Resistance: {self.resistance:.3f} ohms")
        print(f"Power: {self.power:.3f} W")
        print(f"Current: {self.current:.3f} A")
        print(f"Voltage: {self.voltage:.3f} V")
        print(f"Net Turns: {self.net_turns}")
        print(f"Geometry 1 Turns: {self.geometry_1_turns}")
        print(f"Geometry 2 Turns: {self.geometry_2_turns}")
        

    def add_turn(self, center, radius, geometry_index=1, points_per_turn=1000):
        """
        Add one turn worth of points to coil geometry array.
        
        Creates a circular turn in the xy-plane at the specified z-coordinate.
        
        Parameters
        ----------
        center : ndarray, shape (3,)
            Center of the turn in x, y, z coordinates
        radius : float
            Radius of the turn in meters
        geometry_index : int, optional
            Index of the geometry to add the turn to (1 for geometry_1, 2 for geometry_2).
            Default is 1.
        points_per_turn : int, optional
            Number of points to use to represent the turn, default 1000
            
        Returns
        -------
        None
            
        Notes
        -----
        Modifies self.geometry_1 or self.geometry_2 by appending the new turn coordinates.
        The turn is created in the xy-plane with constant z-coordinate.
        
        Raises
        ------
        ValueError
            If center is not a 3D vector, radius is not positive, 
            points_per_turn is not positive, or geometry_index is not 1 or 2
        """
        if center.shape[0] != 3:
            raise ValueError("add_turn: Center must be a 3D vector")
        
        if radius <= 0:
            raise ValueError("add_turn: Radius must be positive")
        
        if points_per_turn <= 0:
            raise ValueError("add_turn: Points per turn must be positive")
        
        if geometry_index not in [1, 2]:
            raise ValueError("add_turn: geometry_index must be 1 or 2")
        
        # Create array of radians for the turn calculation
        theta = np.linspace(0, 2 * np.pi, points_per_turn, endpoint=False)

        # Create array of x, y, z coordinates for the turn
        x = center[0] + radius * np.cos(theta)
        y = center[1] + radius * np.sin(theta)
        z = np.full(points_per_turn, center[2])
        
        # Combine coordinates into a single array
        turn_coords = np.column_stack((x, y, z))
        
        # Determine which geometry to add to
        if geometry_index == 1:
            # Add to geometry_1
            if self.geometry_1.size == 0:
                self.geometry_1 = turn_coords
            else:
                self.geometry_1 = np.concatenate((self.geometry_1, turn_coords), axis=0)
            self.geometry_1_turns += 1
        else:  # geometry_index == 2
            # Add to geometry_2
            if self.geometry_2.size == 0:
                self.geometry_2 = turn_coords
            else:
                self.geometry_2 = np.concatenate((self.geometry_2, turn_coords), axis=0)
            self.geometry_2_turns += 1
        self.update_parameters()
        self.net_turns = self.geometry_1_turns + self.geometry_2_turns

    def add_mirrored_turns(self, radius, separation, center=np.array([0, 0, 0]),
                        points_per_turn=1000):
        """
        Add two spaced turns in mirrored configuration across xy-plane at
        specified z-coordinate.
        
        Creates two circular turns separated by the specified distance. One turn 
        is added to geometry_1 and one to geometry_2.
        
        Parameters
        ----------
        radius : float
            Radius of each turn in meters
        separation : float
            Distance between the centers of the two turns in meters
        center : ndarray, shape (3,), optional
            Center point between the two turns. Default is [0, 0, 0]
        points_per_turn : int, optional
            Number of points to use for each turn. Default is 1000
            
        Returns
        -------
        None
            
        Notes
        -----
        Modifies self.geometry_1 and self.geometry_2 by adding one turn to each.
        Top turn goes to geometry_1, bottom turn goes to geometry_2.
        
        Raises
        ------
        ValueError
            If radius or separation is not positive
        """
        if radius <= 0:
            raise ValueError("add_mirrored_turns: Radius must be positive")
        
        if separation <= 0:
             raise ValueError("add_mirrored_turns: Separation must be positive")
         
        if points_per_turn <= 0:
            raise ValueError("add_mirrored_turns: Points per turn must be positive")
        
         # Create array of radians for the turn calculation
        theta = np.linspace(0, 2 * np.pi, points_per_turn, endpoint=False)
        
        # Calculate centers for each turn (Helmholtz configuration)
        turn1_center = center + np.array([0, 0, separation/2])
        turn2_center = center + np.array([0, 0, -separation/2])
        
        # Create coordinates for turn 1 (top) - goes to geometry_1
        x1 = turn1_center[0] + radius * np.cos(theta)
        y1 = turn1_center[1] + radius * np.sin(theta)
        z1 = np.full(points_per_turn, turn1_center[2])
        top_coords = np.column_stack((x1, y1, z1))
        
        # Create coordinates for turn 2 (bottom) - goes to geometry_2
        x2 = turn2_center[0] + radius * np.cos(theta)
        y2 = turn2_center[1] + radius * np.sin(theta)
        z2 = np.full(points_per_turn, turn2_center[2])
        bottom_coords = np.column_stack((x2, y2, z2))
        
        # Add top turn to geometry_1
        if self.geometry_1.size == 0:
            self.geometry_1 = top_coords
        else:
            self.geometry_1 = np.concatenate((self.geometry_1, top_coords), axis=0)
        
        # Add bottom turn to geometry_2
        if self.geometry_2.size == 0:
            self.geometry_2 = bottom_coords
        else:
            self.geometry_2 = np.concatenate((self.geometry_2, bottom_coords), axis=0)
        
        self.update_parameters()
        self.geometry_1_turns += 1
        self.geometry_2_turns += 1
        self.helm_turns += 1
        self.net_turns += 2  # Added two turns
   
#VISUALIZATION AND DATA METHODS

    def plot_coil_geometry(self, index=None):
        """
        Plot the coil geometry in 3D.
        """
        # Create 3D figure
        fig = plt.figure(figsize=(10, 8))
        ax = fig.add_subplot(111, projection='3d')
        
        # Plot geometry 1
        if self.geometry_1.size > 0:
            ax.plot(self.geometry_1[:, 0], 
                   self.geometry_1[:, 1], 
                   self.geometry_1[:, 2], 
                   label='Geometry 1', color='blue')
        
        # Plot geometry 2
        if self.geometry_2.size > 0:
            ax.plot(self.geometry_2[:, 0], 
                   self.geometry_2[:, 1], 
                   self.geometry_2[:, 2], 
                   label='Geometry 2', color='red')
        
        # Set labels and title
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('Coil Geometry')
        
        # Add legend
        ax.legend()
        
        # Set equal aspect ratio for better visualization
        ax.set_box_aspect([1, 1, 1])
        
        plt.show()

    def plot_b_field(self, index=None):
        """
        Plot the b-field in 3D along with coil geometry.
        """
        # Create 3D figure
        fig = plt.figure(figsize=(12, 10))
        ax = fig.add_subplot(111, projection='3d')

        # Plot coil geometry first (so it appears behind the field vectors)
        if self.geometry_1.size > 0:
            ax.plot(self.geometry_1[:, 0], 
                   self.geometry_1[:, 1], 
                   self.geometry_1[:, 2], 
                   'b-', linewidth=3, label='Geometry 1', alpha=0.8)
        
        if self.geometry_2.size > 0:
            ax.plot(self.geometry_2[:, 0], 
                   self.geometry_2[:, 1], 
                   self.geometry_2[:, 2], 
                   'r-', linewidth=3, label='Geometry 2', alpha=0.8)

        # Plot b-field vectors
        if self.b_field.size > 0:
            # Reduce density for better visualization (plot every nth vector)
            step = max(1, len(self.field) // 1000)  # Show max 100 vectors
            ax.quiver(self.field[::step, 0], self.field[::step, 1], self.field[::step, 2],
                     self.b_field[::step, 0], self.b_field[::step, 1], self.b_field[::step, 2],
                     length=.05, normalize=True, color='green', alpha=0.6, label='Magnetic Field')
        
        # Set labels and title
        ax.set_xlabel('X (m)')
        ax.set_ylabel('Y (m)')
        ax.set_zlabel('Z (m)')
        ax.set_title('Magnetic Field and Coil Geometry')
        
        # Add legend
        ax.legend()
        
        # Set equal aspect ratio for better visualization
        ax.set_box_aspect([1, 1, 1])
        
        plt.show()
        plt.show()
        

    def get_geometry_info(self):
        """
        Get information about all geometries including their sizes.
        
        Returns
        -------
        dict
            Dictionary containing information about each geometry
        """
        if self.geometry_1.size == 0 and self.geometry_2.size == 0:
            return {"total_geometries": 0, "geometries": []}
        
        info = {
            "total_geometries": self.geometry_1.shape[0] + self.geometry_2.shape[0],
            "geometries": []
        }
        
        # Geometry 1 info
        if self.geometry_1.size > 0:
            geometry_info = {
                "index": 0,
                "name": "Geometry 1",
                "points": self.geometry_1.shape[0],
                "shape": self.geometry_1.shape,
                "bounds": {
                    "x_min": np.min(self.geometry_1[:, 0]),
                    "x_max": np.max(self.geometry_1[:, 0]),
                    "y_min": np.min(self.geometry_1[:, 1]),
                    "y_max": np.max(self.geometry_1[:, 1]),
                    "z_min": np.min(self.geometry_1[:, 2]),
                    "z_max": np.max(self.geometry_1[:, 2])
                }
            }
            info["geometries"].append(geometry_info)
        
        # Geometry 2 info
        if self.geometry_2.size > 0:
            geometry_info = {
                "index": 1,
                "name": "Geometry 2",
                "points": self.geometry_2.shape[0],
                "shape": self.geometry_2.shape,
                "bounds": {
                    "x_min": np.min(self.geometry_2[:, 0]),
                    "x_max": np.max(self.geometry_2[:, 0]),
                    "y_min": np.min(self.geometry_2[:, 1]),
                    "y_max": np.max(self.geometry_2[:, 1]),
                    "z_min": np.min(self.geometry_2[:, 2]),
                    "z_max": np.max(self.geometry_2[:, 2])
                }
            }
            info["geometries"].append(geometry_info)
        
        return info
    
    #FIELD METHODS

    def def_field(self, x_range, y_range, z_range, points_per_axis=20):
        """
        Define the field to analyze.

        Parameters
        ----------
        x_range : tuple
            Tuple of (x_min, x_max)
        y_range : tuple
            Tuple of (y_min, y_max)
        z_range : tuple
            Tuple of (z_min, z_max)
        points_per_axis : int, optional
            Number of points to use for each axis. Default is 20.

        Returns
        -------
        None

        Notes
        -----
        Modifies self.field by setting the field grid to analyze.
        Indexing specific field points from this array is not recommended.
        Every row of the field array is a point in the field grid.
        The field array is not a meshgrid, it is a flat array of points.
        This is done to simplify the numpy broadcasting operations.
        reshape the field array to a meshgrid if needed for plotting.

        Raises
        ------
        ValueError
            If x_range, y_range, or z_range is not a tuple
            If points_per_axis is not a positive integer

        """
        if not isinstance(x_range, tuple) or not isinstance(y_range, tuple) or not isinstance(z_range, tuple):
            raise ValueError("def_field: x_range, y_range, and z_range must be tuples")
        
        if points_per_axis <= 0:
            raise ValueError("def_field: points_per_axis must be a positive integer")
    

        #2d cases
        if x_range[0] == x_range[1]:
            # y-z plane (x = 0)
            y_points = np.linspace(y_range[0], y_range[1], points_per_axis)
            z_points = np.linspace(z_range[0], z_range[1], points_per_axis)
            Y, Z = np.meshgrid(y_points, z_points, indexing='ij')
            X = np.full(Y.flatten().shape, x_range[0])
            self.field = np.column_stack((X, Y.flatten(), Z.flatten()))
            return
        elif y_range[0] == y_range[1]:
            # x-z plane (y = 0)
            x_points = np.linspace(x_range[0], x_range[1], points_per_axis)
            z_points = np.linspace(z_range[0], z_range[1], points_per_axis)
            X, Z = np.meshgrid(x_points, z_points, indexing='ij')
            Y = np.full(X.flatten().shape, y_range[0])
            self.field = np.column_stack((X.flatten(), Y, Z.flatten()))
            return
        elif z_range[0] == z_range[1]:
            # x-y plane (z = 0)
            x_points = np.linspace(x_range[0], x_range[1], points_per_axis)
            y_points = np.linspace(y_range[0], y_range[1], points_per_axis)
            X, Y = np.meshgrid(x_points, y_points, indexing='ij')
            Z = np.full(X.flatten().shape, z_range[0])
            self.field = np.column_stack((X.flatten(), Y.flatten(), Z))
            return

        #3d case
        
        x_points = np.linspace(x_range[0], x_range[1], points_per_axis)
        y_points = np.linspace(y_range[0], y_range[1], points_per_axis)
        z_points = np.linspace(z_range[0], z_range[1], points_per_axis)
        
        # Create meshgrid and flatten to get all coordinate combinations
        X, Y, Z = np.meshgrid(x_points, y_points, z_points, indexing='ij')
        self.field = np.column_stack((X.flatten(), Y.flatten(), Z.flatten()))

    def calc_displacement_vectors(self):
        """
        Calculate displacement vectors from every coil point to every field point.
        
        Returns
        -------
        tuple
            (displacements_1, displacements_2) where each is a 3D array of shape
            (n_field_points, n_coil_points, 3) containing displacement vectors.
            
        Notes
        -----
        Uses broadcasting for efficient computation.
        np.newaxis is used to expand dimensions of arrays to match shapes for
        broadcasting.
        displacements[i, j, :] = field_point[i] - coil_point[j]
        The index i is the field point, and the index j is the coil point.
        The last three dimensions are the x, y, and z components of the displacement vector.
        The result is a 3D array of shape (n_field_points, n_coil_points, 3)
        """
        if self.field.size == 0:
            raise ValueError("calc_displacement_vectors: field must be defined")
        
        # Get field points: shape (n_field_points, 3)
        field_points = self.field
        
        # Calculate displacements for geometry_1
        if self.geometry_1.size > 0:
            # field_points[:, np.newaxis, :] has shape (n_field_points, 1, 3)
            # self.geometry_1[np.newaxis, :, :] has shape (1, n_coil_points, 3)
            # Result: (n_field_points, n_coil_points, 3)
            displacements_1 = field_points[:, np.newaxis, :] - self.geometry_1[np.newaxis, :, :]
        else:
            displacements_1 = np.array([])
        
        # Calculate displacements for geometry_2
        if self.geometry_2.size > 0:
            displacements_2 = field_points[:, np.newaxis, :] - self.geometry_2[np.newaxis, :, :]
        else:
            displacements_2 = np.array([])
        
        return displacements_1, displacements_2

    def calc_displacement_mag(self, displacements):
        """
        Calculate the magnitude of the displacement vectors.
        """
        return np.linalg.norm(displacements, axis=2)
    
    def calc_displacement_unit(self, displacements):
        """
        Calculate the unit vector of the displacement vectors.
        """
        return displacements / self.calc_displacement_mag(displacements)[:, :, np.newaxis]

    def calc_b_field(self):
        """
        Analyzes the Magnetic Field of the Electromagnet. Using the Biot-Savart law.

        Parameters
        ----------
        None

        Returns
        -------
        None
        
        Notes
        -----
        Modifies self.b_field by calculating the magnetic field at each point in the field grid.

        Raises
        ------
        ValueError
            If self.field is not defined

        """
        if self.field.size == 0:
            raise ValueError("def_b_field: field must be defined")
        
        # Calculate the magnetic field at each point in the field grid
        # b_field shape: (n_field_points, 3) for Bx, By, Bz components at each point
        self.b_field = np.zeros((self.field.shape[0], 3))
        mu_0 = 4 * np.pi * 10**-7
        const = (mu_0*self.current) / (4 * np.pi)
        
        # Calculate delta vectors first
        self.calc_delta()
        displacements_1, displacements_2 = self.calc_displacement_vectors()
        
        # Calculate cross products with proper broadcasting
        # delta_geometry_1: (n_coil_points_1, 3)
        # displacements_1: (n_field_points, n_coil_points_1, 3)
        # Need to broadcast delta_geometry_1 to match displacements_1
        if self.geometry_1.size > 0 and self.delta_geometry_1.size > 0:
            # Reshape delta_geometry_1 to (1, n_coil_points_1, 3) for broadcasting
            # this is done so every field point is crossed with every coil point
            # the result is a 3D array of shape (n_field_points, n_coil_points_1, 3)
            # broadcasting with np.newaxis saves memory and time by not creating
            # a new axis filled with identical values for each field point
            
            # This line is done to make the delta_geometry_1 array the same length as the displacements_1 array
            # It is assumed that the last dl will be the same as the next to last
            delta_1_broadcast = np.concatenate((self.delta_geometry_1, self.delta_geometry_1[-1:,:]))
            delta_1_broadcast = delta_1_broadcast[np.newaxis, :, :]

            geometry_1_cross = np.cross(delta_1_broadcast, displacements_1, axis=2)
            # Calculate Biot-Savart contribution for each coil segment
            # geometry_1_cross: (n_field_points, n_coil_points_1, 3)    
            # displacements_1_mag_cubed: (n_field_points, n_coil_points_1)
            displacements_1_mag_cubed = self.calc_displacement_mag(displacements_1)**3
            
            # Biot-Savart: B = (μ₀I/4π) * ∫ (dl × r̂) / r³
            # Use np.trapz to integrate over wire segments with dx = delta magnitudes
            integrand = geometry_1_cross / displacements_1_mag_cubed[:, :, np.newaxis]
            geometry_1_b_field = const * np.sum(integrand, axis=1)
        else:
            geometry_1_cross = np.array([])
            geometry_1_mag_cubed = np.array([])
            geometry_1_b_field = np.zeros((self.field.shape[0], 3))
        
        if self.geometry_2.size > 0 and self.delta_geometry_2.size > 0:
            # Reshape delta_geometry_2 to (1, n_coil_points_2, 3) for broadcasting
            # This line is done to make the delta_geometry_2 array the same length as the displacements_2 array
            # It is assumed that the last dl will be the same as the next to last
            delta_2_broadcast = np.concatenate((self.delta_geometry_2, self.delta_geometry_2[-1:,:]))
            delta_2_broadcast = delta_2_broadcast[np.newaxis, :, :]
            
            geometry_2_cross = np.cross(delta_2_broadcast, displacements_2, axis=2)
            
            # Calculate Biot-Savart contribution for geometry_2
            displacements_2_mag_cubed = self.calc_displacement_mag(displacements_2)**3
            
            # Use np.trapz to integrate over wire segments with dx = delta magnitudes
            integrand = geometry_2_cross / displacements_2_mag_cubed[:, :, np.newaxis]  # (n_field_points, n_coil_points_2, 3)
            geometry_2_b_field = const * np.sum(integrand, axis=1)
        else:
            geometry_2_cross = np.array([])
            geometry_2_mag_cubed = np.array([])
            geometry_2_b_field = np.zeros((self.field.shape[0], 3))
        
        # Combine contributions from both geometries
        if self.geometry_1.size > 0 and self.geometry_2.size > 0:
            self.b_field = geometry_1_b_field + geometry_2_b_field
        elif self.geometry_1.size > 0:
            self.b_field = geometry_1_b_field
        elif self.geometry_2.size > 0:
            self.b_field = geometry_2_b_field
        else:
            self.b_field = np.zeros((self.field.shape[0], 3))
       
    def calc_b_field_trapz(self):
        """
        Analyzes the Magnetic Field of the Electromagnet using the Biot-Savart law with np.trapz integration.

        Parameters
        ----------
        None

        Returns
        -------
        ndarray
            Magnetic field array with shape (n_field_points, 3) for Bx, By, Bz components
        
        Notes
        -----
        Uses np.trapz for numerical integration instead of np.sum.
        This method returns the b_field array but doesn't modify self.b_field.

        Raises
        ------
        ValueError
            If self.field is not defined

        """
        if self.field.size == 0:
            raise ValueError("calc_b_field_trapz: field must be defined")
        
        # Calculate the magnetic field at each point in the field grid
        # b_field shape: (n_field_points, 3) for Bx, By, Bz components at each point
        b_field = np.zeros((self.field.shape[0], 3))
        mu_0 = 4 * np.pi * 10**-7
        const = (mu_0*self.current) / (4 * np.pi)
        
        # Calculate delta vectors first
        self.calc_delta()
        displacements_1, displacements_2 = self.calc_displacement_vectors()
        
        # Calculate cross products with proper broadcasting
        if self.geometry_1.size > 0 and self.delta_geometry_1.size > 0:
            # Reshape delta_geometry_1 to (1, n_coil_points_1, 3) for broadcasting
            delta_1_broadcast = np.concatenate((self.delta_geometry_1, self.delta_geometry_1[-1:,:]))
            delta_1_broadcast = delta_1_broadcast[np.newaxis, :, :]

            geometry_1_cross = np.cross(delta_1_broadcast, displacements_1, axis=2)
            # Calculate Biot-Savart contribution for each coil segment
            displacements_1_mag_cubed = self.calc_displacement_mag(displacements_1)**3
            
            # Biot-Savart: B = (μ₀I/4π) * ∫ (dl × r̂) / r³
            # Use np.trapz to integrate over wire segments
            integrand = geometry_1_cross / displacements_1_mag_cubed[:, :, np.newaxis]
            
            # Integrate using np.trapz for each field point and component
            geometry_1_b_field = np.zeros((self.field.shape[0], 3))
            for i in range(self.field.shape[0]):
                for j in range(3):
                    geometry_1_b_field[i, j] = np.trapz(integrand[i, :, j])
            
            geometry_1_b_field = const * geometry_1_b_field
        else:
            geometry_1_b_field = np.zeros((self.field.shape[0], 3))
        
        if self.geometry_2.size > 0 and self.delta_geometry_2.size > 0:
            # Reshape delta_geometry_2 to (1, n_coil_points_2, 3) for broadcasting
            delta_2_broadcast = np.concatenate((self.delta_geometry_2, self.delta_geometry_2[-1:,:]))
            delta_2_broadcast = delta_2_broadcast[np.newaxis, :, :]
            
            geometry_2_cross = np.cross(delta_2_broadcast, displacements_2, axis=2)
            
            # Calculate Biot-Savart contribution for geometry_2
            displacements_2_mag_cubed = self.calc_displacement_mag(displacements_2)**3
            
            # Use np.trapz to integrate over wire segments
            integrand = geometry_2_cross / displacements_2_mag_cubed[:, :, np.newaxis]
            
            # Integrate using np.trapz for each field point and component
            geometry_2_b_field = np.zeros((self.field.shape[0], 3))
            for i in range(self.field.shape[0]):
                for j in range(3):
                    geometry_2_b_field[i, j] = np.trapz(integrand[i, :, j])
            
            geometry_2_b_field = const * geometry_2_b_field
        else:
            geometry_2_b_field = np.zeros((self.field.shape[0], 3))
        
        # Combine contributions from both geometries
        if self.geometry_1.size > 0 and self.geometry_2.size > 0:
            b_field = geometry_1_b_field + geometry_2_b_field
        elif self.geometry_1.size > 0:
            b_field = geometry_1_b_field
        elif self.geometry_2.size > 0:
            b_field = geometry_2_b_field
        else:
            b_field = np.zeros((self.field.shape[0], 3))
        
        return b_field
       
    def reshape_field_for_plotting(self, points_per_axis):
        """
        Reshape the field array back to meshgrid format for plotting.
        
        Parameters
        ----------
        points_per_axis : int
            Number of points used per axis when creating the field
            
        Returns
        -------
        tuple
            (X, Y, Z) meshgrid arrays for plotting
            
        Notes
        -----
        Converts the flattened field array back to 3D meshgrid format.
        This is useful for plotting field data with plt.contour, plt.contourf, etc.
        """
        if self.field.size == 0:
            raise ValueError("reshape_field_for_plotting: field must be defined")
        
        # Extract x, y, z coordinates from flattened field
        x_coords = self.field[:, 0]
        y_coords = self.field[:, 1]
        z_coords = self.field[:, 2]
        
        # Reshape to 3D arrays
        X = x_coords.reshape(points_per_axis, points_per_axis, points_per_axis)
        Y = y_coords.reshape(points_per_axis, points_per_axis, points_per_axis)
        Z = z_coords.reshape(points_per_axis, points_per_axis, points_per_axis)
        
        return X, Y, Z
    
    def reshape_b_field_for_plotting(self, points_per_axis):
        """
        Reshape the b_field array back to meshgrid format for plotting.
        
        Parameters
        ----------
        points_per_axis : int
            Number of points used per axis when creating the field
            
        Returns
        -------
        tuple
            (Bx, By, Bz) meshgrid arrays for plotting magnetic field components
            
        Notes
        -----
        Converts the flattened b_field array back to 3D meshgrid format.
        Each component (Bx, By, Bz) is reshaped separately.
        """
        if self.b_field.size == 0:
            raise ValueError("reshape_b_field_for_plotting: b_field must be defined")
        
        # Extract Bx, By, Bz components from flattened b_field
        Bx = self.b_field[:, 0]
        By = self.b_field[:, 1]
        Bz = self.b_field[:, 2]
        
        # Reshape to 3D arrays
        Bx_3d = Bx.reshape(points_per_axis, points_per_axis, points_per_axis)
        By_3d = By.reshape(points_per_axis, points_per_axis, points_per_axis)
        Bz_3d = Bz.reshape(points_per_axis, points_per_axis, points_per_axis)
        
        return Bx_3d, By_3d, Bz_3d
        
        
        