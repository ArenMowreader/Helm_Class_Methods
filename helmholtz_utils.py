import numpy as np
import pandas as pd

def b_helm_current(I, n_single_coil, r, mu_0 = 4 * np.pi * 10**-7):
    """
    Calculate the magnetic field strength of a Helmholtz coil based on current.
    
    Parameters:
    -----------
    I : float
        Current in amperes
    n_both_coils : int or array
        Number of turns on both coils
    r : float
        Radius of the coil in meters
    mu_0 : float, optional
        Permeability of free space (default: 4π × 10^-7 H/m)
        
    Returns:
    --------
    float or array
        Magnetic field strength in Tesla
    """
    return (8*mu_0*I*n_single_coil)/(r*np.sqrt(125))

def b_helm_power(power, n_single_coil, r, ohms_per_meter, mu_0 = 4 * np.pi * 10**-7):
    """
    Calculate the magnetic field strength of a Helmholtz coil based on power.
    
    Parameters:
    -----------
    power : float
        Power in watts
    n_both_coils : int or array
        Number of turns on both coils
    r : float
        Radius of the coil in meters
    ohms_per_meter : float
        Resistance per meter of wire
    mu_0 : float, optional
        Permeability of free space (default: 4π × 10^-7 H/m)
        
    Returns:
    --------
    float or array
        Magnetic field strength in Tesla
    """
    return ((8 * mu_0) / r**(3/2)) * ((n_single_coil * power) / (np.pi * ohms_per_meter * 500))**(1/2)

def n_both_coils_crit_Imax(power_supply, ohms_per_meter, r):
    """
    Calculate the critical number of turns for a power supply based on the 
    maximum current.

    IMPORTANT: This is the number of turns for BOTH coils as it is assumed the
    coils are run in series. Odd numbers will break the symetry and do not make
    sense to construct a coil.
    
    Parameters:
    -----------
    power_supply : dict
        Dictionary containing power supply specifications (Pmax, Imax, Vmax)
    ohms_per_meter : float
        Resistance per meter of wire
    r : float
        Radius of the coil in meters
        
    Returns:
    --------
    int
        Critical number of turns    
    """
    return round(power_supply['Pmax'] / (2 * np.pi * r * ohms_per_meter * power_supply['Imax']**2))

def n_both_coils_crit_Vmax(power_supply, ohms_per_meter, r):
    """
    Calculate the critical number of turns for a power supply based on the
    maximum voltage.

    IMPORTANT: This is the number of turns for BOTH coils as it is assumed the
    coils are run in series. Odd numbers will break the symetry and do not make
    sense to construct a helmholtz coil.
    
    Parameters:
    -----------
    power_supply : dict
        Dictionary containing power supply specifications (Pmax, Imax, Vmax)
    ohms_per_meter : float
        Resistance per meter of wire
    r : float
        Radius of the coil in meters
        
    Returns:
    --------
    int
        Critical number of turns
    """
    return round(power_supply['Vmax']**2 / (power_supply['Pmax'] * 2 * np.pi * r * ohms_per_meter))

def analyze_power_supply(power_supply, ohms_per_meter, n_single_coil, r):
    """
    Analyze power supply performance for different numbers of turns.
    
    Parameters:
    -----------
    power_supply : dict or pd.Series
        Dictionary containing power supply specifications (Pmax, Imax, Vmax)
        or pd.Series containing power supply specifications (Pmax, Imax, Vmax)
    ohms_per_meter : float
        Resistance per meter of wire
    n : np.array
        Array of number of turns to analyze.

        IMPORTANT: This is the number of turns on a SINGLE coil.
    r : float
        Radius of the coil in meters
        
    Returns:
    --------
    pd.DataFrame
        DataFrame with (V, I, B, power) arrays for each number of turns (even only)
    """
    # Initialize arrays
    power = np.zeros_like(n_single_coil, dtype=float)
    I = np.zeros_like(n_single_coil, dtype=float)
    V = np.zeros_like(n_single_coil, dtype=float)
    
    # Calculate critical points
    n_crit1 = round(power_supply['Pmax'] / (4 * np.pi * r * ohms_per_meter * 
                                            power_supply['Imax']**2))
    n_crit2 = round(power_supply['Vmax']**2 / (power_supply['Pmax'] * 4 * np.pi *
                                                r * ohms_per_meter))
    
    power[0:n_crit1] = power_supply['Imax']**2 * (4 * np.pi * r * ohms_per_meter 
                                                  * n_single_coil[0:n_crit1])
    I[0:n_crit1] = power_supply['Imax']
    V[0:n_crit1] = I[0:n_crit1] * (4 * np.pi * r * ohms_per_meter * n_single_coil[0:n_crit1])

    power[n_crit1:n_crit2] = power_supply['Pmax'] 
    I[n_crit1:n_crit2] = np.sqrt(power[n_crit1:n_crit2] / (4 * np.pi * r *
                                                            ohms_per_meter * 
                                                            n_single_coil[n_crit1:n_crit2]))
    V[n_crit1:n_crit2] = np.sqrt(power[n_crit1:n_crit2] * (4 * np.pi * r * 
                                                           ohms_per_meter * 
                                                           n_single_coil[n_crit1:n_crit2]))

    power[n_crit2:] = power_supply['Vmax']**2 / (4 * np.pi * r * ohms_per_meter *
                                                  n_single_coil[n_crit2:])
    I[n_crit2:] = power_supply['Vmax'] / (4 * np.pi * r * ohms_per_meter *
                                           n_single_coil[n_crit2:])
    V[n_crit2:] = power_supply['Vmax']

    # Calculate B field
    B_current = b_helm_current(I, n_single_coil, r)
    B_power = b_helm_power(power, n_single_coil, r, ohms_per_meter)
    
    return pd.DataFrame({
        'V': V,
        'I': I, 
        'B_current': B_current,
        'B_power': B_power,
        'power': power
    })

def weight(n_one_coil, r, lbs_per_meter):
    """
    Calculate the weight of one coil.
    
    IMPORTANT: This is the weight of ONE coil.
    If using the n that has been calculated from something like the current,
    voltage, or power, the input n should be divided by 2. The n from those
    calculations is the turns on both coils as it is assumed the coils are 
    run in series.
    Parameters:
    -----------
    n : int
        Number of turns
    r : float
        Radius of the coil in meters
    lbs_per_meter : float
        Weight per meter of wire in pounds
        
    Returns:
    --------
    float
        Total weight of the coil in pounds
    """
    return lbs_per_meter * n_one_coil * 2 * np.pi * r

def weight_limit_n_one_coil(limit, r, lbs_per_meter):
    """
    Calculate the maximum number of turns before exceeding weight limit.
    IMPORTANT: This is the weight of ONE coil
    If using the n that has been calculated from something like the current,
    voltage, or power, the input n should be divided by 2. The n from those
    calculations is the turns on both coils as it is assumed the coils are 
    run in series.
        
    Parameters:
    -----------
    limit : float
        Weight limit in pounds
    r : float
        Radius of the coil in meters
    lbs_per_meter : float
        Weight per meter of wire in pounds
        
    Returns:
    --------
    int
        Maximum number of turns before exceeding weight limit
    """
    n = 1
    while weight(n, r, lbs_per_meter) < limit:
        n += 1
    return n

def wire_length(n, r):
    """
    Calculate the total length of wire needed.
    
    Parameters:
    -----------
    n : int
        Number of turns
    r : float
        Radius of the coil in meters
        
    Returns:
    --------
    float
        Total length of wire in meters
    """
    return n * 2 * np.pi * r

def wire_length_advanced(radius, turns_per_layer, layers, wire_diameter):
    """
    Calculate the total length of wire needed for variable dimension coil
    IMPORTANT: This is the length of BOTH coil.
    
    Parameters:
    -----------
    radius: float
        Radius of the coil in meters
    turns_per_layer: int
        Number of turns in one layer of the spiral
    layers: int
        Number of layers in the total coil
    wire_diameter: float
        Diameter of the wire in meters
    """
    for i in range(layers):
        length = 2 * np.pi * turns_per_layer * (radius + wire_diameter * i)
        if i == 0:
            length_total = length
        else:
            length_total += length
    return length_total*2

def wire_length_square(radius, turns, wire_diameter):
    """
    Calculate the total length of wire assuming square coil geometry.
    
    turns per layer = layers


    
    Parameters:
    -----------
    radius: float
        Radius of the coil in meters
    turns: int
        Number of turns
    wire_diameter: float
        Diameter of the wire in meters

    Returns:
    --------
    float
        Total length of wire in meters
    """
    dimension = np.round(np.sqrt(turns)) #round to the nearest integer
    for i in range(dimension):
        length = 2 * np.pi * (radius + wire_diameter * i) * dimension
        if i == 0:
            length_total = length
        else:
            length_total += length
    return length_total*2

#def wire_length_advanced(radius, turns, wire_diameter):

def coil_layout(turns):
    """Calculate the dimensions of Helmholtz coil trying to maintain as square
    of a geometry as possible
    
    Parameters:
    -----------
    radius: float
        Radius of the coil in meters
    turns: int
        Number of turns
    wire_diameter: float
        Diameter of the wire in meters

    Returns: 
    --------
    turns_per_layer: int
        Number of turns in one layer of the spiral
    layers: int
        Number of layers in the total coil
    remainder: int
        Remaining turns for the last layer
    """
    turns_per_layer = np.floor(np.sqrt(turns)) #round down to the nearest integer
    layers = np.floor(turns / turns_per_layer)
    remainder = turns - turns_per_layer * layers
    return turns_per_layer, layers, remainder

def power_at_length(Pmax, Imax, Vmax, length, ohms_per_meter):
    """
    Calculate the power at a given length of wire.
    Parameters:
    -----------
    power: float
        Power in watts
    I: float
        Current in amperes
    V: float
        Voltage in volts
    length: float
        Length of wire in meters
    ohms_per_meter: float
        Resistance per meter of wire
        
    Returns:
    --------
    float
        Power in watts
    """
    #calculate the critical lengths for the power supply
    length_crit1 = Pmax / (ohms_per_meter * Imax**2)
    length_crit2 = Vmax**2 / (ohms_per_meter * Pmax)
    #calculate the voltage drop at the length

    if length < length_crit1:
        return Imax**2 * length * ohms_per_meter
    elif length < length_crit2:
        return Pmax
    else:
        return Vmax**2 / (ohms_per_meter * length)

def resistance_at_length(ohms_per_meter, length):
    """
    Calculate the resistance at a given length of wire.
    """
    return ohms_per_meter * length

def get_wire_data(data, start_size, end_size, step_size=1):

    """
    Get wire data for a range of AWG sizes.
    
    Parameters:
    -----------
    data : pd.DataFrame
        DataFrame containing wire data
    start_size : int
        Starting AWG size
    end_size : int
        Ending AWG size
    step_size : int
        Step size between AWG sizes
        
    Returns:
    --------
    tuple
        (ohms_per_meter, lbs_per_meter, wire_diameter) arrays for the specified AWG sizes
    """
    n = 0
    ohms_per_meter = np.zeros(len(np.arange(start_size, end_size+step_size, 
                                            step_size)))
    lbs_per_meter = np.zeros(len(np.arange(start_size, end_size+step_size, 
                                           step_size)))
    wire_diameter = np.zeros(len(np.arange(start_size, end_size+step_size, 
                                           step_size)))
    for i in np.arange(start_size, end_size+step_size, step_size):
        ohms_per_meter[n] = data.loc[data['SIZE_AWG'] == i, 'ohms_per_meter'].iloc[0]
        lbs_per_meter[n] = data.loc[data['SIZE_AWG'] == i, 'lbs_per_meter'].iloc[0]
        wire_diameter[n] = data.loc[data['SIZE_AWG'] == i, 'DIAMETER_NOM'].iloc[0]
        n += 1
    return ohms_per_meter, lbs_per_meter , wire_diameter


def spiral(radius, wire_diameter, turns_per_layer, layers, dt, z_offset):
    '''
    Function that paramaterizes a spiral.
    Inputs:
        radius: radius of the spiral
    
        wire_diameter: Diameter of the wire to be used.
    
        turns_per_layer: number of turns in one layer of the spiral
    
        layers: number of layers in the total coil. Input the actual
            number of layers, not the number of layers minus 1.
    
        dt: number of points to paramaterize the spiral. More points will
            make the spiral smoother and the total approximation more accurate.
            But also more computationally expensive.
    
        z_offset: offset of the spiral in the z direction.
    
    Outputs:
        array of x, y, z values of the spiral. 
        Dimensions: (dt*layers*turns_per_layer, 3)
    '''
    for i in range(layers):
        t = np.linspace(0, 2*np.pi*turns_per_layer, dt) # from 0 to 2pi*turns_per_layer with dt points
        x = (radius+(wire_diameter*i))*np.cos(t)
        y = (radius+(wire_diameter*i))*np.sin(t)
        z = z_offset+t*(wire_diameter/(2*np.pi))
        if i == 0:
            x_total = x
            y_total = y
            z_total = z
        else: # add the new layer to the total array
            x_total = np.concatenate((x_total, x))
            y_total = np.concatenate((y_total, y))
            z_total = np.concatenate((z_total, z))
    return np.column_stack((x_total, y_total, z_total))

#Paramaterize dl of the spiral.
def dl(radius, wire_diameter, turns_per_layer, layers, dt):
    '''
    Function that paramaterizes dl of the spiral.
    Inputs:
        radius: radius of the spiral
    
        wire_diameter: Diameter of the wire to be used.

        turns_per_layer: number of turns in one layer of the spiral
        
        layers: number of layers in the total coil. Input the actual
            number of layers.
        
        dt: number of points to paramaterize the spiral

    Outputs:
        array of dl values of the spiral. Array corresponds to spiral array 
        function. spiral[0] is the position vector while I*dl[0] is equal to 
        the velocity of one small charge element(dq) in the spiral at that 
        position. 
        
        Dimensions: (dt*layers*turns_per_layer, 3)
    '''
    for i in range(layers):
        t = np.linspace(0, 2*np.pi*turns_per_layer, dt) 
        x = -(radius+(wire_diameter*i))*np.sin(t)  # Fixed: use i instead of layers
        y = (radius+(wire_diameter*i))*np.cos(t)   # Fixed: use i instead of layers
        #array filled with wire_diameter/(2*pi) to match the shape of x and y
        z = np.full(dt, wire_diameter/(2*np.pi))
        if i == 0:
            x_total = x
            y_total = y
            z_total = z
        else:
            x_total = np.concatenate((x_total, x))
            y_total = np.concatenate((y_total, y))
            z_total = np.concatenate((z_total, z))
    return np.column_stack((x_total, y_total, z_total))

def displacement_vector(measured_position, charge_position):
    '''
    Function that calculates the displacement vector from a position to all 
    points on the spiral.
    Inputs:
        spiral: array of x, y, z values of the spiral. 
        Dimensions: (n, 3)
        
        position: position of the point being measured. Can be single point or
        array of points.
        Dimensions: (n, 3)
    
    Outputs: displacement vector from the spiral to all points in the grid.
    Dimensions: 
    '''
    #check if the input is of the form (n, 3)
    if measured_position.shape[1] != 3:
        raise ValueError('Input array measured_position must have three columns')
    if charge_position.shape[1] != 3:
        raise ValueError('Input array charge_position must have three columns')
    for i in range(len(measured_position)):
        x = measured_position[i,0] - charge_position[:,0] #
        y = measured_position[i,1] - charge_position[:,1]
        z = measured_position[i,2] - charge_position[:,2]
        if i == 0:
            x_total = x
            y_total = y
            z_total = z
        else:
            x_total = np.concatenate((x_total, x))
            y_total = np.concatenate((y_total, y))
            z_total = np.concatenate((z_total, z))
    return np.column_stack((x_total, y_total, z_total))

def vector_magnitude_3d(vector):
    '''
    Function that calculates the magnitude of any given three dimensional vector.
    Inputs:
        vector: array of x, y, z values of the vector.
        Dimensions: (n, 3)
        
    
    Outputs: magnitude of the vector or set of vectors.
    Dimensions: (n, 1) This will correspond to the displacement vector array.
    index i of output is equivalent to the magnitude of three dimensional vector 
    at index i of the input array.
    '''
    #check if the input is of the form (n, 3)
    if vector.shape[1] != 3:
        raise ValueError('Input array must have three columns')
    return np.sqrt(vector[:, 0]**2 + 
        vector[:, 1]**2 + vector[:, 2]**2) # shape 


def Helm(radius, wire_diameter, turns_per_layer, layers, dt, z_offset,
          grid_positions, I, mu_0 = 4 * np.pi * 10**-7):
    '''Desription: Function that calculates the magnetic field produced by two coils
    in a helmholtz configuration. The magnetic field is calculated at an array of points
    in space with the Biot-Savart law. Is possible to that it will calculate field given
    any line of current and its dl vector. Further testing is needed to confirm this.
    Inputs:
        radius: radius of the spiral

        wire_diameter: Diameter of the wire to be used.

        turns_per_layer: number of turns in one layer of the spiral

        layers: number of layers in the total coil. Input the actual
            number of layers, not the number of layers minus 1.

        dt: number of points to paramaterize the spiral. More points will
            make the spiral smoother and the total approximation more accurate.
            But also more computationally expensive.

        z_offset: offset of the spiral in the z direction.

        measured_position: array of x, y, z values of the points in space where the
        magnetic field is to be calculated.

        I: current in the wire (A)

    Returns:
        B_vec: array of x, y, z values of the magnetic field at the points in space.
        Indexing corresponds to the index of the measured_position array.
        Both will need to be reshaped to the shape of the grid to be plotted. 
        Dimensions: (m, 3)
    '''
    constants = I*mu_0/(4*np.pi) 
  

    #spiral1: array of x, y, z values of the first spiral.
    #Dimensions: (n, 3) VERY IMPOTANT that all n values are the same for all inputs.
    spiral1 = spiral(radius, wire_diameter, turns_per_layer, layers, dt, z_offset)
    
    #spiral2: array of x, y, z values of the second spiral.
    #Dimensions: (n, 3)
    spiral2 = spiral(radius, wire_diameter, -turns_per_layer, layers, dt, -z_offset)

    #dl1: array of x, y, z values of the dl vector of the first spiral. Derivative
    #of the position vector of spiral1.
    #Dimensions: (n, 3)
    dl1 = dl(radius, wire_diameter, turns_per_layer, layers, dt)

    #dl2: array of x, y, z values of the dl vector of the second spiral. Derivative
    #of the position vector of spiral2.
    #Dimensions: (n, 3)
    dl2 = dl(radius, wire_diameter, -turns_per_layer, layers, dt)

    # Calculate delta for both spirals (they should be identical due to same geometry)
    # Use dl vectors minus the last position (since last point has no "next" point)
    # This gives us the infinitesimal displacement vectors for integration
    delta = dl1[:-1]  # Remove last element since it has no next point
    
    #initialize the array to store the magnetic field at each point.
    # this index corresponds to the index of the measured_position array.
    # Stores vector components of the magnetic field at each point.
    B_vec = np.zeros((grid_positions.shape[0], 3))

    #for loop that goes through all the points in the measured_position array
    # and calculates the magnetic field at that point
    for i in range(len(grid_positions)):
        #takes the first point in the measured_position array. .reshape(1,-1) is used
        #to make the array 2d. (1,-1) means one row and as many columns as needed.
        #this satisfies the input requirements of the displacement_vector function.
        displacement1 = displacement_vector(grid_positions[i].reshape(1,-1), spiral1)
        displacement2 = displacement_vector(grid_positions[i].reshape(1,-1), spiral2)
        # calculate the magnitude cubed of the displacement vectors.
        # it is uneccessary to store the magnitude as it is a intermediate step.
        # this index corresponds to the index of the displacement vector array,
        # and the position vector array for the spiral. It should be 1d though as it is
        # a scalar value. When dividing by remember to do to each element of vector.

        mag_cubed1 = (vector_magnitude_3d(displacement1))**3
        mag_cubed2 = (vector_magnitude_3d(displacement2))**3
        # transform mag_cubed to a 2d (n,3) array to match the shape of the displacement vector
        # array. This is done to make the division of the displacement vector by the magnitude
        # cubed easier.
        mag_cubed1 = np.column_stack((mag_cubed1, mag_cubed1, mag_cubed1))
        mag_cubed2 = np.column_stack((mag_cubed2, mag_cubed2, mag_cubed2))
        


        # calculate the cross product of dl and the displacement vector
        # this index corresponds to the index of the displacement vector array,
        # and the position vector array for the spiral. It is the cross product for each
        # point in the spiral with the point currently being measured.
        
        cross1 = np.cross(dl1, displacement1)
        cross2 = np.cross(dl2, displacement2)
        
        # calculate the magnetic field at the point being measured.
        # We need to integrate along the wire path (axis=0)

        
        # Integrate along the wire path component-wise
        # cross1/mag_cubed1 has shape (n, 3), we integrate each component separately
        # to get a 3-element vector for the B-field at this point
        B1 = np.zeros(3)
        B2 = np.zeros(3)
        
        for component in range(3):
            B1[component] = constants * np.trapz(cross1[:, component]/
                                                 mag_cubed1[:, component],
                                                   dx=delta[:, component])
            B2[component] = constants * np.trapz(cross2[:, component]/
                                                 mag_cubed2[:, component], 
                                                 dx=delta[:, component])
        B_vec[i] = B1 + B2

    return B_vec

def grid_positions_xz(width, height, spacing=.05):
    """
    Creates grid positions for use in the Helm() function.
    Inputs:
        width: width of the grid in the x direction
        height: height of the grid in the z direction
        spacing: spacing of the grid in the x and z directions
    Outputs:
        X: 2D array of x coordinates (for reshaping)
        Z: 2D array of z coordinates (for reshaping)
        grid_positions: flattened array of x, y, z values for Helm() function
    """
    x = np.arange(-width/2, width/2, spacing)
    z = np.arange(-height/2, height/2, spacing)
    X, Z = np.meshgrid(x, z)
    grid_positions = np.column_stack((X.flatten(), np.full(X.size, 0), Z.flatten()))
    return X, Z, grid_positions

