import numpy as np
import sys
sys.path.insert(0, '../code' )      # Add the SpinFlowSim folder where the syn.py and visu.py files are stored
import syn
import visu


def calcMicroPar(net):
    '''
    Computes the mean and standard deviation of velocities, as well as the mean values of radii and segment lengths, for a given vascular network. 

    MANDATORY INPUT PARAMETERS:
    - net: instance of pipenet class, storing an initialised vascular network.

    RETURNS:
    - vm: mean of the velocity distribution across capillary segments (units: mm/s).
    - vs: standard deviation of the velocity distribution across capillary segments (units: mm/s).
    - lm: mean of lengths of the network segments (units: mm).
    - rm: mean of the radii of the network segments (units: mm).
    '''
    connectivity = np.abs(net.connmat)  # connectivity matrix displaying how the nodes are connected

    # Access velocities array
    v = np.abs(net.velmat)
    varray = v[connectivity==1]         # exclude nodes with -1 connectivity due to symmetry
    vm = np.mean(varray)                # calculate mean velocity  (mm/s)
    vs = np.std(varray)                 # calculate standard deviation of velocity  (mm/s)

    # Access matrix of segment radii
    r = np.abs(net.radmat)
    marray = r[connectivity==1]         # exclude nodes with -1 connectivity due to symmetry
    rm = np.mean(marray)                # calculate mean radius for network (mm)

    # Access matrix of segement lengths 
    l = np.abs(net.lengthmat)
    larray = l[connectivity==1]         # exclude nodes with -1 connectivity due to symmetry
    lm = np.mean(larray)                # calculate mean segment length for network (mm)

    return vm, vs, lm, rm



def anaNetBran(rmat, nspins):
    '''
    Analyze the changes in x-, y-, z- trajectories of spins during a simulation to calculate the apparent network branching of a network.

    This function:
    1. Extracts the the x-, y-, z- components from a matrix stroing spin trajectories (units: mm). 
       This matrix is obtained by simulating a diffusion-weighted MRI measurement for a single b-value from a given number of protons flowing within a vascular network.
    2. Calculates the trajectory differences between successive timesteps.
    3. Counts the number of times each spin changes direction.
    4. Returns the average number of changes in direction a spin experiences during simulation time.


    MANDATORY INPUT PARAMETERS:
    - rmat: matrix storing the trajectories of the spins (units: mm). 
            It has size 3 x Ntsteps x Nspins , where Nspins is the number of requested spins,
            and Ntsetps is the number of simulation time points at a temporal resolution equal to dt.
    - nspins: number of flowing spins used for simulation.


    RETURNS:
    - anb: measures, on average, the number of segments spins travel through during a simulation.
    '''
    
    # Extract x, y, z trajectories from rmat
    x_traj = rmat[0, :, :]  # x-component of spins trajectory for all time steps
    y_traj = rmat[1, :, :]  # y-component of spins trajectory for all time steps
    z_traj = rmat[2, :, :]  # z-component of spins trajectory for all time steps
    
    # Transpose x-, y-, z- component trajectory matrices for looping
    x_traj_t, y_traj_t, z_traj_t = [np.transpose(traj) for traj in [x_traj, y_traj, z_traj]]


    def get_difference(positions):
        '''
        Calculate the trajectory difference between successive timesteps for all spin trajectories.
        '''
        # Calculate differences between successive timesteps
        return [positions[i, j + 1] - positions[i, j] 
                for i in range(positions.shape[0]) # Iterate over each spin
                for j in range(positions.shape[1] - 1)] # Iterate over each timestep 

    x_diff = get_difference(x_traj_t) # list of changes along x between consecutive timesteps for all spins
    y_diff = get_difference(y_traj_t) # list of changes along y between consecutive timesteps for all spins
    z_diff = get_difference(z_traj_t) # list of changes along z between consecutive timesteps for all spins


    #Count the number of times the spins changes direction
    changes = 0
    for t in range(1, len(x_diff)):
        previous = (x_diff[t - 1], y_diff[t - 1], z_diff[t - 1])  # Previous positions -x, -y, -z
        current = (x_diff[t], y_diff[t], z_diff[t])  # Current positions in -x, -y, -z
        if current != previous:  # Check if current positions differ from the previous ones
            changes += 1  # Increment the counter when a change is detected
        
    # Count the average number of changes for all spins across -x, -y, -z directions
    anb = changes / nspins  

    return anb
