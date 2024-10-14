import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib


def calculate_velocity(trajectories):
    """
    Function to calculate spin velocities from their trajectories of at each time step 

    """
    num_dimensions = trajectories.shape[0]
    num_tsteps = trajectories.shape[1]
    num_spins = trajectories.shape[2]
    
    velocities = np.zeros((num_dimensions, num_tsteps - 1, num_spins))

    for dim in range(num_dimensions):
        for step in range(1, num_tsteps):
            velocities[dim, step - 1] = (trajectories[dim, step] - trajectories[dim, step - 1]) / 1.0  # Assuming unit time steps
    
    return velocities

def update_point_pos(fnum, x, y, z, vcolor, scatter_points, disp_fnum):
    """
    Function for repeatedly updating all the elements of the plot at each time step.

    Parameters
    ----------
    fnum : total number of frames in the animation, equal to (Nsteps - 1)
    x : x-cooordinates of spins
    y : y-cooordinates of spins
    z : z-cooordinates of spins
    vcolor: color value representing velocity magnitude of spin
    scatter_points : object containing collection of data points in a scatter plot, displaying the plotted spins
    disp_fnum :  displays the current frame number on the plot 


    Returns
    -------
    scatter_points : Matplotlib's PathCollection object
        object storing collection of points in a scatter plot with updated spin trajectories and velocity colors
    disp_fnum : int
        current frame number (time step) of the plot 

    """
    disp_fnum.set_text('T={:d}'.format(fnum))
    new_x = x[fnum]
    new_y = y[fnum]
    new_z = z[fnum]
    new_cols = vcolor[fnum]
    scatter_points._offsets3d = (new_x, new_y, new_z)
    scatter_points._facecolors = new_cols
    return scatter_points, disp_fnum



def spin_animation_setup(array):
    """
    Extract data from trajectories to return a separate array for each coordinate

    Parameters
    ----------
    path_to_npy : str
        Path to the npy file storing trajectos

    Returns
    -------
    x : numpy.array
        Array containing all the x trajectories with shape (Nsteps, Nspins)
    y : numpy.array
        Array containing all the y trajectories with shape (Nsteps, Nspins)
    z : numpy.array
        Array containing all the z trajectories with shape (Nsteps, Nspins)
    Nsteps : int
        Number of steps of the simulation

    """
    Nsteps = array.shape[1]
    Nspins = array.shape[2]
    x = []
    y = []
    z = []
    for i in range(Nsteps):
        x.append(array[0,i,:])
        y.append(array[1,i,:])
        z.append(array[2,i,:])
    x = np.array(x) # e.g. x[0] trajectories of all spins for timestep 0
    y = np.array(y)
    z = np.array(z)
    
    return x, y, z, Nsteps, Nspins
   


def spin_animation(rpos,outvideo=None):
    """
    Animation function for arrays with the structure:
    3 (xyz) x Nsteps x Nspins

    spin_animation(rpos)    
    spin_animation(rpos,outvideo=<STRING>)

    MANDATORY INPUT PARAMETERS
    ----------
    rpos :     Numpy array
               storing the trajectories of Nspins spins
               for Nsteps time steps. It must be an 
               array of size
               3 x Nsteps x Nspins where
               rpos[0,:,:] stores the x-components
               rpos[1,:,:] stores the y-components
               rpos[2,:,:] stores the z-components
               
               For example, rpos[0,t,s] stores the
               x-component of the trajectory of the
               s-th spin at time step t-th 
    
    OPTIONAL INPUT PARAMETERS
    ----------    
    outvideo:  string indicating the path of a .mp4 
               where the animation will be saved
               as a video (default: None - no video saved)
        
    RETURNS
    -------
    None.

    """

    trajectories = rpos                         # load numpy array storing the trajectories of the spins with size 3 x Ntsteps x Nspins,
                                                # where trajectories[0,t,s] stores the x-component of the trajectory of the s-th spin at time step t-th 
    velocities = calculate_velocity(trajectories)
    tot_velocities = np.linalg.norm(velocities,axis=0) # compute the magnitude of all velocity vectors across all components (x, y, and z).
    np.save("spin_velocities.npy", tot_velocities)
    np.save("spin_trajectories.npy", trajectories)
    
    x, y, z, Nsteps, Nspins = spin_animation_setup(trajectories) 

    # Initalize scatter plot with spin positions color coded by velocity magnitude
    fig = plt.figure(figsize=(16, 9), dpi=1920 / 16)
    ax = fig.add_subplot(111, projection='3d')
    cmap = matplotlib.colormaps['viridis']
    vcolor = cmap(tot_velocities*20000) # Scale values to highlight differences in velocity magnitude
    scatter_points = ax.scatter(x[0], y[0], z[0], 'o', facecolors = vcolor[0], s=10, alpha=1) # Initial spin positions and velocities
    x = x[1:]
    y = y[1:]
    z = z[1:]
    disp_fnum = fig.suptitle('') # Initialize title at top of figure

    #   Create animation by repeatedly running function for a number of frames
    ani = animation.FuncAnimation(
        fig, update_point_pos, frames= Nsteps - 1, fargs=(x, y, z, vcolor, scatter_points, disp_fnum)) 
    ax.set_xlabel("x-position [mm]")
    ax.set_ylabel("y-position [mm]")
    ax.set_zlabel('z-position [mm]')
    plt.show()
    # Save a video with the animation if required
    if(outvideo is not None):
        ani.save(outvideo)

    
