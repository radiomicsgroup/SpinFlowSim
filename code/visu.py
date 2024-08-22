import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib


def _ptupdate(num,x,y,z,graph,txt):
    txt.set_text('T={:d}'.format(num))
    new_x = x[num]
    new_y = y[num]
    new_z = z[num]
    # update properties
    graph._offsets3d = (new_x,new_y,new_z)
    # return modified artists
    return graph,txt


def _spin_animation_setup(array):
    """
    Sets up the npy array for animation

    Parameters
    ----------
    path_to_npy : str
        Path to the npy file

    Returns
    -------
    x : numpy.array
        Array containing all the x positions with shape (Nsteps, Nspins)
    y : numpy.array
        Array containing all the y positions with shape (Nsteps, Nspins)
    z : numpy.array
        Array containing all the z positions with shape (Nsteps, Nspins)
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
    x = np.array(x) # e.g. x[0] positions of all spins for timestep 0
    y = np.array(y)
    z = np.array(z)
    
    return x,y,z,Nsteps,Nspins

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
    outvideo:  string indicating the path of a GIF 
               where the animation will be saved
               as a video (default: None - no video saved)
        
    RETURNS
    -------
    None.

    """
    
    x,y,z,Nsteps,Nspins = _spin_animation_setup(rpos)
    
    # Set up the figure
    fig = plt.figure(figsize=(16, 9), dpi=1920/16)
    ax = fig.add_subplot(111, projection='3d')
    # Colors
    cmap = matplotlib.cm.get_cmap('viridis')
    np.random.seed(123) # To fix desired colors
    rnums = np.random.randn(Nspins)
    cols = cmap(rnums)
    cols = cols[:,:3]
    # Graph
    graph = ax.scatter(x[0], y[0], z[0], 'o',c = cols,s=20, alpha = 1)
    x = x[1:]
    y = y[1:]
    z = z[1:]
    txt = fig.suptitle('')
    
    # Make an animation
    ani=animation.FuncAnimation(fig, _ptupdate, frames=Nsteps-1, fargs=(x,y,z,graph,txt))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_zlabel('z')
    plt.show()
    
    # Save a video with the animation if required
    if(outvideo is not None):
    	ani.save(outvideo)
     
     


