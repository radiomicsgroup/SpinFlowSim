This is the manual of the SpiFlowSim _visu_ module, which can be found in the [visu.py](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/visu.py) file. 

```
Help on module visu:

NAME
    visu

FUNCTIONS
    spin_animation(rpos, outvideo=None)
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
        outvideo:  string indicating the path of an .mp4 
                   where the animation will be saved
                   as a video (default: None - no video saved)
            
        RETURNS
        -------
        None.
```  
