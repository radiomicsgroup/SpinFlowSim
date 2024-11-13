import numpy as np
import time
import graph_tool as gt
from graph_tool import draw as gdraw
from graph_tool import topology as gtop
import lcapy
import PySpice.Logging.Logging as Logging
logger = Logging.setup_logging()
from PySpice.Spice.Netlist import Circuit
from PySpice.Unit import *

class pipenet():
    ''' Class to define a vascular network and use it to simulate IVIM diffusion MRI signals

        obj = pipenet(<parameters>) initialises a vascular network made of connected pipes

        MANDATORY INPUT PARAMETERS
        * nodes:     matrix of node positions (in mm), of size 3 x Nnodes (rows: xpos, ypos, zpos in mm);
                     a node is defined as the input (or as the output) of a pipe
        * radii:     matrix of pipe radii (in mm), of size Nnodes x Nnodes (element (i,j) stores the radius
                     of the pipe connecting node i with node j)
        * qin:       input volumetric flow rate in mm3/s
        * idxin:     index of the input node in the node matrix (nodes[:,idxin] provides the x,y,z coordinates
                     of such a input node)
        * idxout:    index of the output node in the node matrix (nodes[:,idxout] provides the x,y,z coordinates
                     of such an input node)

        OPTIONAL INPUT PARAMETERS
        * radiusin:  radius of the pipe bringing fluid to the input node in mm (default: None).
                     If None, the this will be calculated as the mean of all the
                     radii for all pipes that emanate from the input node
        * visc:      ideal dinaymic viscosity of the fluid in mPa x s (default: 1.20 mPa x s).
                     If blood is modelled, this value should be set to:
                    * 4 x the viscosity of pure blood plasma (when no RBCs are present), for flowmodel
                      "ideal" (the viscosity of pure blood plasma is approx. 1.20 mPa x s, so
                      4 x 1.20 = 4.80 mPa x s)
                    * exactly the viscosity of pure blood plasma (when no RBCs are present;
                      approx. 1.20 mPa x s) for model "blinder"
        * flowmodel: string indicating the model that should be used to calculate flow resistance
                     in each pipe/segment of the vascular network. Values could be:
                     - "ideal":   to use the Hagen-Poiseuille law for ideal Netwonian fluids;
                     - "blinder": to use the modified Hagen-Poiseuille law as in Blinder et al,
                                  Nat Neurosci 2013, 16(7): 889-897, doi: 10.1038/nn.3426. In this model,
                                  the viscosity of the blood is a function of the vessel radius; for very
                                  large vessels, the viscosity approaches 4 x the viscosity of pure blood
                                  plasma when no red blood cells (RBCs) are present;
                     Default value is "blinder"
        * solver:   string indicating the library package that should be used to calculate the flow matrix
                    of the vascular network. Values could be: 
                     - "numerical": to use the PySpice package;
                     - "symbolic": to use the lcapy symbolic computation package;
                     Default value is "numerical"

        RETURNS
        * obj:       instance of class pipenet, storing an initialised vascular network.
                     It features the following attributes:
                     - flowmodel: a string indicating the hemodynamics model (e.g., "ideal", "blinder").
                                  It is a copy of the input parameter "flowmodel"
                     - nodepos: a 3 x Nnodes matrix storing the positions of the nodes (in mm); a node is defined
                                as the input (or as the output) of a pipe. It is a copy of the
                                input parameter "nodes"
                                nodepos[0,nn] stores the x-component of the nn-th node (units: mm),
                                nodepos[1,nn] stores the y-component of the nn-th node (units: mm),
                                nodepos[2,nn] stores the z-component of the nn-th node (units: mm).
                     - totnodes: total number of nodes (Nnodes)
                     - nodein:  index of the input node (indices start from 0). It stores the value passed
                                with input parameter "idxin"
                     - nodeout: index of the output node (indices start from 0). It stores the value passed
                                with input parameter "idxout"
                     - connmat: connectivity matrix C = [cij] of size Nnodes x Nnodes; element C[i,j] is +1/-1
                                if two nodes are connected, 0 otherwise, with NaNs being stored along the diagonal;
                                cij = 1 --> flow from i to j; cij = -1 --> flow from j to i;
                                cij = -cji
                     - flowmat: volumetric flow rate matrix Q = [qij] of size Nnodes x Nnodes; element Q[i,j] stores
                                the volumetric flow rate qij between node i and node j (units: mm3/s);
                                NaNs are stored along the diagonal;
                                qij > 0 --> flow from i to j; qij = < 0 --> flow from j to i;
                                qij = -qji
                     - velmat:  flow velocity matrix V = [vij] of size Nnodes x Nnodes; element V[i,j] stores
                                the velocity of the flow vij between node i and node j (units: mm/s);
                                NaNs are stored along the diagonal;
                                vij > 0 --> flow from i to j; vij = < 0 --> flow from j to i;
                                vij = -vji
                     - qin:     input volumetric flow rate (units: mm3/s) arriving at the input node (see nodein)
                     - radiusin: radius of the pipe bringing fluid to the input node (units: mm). It is a copy
                                of input parameter "radiusin" if it was provided; if not, radiusin is set as the
                                mean of the radii of all branches emanating from the input node
                     - radmat:  matrix of pipe radii R = [rij] of size Nnodes x Nnodes. rij stores the radius
                                of the pipe connecting node i with node j, s.t. rij = rji (units: mm).
                                NaNs are stored along the diagonal. It is essentially a copy of the input
                                parameter "radii", with NaNs along the diagonal
                    - lengthmat: matrix of pipe lengths L = [lij] of size Nnodes x Nnodes. lij stores the length
                                 radius of the pipe connecting node i with node j, s.t. lij = lji (units: mm).
                                 NaNs are stored along the diagonal
                    - resmat:   matrix of pipe resistances Z = [zij] of size Nnodes x Nnodes. zij stores the
                                resistance experienced by the flow between nodes i and j, s.t. zij = zji
                                (units: mPa x s/mm3). The way resmat is computed depends on the adopted
                                hemodynamics model (see flowmodel above)
                    - iopaths:  list of P elements containing all possible paths from the input node to the output
                                node. Each element iopaths[p] contains a description for the p-th path, with
                                p = 0, ..., P - 1, in the form of a list storing the sequence of nodes that
                                are encountered by a molecule following path p-th. For example, assume iopaths
                                stores something like [[1,2], [1,0,2]]; this would imply that there are two possible
                                paths between the input node 1 and the output node 2, and these would be [1,2]
                                (path number 0) and [1,0,2] (path number 1)
                    - gnet:     a directed Graph object from the graph_tool module populated with vertices ('node_list')
                                and edges ('edge_list') from the connectivity matrix C.
                    - visco:    dynamic viscosity of the fluid (units: mPa x s). It is a copy of input parameter "visc"
                    - solver:   a string indicating the computation package (e.g., "numerical", "symbolic").
                                  It is a copy of the input parameter "solver"

        DEPENDENCIES
        To use pipenet you need:
        * numpy (developed with numpy 1.23.5; https://numpy.org)
        * lcapy (developed with lcapy 1.10; https://lcapy.readthedocs.io/en/latest)
        * PySpice (developed with PySpice 1.5; https://github.com/FabriceSalvaire/PySpice)
        * graph-tool (developed with graphtool 2.45, commit b1a649d8; https://graph-tool.skewed.de)
        * time (python standard library)
        pipenet was developed with python 3.10.8
    '''
    def __init__(self,nodes,radii,qin,idxin,idxout,radiusin=None,visc=1.20,flowmodel='blinder',solver='numerical'):
        '''
        pipenet object constructor -- type help(syn.pipenet) for the help manual
        '''

        # Check for consistency of input parameters
        if( nodes.shape[0]!=3):
            raise RuntimeError('ERROR. The matrix of nodes must have size 3 x number of nodes.')
        Nnodes = nodes.shape[1]

        if( (radii.shape[0]!=radii.shape[1]) or (radii.shape[0]!=Nnodes) ):
            raise RuntimeError('ERROR. The size of the matrix of pipe radii does not match the number of nodes.')

        for ii in range(0,Nnodes):
            radii[ii,ii] = np.nan           # Make sure the radius matrix stores NaN along the diagonal

        if( (idxin > (Nnodes-1)) or (idxout > (Nnodes-1)) or (idxin < 0) or (idxout < 0) ):
            raise RuntimeError('ERROR. Invalid index of the input and/or output output nodes.')

        if( (flowmodel=='ideal') or (flowmodel=='blinder') ):
            pass
        else:
            raise RuntimeError('ERROR. The flow model {} does not exist - use one of ideal or blinder'.format(flowmodel))
        
        if( (solver=='numerical') or (solver=='symbolic') ):
            pass
        else:
            raise RuntimeError('ERROR. The solver {} does not exist - use one of numerical or symbolic'.format(solver))

        # Calculate matrix of pipe lengths
        lengths = np.zeros((Nnodes,Nnodes))
        for ii in range(0,Nnodes):
            for jj in range(0,Nnodes):
                if radii[ii,jj]!=0:
                    lengths[ii,jj] = np.sqrt( (nodes[0,ii] - nodes[0,jj])**2 + (nodes[1,ii] - nodes[1,jj])**2 + (nodes[2,ii] - nodes[2,jj])**2 )
        for ii in range(0,Nnodes):
            lengths[ii,ii] = np.nan

        # Store vascular network parameters
        self.flowmodel = flowmodel            # Hemodynamics model
        self.nodepos = np.copy(nodes)         # Matrix of node positions (in mm), of size 3 x Nnodes (rows: xpos, ypos, zpos in mm)
        self.totnodes = Nnodes                # Number of nodes in the vascular network
        self.radmat = np.copy(radii)          # Matrix of pipe radii (in mm), of size Nnodes x Nnodes (element (i,j) stores the radius of the pipe connecting node i with node j)
        self.nodein = idxin                   # Index of input node
        self.nodeout = idxout                 # Index of output node
        self.visco = visc                     # Dynamic viscosity of the fluid in mPa x s (if blood is modelled, this is the viscosity of pure plasma with no red blood cells)
        self.qin = qin                        # Input volumetric flow rate in mm3/s
        self.lengthmat = np.copy(lengths)     # Matrix of pipe lengths (in mm), of size Nnodes x Nnodes (element (i,j) stores the length of the pipe connecting node i with node j)
        self.solver= solver                   # Computation package
        # Calculate and store the radius of the input node if needed
        if(radiusin is None):
            nonzero_count =np.count_nonzero(~np.isnan(self.radmat[self.nodein, :]) & (self.radmat[self.nodein, :] != 0))
            radius_sum = np.nansum(self.radmat[self.nodein,:])
            radiusin = radius_sum / nonzero_count  
        self.radiusin = radiusin

        # Calculate and store resistances, flows and velocities
        if(self.flowmodel=='ideal'):

            # Resistances in mPa x s/mm3
            resmat = ComputeResHagPois(self.radmat, self.lengthmat, self.visco)
            self.resmat = np.copy(resmat)    # Matrix of resistances (in mPa x s/mm3), of size Nnodes x Nnodes (element (i,j) stores resitance between i and j)

            # Volumetric flow in mm3/s
            flowmat = ComputeFlow(self.radmat, self.resmat, self.nodein, self.nodeout, self.qin, self.solver)
            self.flowmat = np.copy(flowmat)  # Matrix of flows (in mm3/s), of size Nnodes x Nnodes (element (i,j) stores flow qij between i and  j, s.t., qij = - qji)

            # Velocities in mm/s
            areas = np.pi*self.radmat*self.radmat
            velmat = flowmat/areas
            self.velmat = np.copy(velmat)    # Matrix of velocities (in mm/s), of size Nnodes x Nnodes (element (i,j) stores velocity vij between i and  j, s.t., vij = - vji)

        elif(self.flowmodel=='blinder'):

            # Resistances in mPa x s/mm3
            resmat = ComputeResModHagPois(self.radmat, self.lengthmat, self.visco)
            self.resmat = np.copy(resmat)    # Matrix of resistances (in mPa x s/mm3), of size Nnodes x Nnodes (element (i,j) stores resitance between i and j)

            # Volumetric flow in mm3/s
            flowmat = ComputeFlow(self.radmat, self.resmat, self.nodein, self.nodeout, self.qin,self.solver)
            self.flowmat = np.copy(flowmat)  # Matrix of flows (in mm3/s), of size Nnodes x Nnodes (element (i,j) stores flow qij between i and  j, s.t., qij = - qji)

            # Velocities in mm/s
            areas = np.pi*self.radmat*self.radmat
            velmat = flowmat/areas
            self.velmat = np.copy(velmat)    # Matrix of velocities (in mm/s), of size Nnodes x Nnodes (element (i,j) stores velocity vij between i and  j, s.t., vij = - vji)

        else:
            raise RuntimeError('ERROR. The flow model {} does not exist - use one of "ideal" or "blinder"'.format(self.flowmodel))

        # Evaluate and store connectivity matrix based on the flow matrix
        connmat = np.copy(self.flowmat)  # Connectivity matrix
        connmat[connmat>0] = 1
        connmat[connmat<0] = -1
        self.connmat = np.copy(connmat)  # Connectivity matrix size Nnodes x Nnodes (element (i,j) is +1/-1 if two nodes are connected, 0 otherwise), s.t.
                                         # i) connmat[i,j] = - connmat[j,i];
                                         # ii) connmat[i,j] = 1 --> flow from i to j; connmat[i,j] = -1 --> flow from j to i

        # Calculate and store all possible paths from input to output node. Construct directed graph with vertices and edges from connectivity matrix.
        iopaths, gnet = ComputeAllPaths(self.connmat, self.nodein, self.nodeout)
        self.iopaths = iopaths
        self.gnet = gnet



    def dMRISynMea(self,bval,Gdur,Gsep,Gdir,Nspins=2000,dt=2e-5,Gx=None,Gy=None,Gz=None,rndseed=20181019,bcon='periodic'):
        '''
        Synthesises a diffusion-weighted MRI measurement from protons flowing within a vascular network

        magtot, phasetot, phasespin, rmat, tarray, grx, gry, grz = obj.dMRISynMea(<parameters>)

        MANDATORY INPUT PARAMETERS
        * bval: b-value for a PGSE sequence (scalar; units: s/mm2);
        * Gdur: gradient duration (small delta) for a PGSE sequence (scalar; units: s); if optional 
                Gx, Gy, Gz gradient waveforms are provided, this parameter is ignored;
        * Gsep: gradient separation (large Delta) for a PGSE sequence (scalar; units: s); if optional
                Gx, Gy, Gz gradient waveforms are provided, this parameter is ignored
        * Gdir: gradient direction for a PGSE sequence; it must be a 1D numpy array of size 3, storing the
                x, y and z components of the gradient direction; ; if optional
                Gx, Gy, Gz gradient waveforms are provided, this parameter is ignored
                
        **** NOTE THAT ALL MANDATORY INPUT PARAMETERS ABOUT b-value, GRADIENT DURATION/SEPARATION/DIRECTION
        ARE IGNORED IF CUSTOM GRADIENT WAVEFORMS ARE PROVIDED VIA OPTIONAL INPUT PARAMETERS Gx, Gy, Gz                

        OPTIONAL INPUT PARAMETERS
        * Nspins: number of spins to simulate (default: 2000). Note that in practice a slighter smaller
                  number of particles will be simulated, as some particles will leave the vascular network
                  during the simulation and will not be used to synthesise the MRI signal, as in
                  Phi Van V et al, NMR in Biomedicine 2021, 34(7):e4528, doi: 10.1002/nbm.4528
        * dt:     temporal resolution of the simulation (scalar; units: s; default: 2e-5, i.e., 20 us)
        * Gx:     waveform for x-component of the diffusion encoding gradient, expressed in T/mm.
                  If specified, previous parameters bval, Gdur, Gsep, Gdir will be ignored.
                  Gx must be a 1D array storing the values of Gx over time (Gx(t=0), Gx(t=dt), Gx(t=2*dt), ....).
                  If Gx is specified, Gy and Gz must be specified too. Gx, Gy and Gz must have the same length.
                  Default value for Gx: None
        * Gy:     waveform for y-component of the diffusion encoding gradient, expressed in T/mm.
                  If specified, previous parameters bval, Gdur, Gsep, Gdir will be ignored.
                  Gy must be a 1D array storing the values of Gy over time (Gy(t=0), Gy(t=dt), Gy(t=2*dt), ....).
                  If Gy is specified, Gx and Gz must be specified too. Gx, Gy and Gz must have the same length.
                  Default value for Gy: None
        * Gz:     waveform for z-component of the diffusion encoding gradient, expressed in T/mm.
                  If specified, previous parameters bval, Gdur, Gsep, Gdir will be ignored.
                  Gz must be a 1D array storing the values of Gz over time (Gz(t=0), Gz(t=dt), Gz(t=2*dt), ....).
                  If Gz is specified, Gx and Gy must be specified too. Gx, Gy and Gz must have the same length.
                  Default value for Gz: None
        *rndseed: Random seed used for spin position initialisation in the network. Spin initialisation is
                  performed by sampling at uniform spatial locations consecutive waves of spins that flow
                  from the input node to the output node, making sure that the number of spins in each pipe at
                  t = 0 is proportional to the flow in that pipe (conservation of mass).
                  Default: 20181019
        * bcon:   string specifying the boundary condition for particles that reach the output node.
                  It can be set to "periodic" or "feedback". Default: "periodic".
                  - "periodic" -> a particle that reaches the output node starts again
                    its trajectory, but starting from the output node, rather than from the input
                    node. In practice this makes the vascular network periodic, and a "copy" of the
                    network is "attached" to the output node. A particle reaching the output node
                    will continue travelling through this "copy" of the network.
                  - "feedback" -> a spin that reaches the output node is fed back to the input node,
                    and it will follow again the same trajectory it had taken before. Note that
                    in this case the particle positions and velocities over time will exhibit a
                    discontinuity when the "jump" from output node to input node occurs. This could cause
                    some unexpected signal behaviours: for example, the signal from a
                    straight pipe (where all spins travel in parallel at the same velocity) will show
                    magnitude MRI signal attenuation, while in theory should show no signal
                    attenuation but only a phase shift. Also, the signal decay could be overstimated when
                    the diffusion encoding gradient has a strong component along the direction of
                    such a "jump", and a considerable number of particles "jump" during the simulation.

        RETURNS
        * magtot:    magnitude of the synthesised diffusion MRI signal (note that in absence of diffusion-weighting,
                     a magnitude of 1.0 will be observed)
        * phasetot:  phase of the synthesised diffusion MRI signal, defined in the range [-pi; pi]
        * phasespin: 1D array storing the phase accrual of all individual spins. The array has size Nspins,
                     so that phasespin[s] stores the phase accrual of the s-th spin, for s = 0, ..., Nspins - 1.
                     It ranges between ]-inf; +inf[ so that one can study the number of phase wraps caused
                     by diffusion-weighting
        * rmat:     matrix storing the trajectories of the spins (units: mm). It has size
                    3 x Ntsteps x Nspins , where Nspins is the number of requested spins,
                    and Ntsetps is the number of simulation time points at a temporal resolution equal to dt.
                    Note that for conservation of total mass the vascular network a boundary condition is
                    applied depending on the value of the input "bcon". This means that spins that reach
                    the output node before the simulation is finshed, can either i) start again their
                    trajectory in a "copy" of the network attached to the output node
                    (if bcon is set to "periodic"), or ii) be fed back to the input node
                    (if bcon is set to "feedback). A detailed description of rmat is
                    - rmat[0,k,s] -> x-component of the trajectory of the s-th spin, at time step k-th
                    - rmat[1,k,s] -> y-component of the trajectory of the s-th spin, at time step k-th
                    - rmat[2,k,s] -> z-component of the trajectory of the s-th spin, at time step k-th
        * tarray:   time array (units: s) of length Ntsteps. tarray[k] = k*dt stores the time corresponding to
                    the k-th simulation step
        * grx:      waveform for the x-component of the diffusion encoding gradient (Gx), expressed in T/mm. It has size
                    Ntsteps, and grx[k] is the value of Gx at time instant tarray[k] = k*dt. If optional input Gx was
                    provided, than grx is just a copy of Gx
        * gry:      waveform for the y-component of the diffusion encoding gradient (Gy), expressed in T/mm. It has size
                    Ntsteps, and gry[k] is the value of Gy at time instant tarray[k] = k*dt. If optional input Gy was
                    provided, than gry is just a copy of Gy
        * grz:      waveform for the z-component of the diffusion encoding gradient (Gz), expressed in T/mm. It has size
                    Ntsteps, and grz[k] is the value of Gz at time instant tarray[k] = k*dt. If optional input Gz was
                    provided, than grz is just a copy of Gz
        '''
        ### Deal with inputs
        if( not np.isscalar(Nspins) ):
            raise RuntimeError('The number of molecules must be scalar integer > 0')
        if( Nspins < 0 ):
            raise RuntimeError('The number of molecules must be scalar integer > 0')
        if( not isinstance(Nspins, int) ):
            raise RuntimeError('ERROR. The number of molecules must be a scalar integer > 0')
        if( not np.isscalar(dt) ):
            raise RuntimeError('The temporal resolution of the simulation must be a scalar > 0 expressed in seconds')
        if( dt < 0.0 ):
            raise RuntimeError('The temporal resolution of the simulation must be a scalar > 0 expressed in seconds')

        # Check whether the user provided an actual gradient waveforms or whether they are requesting a simple PGSE
        if( (Gx is not None) or (Gy is not None) or (Gz is not None) ):

            # If the user passed Gx(t), then Gy(t) and Gz(t) are also required compulsory - not that they are in T/mm, e.g., 65 mT/m = 65*1e-6 T/mm
            if(Gx is None):
                raise RuntimeError('ERROR. All 3 components of the diffusion encoding gradient (x,y,z) must be specified. At least one of the three is missing.')
            if(Gy is None):
                raise RuntimeError('ERROR. All 3 components of the diffusion encoding gradient (x,y,z) must be specified. At least one of the three is missing.')
            if(Gz is None):
                raise RuntimeError('ERROR. All 3 components of the diffusion encoding gradient (x,y,z) must be specified. At least one of the three is missing.')

            # Make sure Gx(t), Gy(t) and Gz(t) are 1D arrays
            if( (Gx.ndim!=1) or (Gy.ndim!=1) or (Gz.ndim!=1) ):
                raise RuntimeError('ERROR. All of Gx, Gy, Gz must be 1D numpy arrays storing the temporal evolution of Gx, Gy, Gz for one MRI measurement')

            # Make sure Gx(t), Gy(t) and Gz(t) have the same length
            if( (Gx.size!=Gy.size) or (Gx.size!=Gz.size) ):
                raise RuntimeError('ERROR. All of Gx, Gy, Gz must have the same duration')

            # Get information on the duration of the simulation and define gamma
            Ntpoints = Gx.size                  # Number of time points of the simulation
            tarray = dt*np.linspace(0,float(Ntpoints)-1,Ntpoints)     # Time array (units: s)

        # The user did not provide the waveforms Gx(t), Gy(t), Gz(t) - we've got to generate them
        else:

            # Make sure the user passed a scalar for bval, Gdur, Gsep
            if( not np.isscalar(bval) ):
                raise RuntimeError('ERROR. The b-value must be a scalar (units: s/mm2), but you passed {} s/mm2'.format(bval))
            if( not np.isscalar(Gdur) ):
                raise RuntimeError('ERROR. The gradient duration must be a scalar (units: s), but you passed {} s'.format(Gdur))
            if( not np.isscalar(Gsep) ):
                raise RuntimeError('ERROR. The gradient separation must be a scalar (units: s), but you passed {} s'.format(Gsep))
            if(bval<0):
                raise RuntimeError('ERROR. A negative b-value does not make sense, but you passed {} s/mm2'.format(bval))
            if(Gdur<0):
                raise RuntimeError('ERROR. A negative gradient duration does not make sense, but you passed {} s'.format(bval))
            if(Gsep<0):
                raise RuntimeError('ERROR. A negative gradient separation does not make sense, but you passed {} s'.format(bval))
            if(Gsep<Gdur):
                raise RuntimeError('ERROR. The gradient separation cannot be shorter than the duration - you passed a duration of {} and a separation of {} s'.format(Gdur,Gsep))

            # Make sure the user passed a 1D array of size 3 for the gradient direction
            if( (Gdir.ndim!=1) or (Gdir.size!=3) ):
                raise RuntimeError('ERROR. The gradient direction must be a 1D numpy array of size 3, but you passed {}'.format(Gdir))

            # Get gradient waveforms
            if( (bval==0.0) or (Gdur==0.0) or (Gsep==0.0) or (np.linalg.norm(Gdir)==0.0) ):

                # Deal with a non-DW measurement
                Gsep = 0.0
                Gdur = 0.0
                bval = 0.0
                Ntpoints = int(np.round(0.040/dt))                         # There is no gradient - just simulate a standard TE of 40 ms
                tarray = dt*np.linspace(0,float(Ntpoints)-1,Ntpoints)      # Time array (units: s)
                Gx = np.zeros((Ntpoints))
                Gy = np.zeros((Ntpoints))
                Gz = np.zeros((Ntpoints))

            else:

                Gx, Gy, Gz,tarray = getGradPGSE(bval, Gdur, Gsep, Gdir, dt)  # Get PGSE gradient waveforms
                Ntpoints = tarray.size

        ### Synthesise the MRI signal given the trajectories of molecules that have been seeded uniformly in the vascular network

        # Get trajectories given the current vascular network
        rmat_final = self.GetTrajUniformSeed(Ntpoints,Nspins,dt,seednumber=rndseed,boundary=bcon)  # Flow of Nspins for Ntpoints points over time; temporal/spatial resolution dt/dr [ms]/[mm]

        # Get MRI signal and phase information
        magtot, phasetot, phasespin = Traj2Signal(rmat_final,Gx,Gy,Gz,dt)

        # Return, in this order: total signal, total phase, each individual spin phase, each individual spin trajectory, t array, Gx(t) array, Gy(t) array, Gz(t) array
        return magtot, phasetot, phasespin, rmat_final, tarray, Gx, Gy, Gz




    def dMRISynProt(self,bvarr,Gdurarr,Gseparr,Gdirarr,Nwater=2000,deltat=2e-5,Gradx=None,Grady=None,Gradz=None,Nrep=1,seed=20181019,boundcon='periodic',vb=True):
        '''
        Synthesises a diffusion-weighted MRI protocol made of multiple measurements from
        protons flowing within a vascular network

        mag, phase = obj.dMRISynProt(<parameters>)

        MANDATORY INPUT PARAMETERS
        * bvarr: array of b-values for a PGSE sequence (units: s/mm2; it must be a 1D array of size Nmeas)
        * Gdurarr: array of gradient durations (small delta) for a PGSE sequence (1D array of size Nmeas; units: s)
        * Gseparr: array of gradient separations (large Delta) for a PGSE sequence (1D array of size Nmeas; units: s)
        * Gdirarr: array of gradient directions for a PGSE sequence; it must be a 2D numpy array of size Nmeasx3,
                   where Gdirarr[m,0] stores the x-component of the gradient direction corresponding to the m-th
                   measurement; Gdirarr[m,1] the y-component; Gdirarr[m,2] the z-component

        ***** NOTE THAT ALL INPUT PARAMETERS ABOUT b-value, GRADIENT DIRECTION/DURATION/SEPARATION ARE IGNORED
              IF CUSTOM GRADIENT WAVEFORMS ARE PROVIDED WITH OPTION INPUT PARAMETERS Gradx, Grady, Gradz
              (SEE BELOW) 

        OPTIONAL INPUT PARAMETERS
        * Nwater: number of spins to simulate (default: 2000). Note that in practice a slighter smaller
                  number of particles will be simulated, as some particles will leave the vascular network
                  during the simulation and will not be used to synthesise the MRI signal, as in
                  Phi Van V et al, NMR in Biomedicine 2021, 34(7):e4528, doi: 10.1002/nbm.4528
        * deltat: temporal resolution of the simulation (scalar; units: s; default: 2e-5, i.e., 20 us)
        * Gradx:  waveforms for x-components (Gx) of the diffusion encoding gradients, expressed in T/mm.
                  If specified, previous parameters bvarr, Gdurarr, Gseparr, Gdirarr will be ignored.
                  Gradx must be a 2D array storing the values of Gx over time (Gx(t=0), Gx(t=dt), Gx(t=2*dt), ....).
                  for all measurements m = 0, ..., Nmeas - 1.
                  Gx(m,k) stores Gx for the m-th protocol measurement at the k-th time step.
                  If Gradx is specified, Grady and Gradz must be specified too.
                  Gradx, Grady and Gradz must have the same size. Default value for Gradx: None
        * Grady:  waveforms for y-components (Gy) of the diffusion encoding gradients, expressed in T/mm.
                  If specified, previous parameters bvarr, Gdurarr, Gseparr, Gdirarr will be ignored.
                  Grady must be a 2D array storing the values of Gy over time (Gy(t=0), Gy(t=dt), Gy(t=2*dt), ....).
                  for all measurements m = 0, ..., Nmeas - 1.
                  Gy(m,k) stores Gy for the m-th protocol measurement at the k-th time step.
                  If Grady is specified, Gradx and Gradz must be specified too.
                  Gradx, Grady and Gradz must have the same size. Default value for Grady: None
        * Gradz:  waveforms for z-components (Gz) of the diffusion encoding gradients, expressed in T/mm.
                  If specified, previous parameters bvarr, Gdurarr, Gseparr, Gdirarr will be ignored.
                  Gradz must be a 2D array storing the values of Gz over time (Gz(t=0), Gz(t=dt), Gz(t=2*dt), ....).
                  for all measurements m = 0, ..., Nmeas - 1.
                  Gz(m,k) stores Gz for the m-th protocol measurement at the k-th time step.
                  If Gradz is specified, Gradx and Grady must be specified too.
                  Gradx, Grady and Gradz must have the same size. Default value for Gradz: None
        * Nrep:   number of repetitions of the simulation. For each repetition, a different seed number is used,
                  resulting in slightly different initial spin positions, and hence MRI signals.
                  Note that the calculation of the spin trajectories (the most costly computational step) is
                  done once per each simulation repetition, and then the calculated trajectories are used
                  to synthesise signals for all the diffusion MRI measurements of the protocol.
        * seed:   Random seed used for spin position initialisation in the network. Spin initialisation is
                  performed by sampling at uniform spatial locations consecutive waves of spins that flow
                  from the input node to the output node, making sure that the number of spins in each pipe at
                  t = 0 is proportional to the flow in that pipe (conservation of mass).
                  Default: 20181019
        * boundcon: string specifying the boundary condition for particles that reach the output node.
                  It can be set to "periodic" or "feedback". Default: "periodic".
                  - "periodic" -> a particle that reaches the output node starts again
                    its trajectory, but starting from the output node, rather than from the input
                    node. In practice this makes the vascular network periodic, and a "copy" of the
                    network is "attached" to the output node. A particle reaching the output node
                    will continue travelling through this "copy" of the network.
                  - "feedback" -> a spin that reaches the output node is fed back to the input node,
                    and it will follow again the same trajectory it had taken before. Note that
                    in this case the particle positions and velocities over time will exhibit a
                    discontinuity when the "jump" from output node to input node occurs. This could cause
                    some unexpected signal behaviours: for example, the signal from a
                    straight pipe (where all spins travel in parallel at the same velocity) will show
                    magnitude MRI signal attenuation, while in theory should show no signal
                    attenuation but only a phase shift. Also, the signal decay could be overstimated when
                    the diffusion encoding gradient has a strong component along the direction of
                    such a "jump", and a considerable number of particles "jump" during the simulation.
        * vb:       verbose. Set this parameter to True (or 1) if you want feedback on the simulation
                    duration to be printed on the standard output. Any other value will be trated as
                    a False and no feedback will be printed. Default: True

        RETURNS
        * mag:    1D or 2D array storing the magnitude of the synthesised MRI signals given the input protocol.
                  If Nrep = 1, then mag will be a 1D array of length Nmeas, with mag[m] storing the
                  magnitude for the m-th protocol measurement. If Nrep > 1, then mag will be a 2D numpy array
                  of size Nrep x Nmeas. In that case, mag[r,m] will store the magnitude for the m-th protocol
                  measurement obtained during the r-th simulation repetition (r = 0, ..., Nrep - 1).
        * phase:  1D or 2D array storing the phase of the synthesised MRI signals given the input protocol.
                  If Nrep = 1, then phase will be a 1D array of length Nmeas, with phase[m] storing the
                  phase for the m-th protocol measurement. If Nrep > 1, then phase will be a 2D numpy array
                  of size Nrep x Nmeas. In that case, phase[r,m] will store the phase for the m-th protocol
                  measurement obtained during the r-th simulation repetition (r = 0, ..., Nrep - 1).
                  phase is defined in [-pi; pi].
        '''

        ### Deal with inputs
        if( not np.isscalar(Nwater) ):
            raise RuntimeError('ERROR. The number of molecules must be a scalar integer > 0')
        if( not isinstance(Nwater, int) ):
            raise RuntimeError('ERROR. The number of molecules must be a scalar integer > 0')
        if( Nwater < 0 ):
            raise RuntimeError('ERROR. The number of molecules must be scalar integer > 0')
        if( not np.isscalar(deltat) ):
            raise RuntimeError('ERROR. The temporal resolution of the simulation must be a scalar > 0 expressed in seconds')
        if( deltat < 0.0 ):
            raise RuntimeError('ERROR. The temporal resolution of the simulation must be a scalar > 0 expressed in seconds')
        if( not np.isscalar(Nrep) ):
            raise RuntimeError('ERROR. The number of simulation repetitions must be a scalar - an integer > 0')
        if( Nrep < 0):
            raise RuntimeError('ERROR. The number of simulation repetitions must be a scalar - an integer > 0')
        if( not isinstance(Nrep, int) ):
            raise RuntimeError('ERROR. The number of simulation repetitions must be a scalar - an integer > 0')
        if(not np.isscalar(vb)):
            raise RuntimeError('ERROR. The verbose parameter must be a scalar')

        # Check whether the user provided an actual gradient waveforms or whether they are requesting a simple PGSE
        if( (Gradx is not None) or (Grady is not None) or (Gradz is not None) ):

            # If the user passed Gx(t), then Gy(t) and Gz(t) are also required compulsory - not that they are in T/mm, e.g., 65 mT/m = 65*1e-6 T/mm
            if(Gradx is None):
                raise RuntimeError('ERROR. All 3 components of the diffusion encoding gradients (x,y,z) must be specified. At least one of the three is missing.')
            if(Grady is None):
                raise RuntimeError('ERROR. All 3 components of the diffusion encoding gradients (x,y,z) must be specified. At least one of the three is missing.')
            if(Gradz is None):
                raise RuntimeError('ERROR. All 3 components of the diffusion encoding gradients (x,y,z) must be specified. At least one of the three is missing.')

            # Make sure Gx(t), Gy(t) and Gz(t) are 2D arrays
            if( (Gradx.ndim!=2) or (Grady.ndim!=2) or (Gradz.ndim!=2) ):
                raise RuntimeError('ERROR. All of Gx, Gy, Gz must be 2D numpy arrays of size Nmeas x Nsteps, with Nmeas being the no. of measurements making up the protocol')

            # Make sure Gx(t), Gy(t) and Gz(t) have the same length
            if( (Gradx.shape[0]!=Grady.shape[0]) or (Gradx.shape[0]!=Gradz.shape[0]) or (Gradx.shape[1]!=Grady.shape[1]) or (Gradx.shape[1]!=Gradz.shape[1]) ):
                raise RuntimeError('ERROR. All of Gx, Gy, Gz must store the same number of diffusion MRI measurements and must have the same temporal duration')

            # Get information on the duration of the simulation and define gamma
            Ndmri = Gradx.shape[0]                # Number of measurement making up the diffusion MRI protocol to be simulated


        # The user did not provide the waveforms Gx(t), Gy(t), Gz(t) - we've got to generate PGSE wave forms
        else:

            # Make sure the user passed an array of b-values, Gdur, Gsep, and gradient directions
            if( bvarr.ndim!=1 ):
                raise RuntimeError('ERROR. The b-value must be an 1D array (units: s/mm2), but you passed {}'.format(bvarr))
            if( Gdurarr.ndim!=1  ):
                raise RuntimeError('ERROR. The gradient duration must be 1D array (units: s), but you passed {}'.format(Gdurarr))
            if( Gseparr.ndim!=1  ):
                raise RuntimeError('ERROR. The gradient separation must be a 1D array (units: s), but you passed {}'.format(Gseparr))
            if( Gdirarr.ndim!=2  ):
                raise RuntimeError('ERROR. The gradient direction array must have size Nmeas x 3, but you passed {}'.format(Gdirarr))

            # Make sure the user passed 1D arrays that make sense
            if( (bvarr.size!=Gdurarr.size) ):
                raise RuntimeError('ERROR. The number of measurements in the b-value and gradient duration arrays do not match')
            if( (bvarr.size!=Gseparr.size) ):
                raise RuntimeError('ERROR. The number of measurements in the b-value and gradient separation array do not match')
            if( (bvarr.size!=Gdirarr.shape[0]) ):
                raise RuntimeError('ERROR. The number of measurements in the b-values and gradient direction array do not match')

            # Make sure that the gradient direction has 3 components (one for x, one for y, one for z)
            if( Gdirarr.shape[1]!=3 ):
                raise RuntimeError('ERROR. The gradient direction array must have size Nmeas x 3')

            # Find out the minimum number of time points required to simulate all the dMRI measurements making up the input protocol - if there are only b = 0, just simulate 40 ms
            Ndmri = bvarr.size             # Number of measurements
            tsim = Gseparr + Gdurarr       # Minimum time duration of the simulation required for each diffusion MRI measurements
            tsim[bvarr==0.0] = np.nan
            tsim[Gseparr==0.0] = np.nan
            tsim[Gdurarr==0.0] = np.nan
            try:
                maxlenid = np.nanargmax(tsim)   # Index of the measurement requireing the longest simulation time - we will use that for all measurements
                dummygrads = getGradPGSE(bvarr[maxlenid], Gdurarr[maxlenid], Gseparr[maxlenid], Gdirarr[maxlenid,:], deltat)  # Get quickly the gradient waveform for the dMRI measurement requiring the longest simulation time
                Nstepsim = dummygrads[0].size     # Total number of simulation time steps required
            except:
                print('WARNING: you are simulating a protocol made of b = 0 measurements only. I will run the simulation for 40 ms')
                Nstepsim = int(np.round(0.040/deltat))    # Number of time points required to simulate trajectories for 40 ms

            # Generate all gradient waveforms for all dMRI measurements
            Gradx = np.zeros((Ndmri,Nstepsim))
            Grady = np.zeros((Ndmri,Nstepsim))
            Gradz = np.zeros((Ndmri,Nstepsim))
            for dd in range(0,Ndmri):
                if(not np.isnan(tsim[dd])):
                    # Get waveform for non-zero b-values
                    glist = getGradPGSE(bvarr[dd],Gdurarr[dd],Gseparr[dd],Gdirarr[dd,:],deltat,Nsample=Nstepsim)
                    Gradx[dd,:] = glist[0]
                    Grady[dd,:] = glist[1]
                    Gradz[dd,:] = glist[2]


        ### Synthesise measurements for current protocol

        # Allocate output signal and phase for the different repetitions
        protmag = np.zeros((Nrep,Ndmri))
        protphase = np.zeros((Nrep,Ndmri))

        # There are no repetitions of the simulation required
        if(vb==True):
            print('')
            print('++++++ dMRISynProt: simulation information')
            print('')

        if(Nrep==1):

            # Get spin trajectories once
            if(vb==True):
                t1 = time.time()
            bufferlist = self.dMRISynMea(None,None,None,None,Nspins=Nwater,dt=deltat,Gx=Gradx[0,:],Gy=Grady[0,:],Gz=Gradz[0,:],rndseed=seed,bcon=boundcon)
            flowwalk = bufferlist[3]    # Trajectories over time of all spins
            if(vb==True):
                t2 = time.time()
                print('   *** It took {} s to generate spin walks'.format(t2-t1))

            # Loop over multiple diffusion MRI measurements
            for dd in range(0,Ndmri):

                # Get MRI signal and phase information
                if(vb==True):
                    t1 = time.time()
                walk2sig = Traj2Signal(flowwalk,Gradx[dd,:],Grady[dd,:],Gradz[dd,:],deltat)   # Function providing MRI signals from spin walks and magnetic field gradients
                if(vb==True):
                    t2 = time.time()
                    print('       * Measurement {}/{}: it took {} s to synthesise signals from spins walks'.format(dd+1,Ndmri,t2-t1))
                protmag[0,dd] = walk2sig[0]      # Signal magnitude
                protphase[0,dd] = walk2sig[1]    # Signal phase

        # There are repetitions of the simulation required
        else:

            # Generate random seeds for the different simulation repetitions
            np.random.seed(seed)
            seed_list = np.random.randint(int(1e7),high=int(1e8)-1,size=Nrep)

            # Loop over simulation repetitions
            for rr in range(0,Nrep):

                # Get spin trajectories once, using a different random seed for each simulation repetition
                if(vb==True):
                    t1 = time.time()
                bufferlist = self.dMRISynMea(None,None,None,None,Nspins=Nwater,dt=deltat,Gx=Gradx[0,:],Gy=Grady[0,:],Gz=Gradz[0,:],rndseed=seed_list[rr])
                flowwalk = bufferlist[3]    # Trajectories over time of all spins
                if(vb==True):
                    t2 = time.time()
                    print('   *** Repetition {}/{}: it took {} s to generate spin walks'.format(rr+1,Nrep,t2-t1))
                # For a fixed random seed (i.e., fixed set of spin trajectories), now simulate multiple diffusion MRI measurements (it will be a lot faster!)
                for dd in range(0,Ndmri):

                    # Get MRI signal and phase information
                    if(vb==True):
                        t1 = time.time()
                    walk2sig = Traj2Signal(flowwalk,Gradx[dd,:],Grady[dd,:],Gradz[dd,:],deltat)   # Function providing MRI signals from spin walks and magnetic field gradients
                    if(vb==True):
                        t2 = time.time()
                        print('       * Measurement {}/{}: it took {} s to synthesise signals from spins walks'.format(dd+1,Ndmri,t2-t1))
                    protmag[rr,dd] = walk2sig[0]      # Signal magnitude
                    protphase[rr,dd] = walk2sig[1]    # Signal phase

        ### Return signal magnitude and signal phase for the complete protocol
        protmag = np.squeeze(protmag)
        protphase = np.squeeze(protphase)
        return protmag, protphase



    def GetTrajUniformSeed(self,Ntpoints,Nspins=2000,dt=2e-5,seednumber=20181019,boundary='periodic'):
        '''
        Calculates spin trajectories for a vascular network and uniformly-seeded spins

        rmat = obj.GetTrajUniformSeed(<parameters>)

        MANDATORY INPUT PARAMETERS
        * Ntpoints: number of time points required (number of discrete time steps at which
                    spin trajectories will be evaluated). The total duration of the simulation
                    will be (Ntpoints-1)*dt, with dt the temporal resolution expressed in s.
                    Note that the default dt is 20 us (see below).

        OPTIONAL INPUT PARAMETERS
        * Nspins:   number of spins to simulate (default: 2000).
        * dt:       temporal resolution of the simulation (scalar; units: s; default: 2e-5,
                    i.e., 20 us)
        *seednumber: Random seed used for spin position initialisation in the network. Spin positions
                     are initialised uniformly within each segment (pipe) of the network. The number
                     of spins placed within each segment its proportional to the volume of the
                     segment, with respecto to the total volume of the network. Default seed: 20181019
        * boundary: string specifying the boundary condition for particles that reach the output node.
                    It can be set to "periodic" or "feedback". Default: "periodic".
                    - "periodic" -> a particle that reaches the output node starts again
                      a new journey through a "copy" of the network that is "attached" to the output node.
                      This effectively makes the network periodic.
                    - "feedback" -> a spin that reaches the output node is fed back to the input node,
                      and starts a new journey through the vascular network. Note that
                      in this case the particle positions and velocities over time will exhibit a
                      discontinuity when the "jump" from output node to input node occurs

        RETURNS
        * rmat:     matrix storing the trajectories of the spins (units: mm). It has size
                    3 x Ntsteps x Nspins , where Nspins is the number of requested spins,
                    and Ntsetps is the number of simulation time points at a temporal resolution equal to dt.
                    Note that for conservation of total mass the vascular network a boundary condition is
                    applied depending on the value of the input "boundary". This means that spins that reach
                    the output node before the simulation is finshed, start a new journey through the network
                    i) in a "copy" of the network that is physically attached to the output node
                       (if boundary is set to "periodic"), or
                    ii) after having "jumped" back to the input node
                    (if boundary is set to "feedback). 
                    A detailed description of rmat is:
                    - rmat[0,k,s] -> x-component of the trajectory of the s-th spin, at time step k-th
                    - rmat[1,k,s] -> y-component of the trajectory of the s-th spin, at time step k-th
                    - rmat[2,k,s] -> z-component of the trajectory of the s-th spin, at time step k-th
        '''

        ### Deal with inputs
        if( not np.isscalar(Ntpoints) ):
            raise RuntimeError('ERROR. The number of time points must be a scalar integer > 0')
        if( not isinstance(Ntpoints, int) ):
            raise RuntimeError('ERROR. The number of time points must be a scalar integer > 0')
        if( Ntpoints < 0 ):
            raise RuntimeError('ERROR. The number of time points must be scalar integer > 0')
        if( not np.isscalar(Nspins) ):
            raise RuntimeError('ERROR. The number of molecules must be a scalar integer > 0')
        if( not isinstance(Nspins, int) ):
            raise RuntimeError('ERROR. The number of molecules must be a scalar integer > 0')
        if( Nspins < 0 ):
            raise RuntimeError('ERROR. The number of molecules must be scalar integer > 0')
        if( not np.isscalar(dt) ):
            raise RuntimeError('ERROR. The temporal resolution of the simulation must be a scalar > 0 expressed in seconds')
        if( dt < 0.0 ):
            raise RuntimeError('ERROR. The temporal resolution of the simulation must be a scalar > 0 expressed in seconds')
        if( not ( (boundary=='periodic') or (boundary=='feedback') ) ):
            raise RuntimeError('ERROR. The boundary condition string can only be either "periodic" or "feedback"')

        # Allocate output storing spin trajectories
        rmat = np.zeros((3,Ntpoints,Nspins))

        # Set random seed
        np.random.seed(seednumber)

        # Get the total volume of segments making up the network
        Vmat = 1000.0*1000.0*1000.0*np.pi*self.lengthmat*self.radmat*self.radmat
        Vmat[np.isnan(Vmat)] = 0.0
        Vmat = np.triu(Vmat)      # Remove elements below the diagonal - this matrix is symmetric
        Vtot = np.nansum(Vmat)    # Total volume of all segments taken together, in um3
        
        # Binary mask indicating the segments to seed
        mymask = np.abs(self.connmat)    # Get connectivity matrix in absolute value
        mymask[np.isnan(mymask)] = 0     # Remove NaNs
        mymask = np.triu(mymask)         # Remove elements below the diagonal - this matrix is anti-symmetric
    
        # Get total number of spins to seed in each segment
        SeedMap = np.ceil( float(Nspins)*Vmat/Vtot ) + 1.0
        SeedMap = SeedMap*mymask
        Nspins_actual = np.nansum(SeedMap)
        # Due to discretisation, we may have allocated too few/many spins. Make sure we allocate the exact number of spins requested
        if(Nspins_actual!=float(Nspins)):

            # Number of spins that have to be added/removed to ensure that the number of spins synthesised matches the requested number of spins
            Nspins_excess = Nspins_actual - Nspins

            # Find segment with the highest number of spins       
            maxidx = np.unravel_index(SeedMap.argmax(), SeedMap.shape)
            Nspins_maxseg_old = SeedMap[maxidx[0],maxidx[1]]
            # Add/remove spins from the branch with the highest number of spins (we are talking about a few tenths of spins over thousands, nothing to worry about)
            Nspins_maxseg_new = Nspins_maxseg_old - Nspins_excess   # If spins have to be removes, Nspins_excess is > 0; if they have to be added, is < 0
            if(Nspins_maxseg_new<=0):
                raise RuntimeError('ERROR. You have requested too few spins for uniform seeding across vascular segments. Please increase Nspins (now {}) slightly.'.format(Nspins))
            SeedMap[maxidx[0],maxidx[1]] = Nspins_maxseg_new


        # Now go on to synthesise spin walks -- loop through all segments (all possible pairs of connected nodes)
        spin_counter = 0   # Spin counter: number of spins that have been simulated
        for ii in range(0,self.totnodes):
            for jj in range(0,self.totnodes):

                # Segment detected: two nodes are connected
                if(mymask[ii,jj]):

                    # Length of current segment in mm
                    Lij = np.abs(self.lengthmat[ii,jj])

                    # No. of spins to be seeded for the current segment, which is proportionaly to its volume. If no spins were seeded in a connected branch, throw an error
                    Nij = int( SeedMap[ii,jj] )
                    if(Nij==0):
                        raise RuntimeError('ERROR. You have requested too few spins for uniform seeding across vascular segments. Please increase Nspins (now {}) slightly.'.format(Nspins))
                                    
                    # Get inlet/outlet for current segment studying the velocity vij between node ii and jj
                    vij = self.velmat[ii,jj]      # vij > 0 implies flow from i to j; vij < 0 implies flow from j to i
                    if(vij>0.0):
                        pijin = np.copy(self.nodepos[:,ii])
                        pijout = np.copy(self.nodepos[:,jj])
                        ij_idxout = jj
                    elif(vij<0.0):
                        pijin = np.copy(self.nodepos[:,jj])
                        pijout = np.copy(self.nodepos[:,ii])
                        vij = np.abs(vij)
                        ij_idxout = ii
                    else:
                        raise RuntimeError('ERROR. The velocity between node i = {} and j = {} must be > 0.0 or < 0.0, but cannot be exactly 0.0'.format(ii,jj))

                    # Get direction of 3D space connecting input and output node (unit norm)
                    dirij = (pijout - pijin)/np.linalg.norm(pijout - pijin)

                    # Loop over spins for the current segment
                    for ww in range(0,Nij):

                        # Number of times the spin has reached the outlet of the whole network - initialise to 0
                        times_reached_outlet = 0

                        # Seed spin position along the segment axis, drawing its from the uniform distribution
                        alphaval = np.random.uniform(low=0.0, high=1.0)
                        pwinit = pijin + alphaval*Lij*dirij

                        # Get trajectory from seed position to the output of the first segment where the spins travels in
                        pw = SegmentIn2Out(pwinit,pijout,vij,dt)     # pw is the trajectory of the w-th spin out of W
                        mycurnode = int(ij_idxout)                   # index of the output we just reached (output of segment ij)

                        # Continue travelling through segments node-by-node unless the trajectory is already longer than the requested duration (pw.shape[1]<=Ntpoints)
                        while(pw.shape[1]<=Ntpoints):
                            
                            # Check if we have reached the output node - if so, let's start travelling from input node again, using an appropriate boundary condition
                            if(mycurnode==self.nodeout):
                                times_reached_outlet = times_reached_outlet + 1     # Increase the counter, counting how many times we reached the output node
                                mycurnode = self.nodein                             # The current node, where we start the travel through a new segment, is again the input node

                            # Get the list of nodes where the spin could potentially travel to
                            candidate_nodes = self.connmat[mycurnode,:]         # Get a line of the connectivity matrix
                            candidate_nodes[np.isnan(candidate_nodes)] = 0      # Get rid of elements along the diagonal of the connectivity matrix
                            candidate_nodes[candidate_nodes<0] = 0              # Get rid of trajectories arrive TO the current node - we care only of those departing FROM it
                            candidate_nodes = np.arange(candidate_nodes.size)*candidate_nodes  # List of available destination nodes 
                            candidate_flows = self.flowmat[mycurnode,:]          # Now get the volumetric flow rates
                            candidate_flows[np.isnan(candidate_flows)] = 0       # Again, get rid fo elements along the diagonal
                            candidate_flows[candidate_flows<0] = 0               # Again, get rid of flows that arrive TO the current node
                            candidate_flows = candidate_flows/np.sum(candidate_flows)    # Get the probabilities of travelling to a specific new node from the current one

                            # Find which node we randomly travel to, with the probability being proportional to the amount of flow exiting towards that node
                            mydicecast = np.random.multinomial(1,candidate_flows)              # Cast a weighted dice with as many faces as potential destination nodes  (multinomial distribution)
                            mynextnode = np.arange(0,candidate_nodes.size)[mydicecast==1][0]   # This is my next node - hurray
                            posoutnew = np.copy(self.nodepos[:,mynextnode])                    # Spatial coordinates of my next node - looking forward to getting there soon
                            vnew = self.velmat[mycurnode,mynextnode]    # velocity of spins flowing from mycurnode to mynextnode - by construction velmat[mycurnode,mynextnode] is > 0, we got rid of negative connection

                            # Get new input/output for the new segment where the spin will travel through
                            posinnew = self.nodepos[:,mycurnode]
                            posoutnew = self.nodepos[:,mynextnode]

                            # If the spin had already reached the outlet, use an appropriate boundary condition
                            if(boundary=='periodic'):    # periodic --> we create virtual copies of the networks attaching their input node to the output node of the previous network
                                # Get vector connecting inlet-to-outlet and add it to the positions ...
                                # .. add it as many times as the spin has reached the outlet already
                                pinlet = self.nodepos[:,self.nodein]
                                poutlet = self.nodepos[:,self.nodeout]
                                inlet_to_outlet = poutlet - pinlet
                                for nbnd in range(0,times_reached_outlet):
                                        posinnew = posinnew + inlet_to_outlet
                                        posoutnew = posoutnew + inlet_to_outlet 
                            elif(boundary=='feedback'):   # feedback --> nothing to do. We jump back to the input node. Be aware that this will create spikes in the velocity profile
                                pass
                            else:
                                raise RuntimeError('ERROR. The requested boundary condition {} is unknown.'.format(boundary))

                            # Get the new part of the trajectory
                            pnew = SegmentIn2Out(posinnew,posoutnew,vnew,dt)   # Trajectory of the spin in the new segment
                            pnew = pnew[:,1:pnew.shape[1]]     # Get rid of the initial position, as it would be a repetition of the last position of the previous segment

                            # Update current node index (variable storing at which node we are at)
                            mycurnode = mynextnode

                            # Concatenate with the previous trajectory and continue generating new segments of trajectory for the current spin
                            pw = np.concatenate((pw,pnew),axis=1)

                        # We have finished synthesing the trajectory for the current spin ...
                        # ... now remove any part of the trajectory synthesis that has exceeded the requested simulation duration
                        rmat[:,:,spin_counter] = pw[:,0:Ntpoints]

                        # Increement the spin counter
                        spin_counter = spin_counter + 1   # We have synthesised a trajectory for one more spin

        ## We have synthesised trajectories for all requested spins - now return the synthesised spin trajectories, job done
        return rmat                
                    




    def GetPathTraj(self,dt=2e-5):
        '''
        Calculates position, velocity and flow experienced by a spin travelling along all possible paths
        that connect the input node of the vascular network to the output node

        rs, vs, qs = obj.getPathTraj()
        rs, vs, qs = obj.getPathTraj(dt)

        OPTIONAL INPUT PARAMETER
        * dt:   temporal resolution used to evaluate trajectories over time
                (scalar; units: s; default: 2e-5, i.e., 20 us)

        RETURNS
        * rs:   list storing the position over time (units: mm) of a spin along all possible P
                paths that connect the input node to the output node. The p-th element rs[p],
                for p = 0, ..., P-1 (where P = len(obj.iopaths)), is a 2D numpy array that stores
                the position of a spin over time along the p-th paths connecting the input node
                to the output node. It has size 3 x Ntpoints, where Ntpoins is the number of time
                points of duration dt required for all the spins injected to the network
                at time t = 0 to reach the output node. When a molecule reaches the output node,
                NaNs are stored in all following time points (the spin has left the network).
                - rs[p][0,k] stores the x-position of a molecule following the p-th path
                  at the k-th time step
                - rs[p][1,k] stores the y-position of a molecule following the p-th path
                  at the k-th time step
                - rs[p][2,k] stores the z-position of a molecule following the p-th path
                  at the k-th time step
        * vs:   list storing the velocity over time (units: mm/s) of a spin along all possible P
                paths that connect the input node to the output node. The p-th element vs[p],
                for p = 0, ..., P-1 (where P = len(obj.iopaths)), is a 2D numpy array that stores
                the velocity of a spin along the p-th paths connecting the input node to the
                output node. It has size 3 x Ntpoints, where Ntpoins is the number of time points of
                duration dt required for all the spins injected to the network
                at time t = 0 to reach the output node. When a molecule reaches the output node,
                NaNs are stored in all following time points (the spin has left the network).
                - vs[p][0,k] stores the x-component of the velocity of a molecule following the
                  p-th path at the k-th time step
                - vs[p][1,k] stores the y-component of the velocity of a molecule following the
                  p-th path at the k-th time step
                - vs[p][2,k] stores the z-component of the velocity of a molecule following the
                  p-th path at the k-th time step
        * qs:   list storing the volumetric flow rate over time (units: mm3/s) experienced by a
                spin while it travels along all possible P paths that connect the input node to
                the output node. The p-th element qs[p], for p = 0, ..., P-1
                (where P = len(obj.iopaths)), is a 1D numpy array that stores
                the flow experienced by a spin along the p-th paths connecting the input node to the
                output node. It has a length equal to Ntpoints elements, where Ntpoins is the
                number of time points of duration dt required for all spins injected to the network
                at time t = 0 to reach the output node. When a molecule reaches the output node,
                NaNs are stored in all following time points (the spin has left the network).
                - qs[p][k] stores the instantaneous flow rate experienced by a molecule following the
                  p-th path at the k-th time step
        '''

        # Check whether inputs make sense
        if( not np.isscalar(dt) ):
            raise RuntimeError('The temporal resolution of the wave front must be a scalar > 0')
        if( dt < 0.0 ):
            raise RuntimeError('The temporal resolution of the wave front must be a scalar > 0')

        # Get velocity matrix in mm/s, flow matrix in mm3/s, input flow in mm3/s and list of all paths from input node to output node
        Vmat = self.velmat        # Velocity matrix in mm/s
        Qmat = self.flowmat       # Flow matrix in mm3/s
        qpaths = self.iopaths     # List of all possible paths from input node to output node
        Npaths = len(qpaths)      # Number of paths

        # Get the trajectories of water moelcules along each path
        Plist = []     # Plist[ii]: position(t) of a molecule completing the ii-th path (self.iopaths[ii]), for ii = 0, ..., Npath - 1 (Npath = len(self.iopaths))
        Vlist = []     # Vlist[ii]: velocity(t) of a molecule completing the ii-th path (self.iopaths[ii]), for ii = 0, ..., Npath - 1
        Qlist = []     # Qlist[ii]: flow q(t) experienced by a molecule completing the ii-th path (self.iopaths[ii]), for ii = 0, ..., Npath - 1

        for nn in range(0,Npaths):

            mypath = qpaths[nn]              # List of nodes the current path passes through
            nsegs = len(mypath) - 1          # Number of segments particles in this path go through

            # Get trajectories of all molecules that complete the journey along all segments within the current path
            Ppath = np.zeros((3,0))   # Ppath size 3 x Nsteps_for_current_path
            Vpath = np.zeros((3,0))   # Vpath size 3 x Nsteps_for_current_path
            Qpath = np.array([])      # Qpath size Nsteps_for_current_path

            # Loop over segments making up the current path
            for bb in range(0,nsegs):

                # Get position of one molecule in the current segment along time
                nodeout = mypath[bb+1]                 # Output node index
                nodein = mypath[bb]                    # Input node index
                posout = self.nodepos[:,nodeout]       # Position of output node [mm]
                posin = self.nodepos[:,nodein]         # Position of input node [mm]
                vinout = Vmat[nodein,nodeout]          # Magnitude of velocity [mm/s] - it can be negative: we use the convention v[i][j] > 0 if i-->j, v[i][j] < 0 if j-->i
                qinout = Qmat[nodein,nodeout]          # Volumetric flow rate [mm3/s] - it can be negative: we use the convention q[i][j] > 0 if i-->j, q[i][j] < 0 if j-->i
                if(vinout<0.0):
                    vinout = (-1.0)*vinout
                if(qinout<0.0):
                    qinout = (-1.0)*qinout
                trajpipe = SegmentIn2Out(posin,posout,vinout,dt)     # 2D array of size 3 x Nsteps_of_current_segment

                # Get velocity of one molecule in the current segment along time
                if(vinout<0):
                    raise RuntimeError('The specified velocity {} is negative - it seems the flow is actually going from output to input. Please switch input/output'.format(vinout))
                vvec = vinout*( (posout - posin)/np.linalg.norm(posout - posin) )   # Vector storing the (x,y,z) components of the velocity (dx/dt,dy/dt,dz/dt)
                vpipe = np.zeros((3,trajpipe.shape[1]))             # 2D array of size 3 x Nsteps_of_current_segment
                vpipe[0,:] = vvec[0]    # x component of velocity vx = dx/dt in mm/s
                vpipe[1,:] = vvec[1]    # y component of velocity vy = dy/dt in mm/s
                vpipe[2,:] = vvec[2]    # z component of velocity vz = dz/dt in mm/s

                # Get volumetric flow rate experienced by one molecule in the current segment along time
                qpipe = qinout*np.ones(trajpipe.shape[1])          # 1D array of size Nsteps_of_current_segment

                # Stack along time trajectory, velocity and flow rate of current segment with those from the previous segment
                Ppath = np.concatenate((Ppath,trajpipe),axis=1)
                Vpath = np.concatenate((Vpath,vpipe),axis=1)
                Qpath = np.concatenate((Qpath,qpipe))

            # Store trajectory r(t) = (x(t),y(t),z(t)) for current path (mm)
            Plist.append(Ppath)  # Plist[ii]: trajectory of a molecule completing the ii-th path (self.iopaths[ii]), for ii = 0, ..., Npath - 1 (Npath = len(self.iopaths))

            # Store velocity v(t) = (dx(t)/dt,dy(t)/dt,dz(t)/dt) = (vx(t),vy(t),vz(t)) for current path (mm/s)
            Vlist.append(Vpath)  # Vlist[ii]: velocity(t) of a molecule completing the ii-th path (self.iopaths[ii]), for ii = 0, ..., Npath - 1

            # Store volumetric flow rate q(t) = dotproduct( v(t), A(r(t)) ) for current path (mm3/s), i.e., volumetric flow rate experienced at any time by the molecule
            Qlist.append(Qpath)  # Qlist[ii]: flow q(t) experienced by a molecule completing the ii-th path (self.iopaths[ii]), for ii = 0, ..., Npath - 1

        # Return, in this order: Plist, Vlist, Qlist
        return Plist, Vlist, Qlist


### Functions outside pipenet() class

def SegmentIn2Out(posin,posout,v,dt):
    '''
    Calculates the trajectory of a spin flowing in a pipe connecting two nodes

    p = SegmentIn2Out(posin,posout,v,dt)

    MANDATORY INPUT PARAMETERS
    * posin:  1D numpy array or a 1D list storing the coordinates of the pipe
              entrance (units: mm).
              posin[0]: x coordinate; posin[1]: y coordinate; posin[2]: z coordinate.
    * posout: 1D numpy array or a 1D list storing the coordinates of the pipe
              exit (units: mm).
              posout[0]: x coordinate; posout[1]: y coordinate; posout[2]: z coordinate.
    * v:      velocity of the spin within the pipe (scalar; units: mm/s)
    * dt:     temporal resolution used to evaluate the spin trajectory over time
              (scalar; units: s)

    RETURNS
    * p:      2D numpy array storing the trajectory of the psin over time. p has size
              3 x Nsteps, where Nsteps is the number of time steps required for the spin
              to reach the exit from the entrance.
              - p[0,k]: x-position of the spin at the k-th time step;
              - p[1,k]: y-position of the spin at the k-th time step;
              - p[2,k]: z-position of the spin at the k-th time step.
    '''

    # Print error if the velocity or the sampling period are not scalar
    if( not np.isscalar(v) ):
        raise RuntimeError('ERROR. The velocity must be a scalar expressed in mm/s')
    if( not np.isscalar(dt) ):
        raise RuntimeError('ERROR. The sampling period must be a scalar expressed in s')

    # Print error if velocity is negative or zero, or if the sampling period makes no sense
    if(v<0):
        raise RuntimeError('The specified velocity {} is negative - it seems the flow is actually going from output to input. Please switch input/output'.format(v))
    if(v==0.0):
        raise RuntimeError('Your velocity is exactly equal to 0.0 and no flow is happening')
    if(dt<=0.0):
        raise RuntimeError('ERROR. The sampling period must be a positive (non-zero) number expressed in s, but you passed {}'.format(dt))

    # Get information regarding the pipe where molecule positions have to be calculated
    pin = np.array(posin)     # Coordinates of pipe entrance
    pout = np.array(posout)   # Coordinates of pipe exit
    if(pin.ndim!=1):
        raise RuntimeError('ERROR. The position of the pipe entrance must be a 1D numpy array or a 1D list made of 3 elements (x,y,z')
    if(pout.ndim!=1):
        raise RuntimeError('ERROR. The position of the pipe exit must be a 1D numpy array or a 1D list made of 3 elements (x,y,z')
    if(pin.size!=3):
        raise RuntimeError('ERROR. The position of the pipe entrance must be a 1D numpy array or a 1D list made of 3 elements (x,y,z')
    if(pout.size!=3):
        raise RuntimeError('ERROR. The position of the pipe exit must be a 1D numpy array or a 1D list made of 3 elements (x,y,z')

    # Get kinematics
    l = np.sqrt( (pout[0]-pin[0])**2 + (pout[1]-pin[1])**2 + (pout[2]-pin[2])**2 )    # Pipe length in mm
    T = l/v      # Time required to travel from input to output node, in seconds (velocity v is expressed in mm/s)
    varray = v*( (pout - pin)/np.linalg.norm(pout - pin) )   # Velocity vector storing the (x,y,z) components of the velocity (vx, vy, vz) = (dx/dt,dy/dt,dz/dt)

    # Allocate useful variables
    Nsteps = int(np.floor(T/dt))   # Number of steps
    traj = np.zeros((3,Nsteps+1))  # x,y,z positions (add an extra step to store the final position with the coordinates of the output node)

    # Initial position
    traj[0,0] = pin[0]
    traj[1,0] = pin[1]
    traj[2,0] = pin[2]

    # Final position
    traj[0,Nsteps] = pout[0]
    traj[1,Nsteps] = pout[1]
    traj[2,Nsteps] = pout[2]

    # Temporal integration to get all intermediate positions
    for nn in range(1,Nsteps):
        traj[0,nn] = traj[0,nn-1] + varray[0]*dt
        traj[1,nn] = traj[1,nn-1] + varray[1]*dt
        traj[2,nn] = traj[2,nn-1] + varray[2]*dt

    # Return array of size 3 x Nsteps, storing (x,y,z) coordinate along dimension 0 (rows) and time steps along dimension 1 (columns)
    return traj




def ComputeResHagPois(piperad,pipelen,mu):
    '''
    Computes flow resistances with the Hagen-Poiseuille law.

    z = ComputeResHagPois(piperad,pipelen,mu)

    MANDATORY INPUT PARAMETERS
    * piperad: pipe radius in (units: mm). It can be a scalar or a 2D numpy array
               of size Nnodes x Nnodes. In this latter case, piperad[i,j] stores
               the radius of the pipe connecting node i with node j. It must
               hold that piperad[i,j] = piperad[j,i]
    * pipelen: pipe length in (units: mm). It can be a scalar or a 2D numpy array
               of size Nnodes x Nnodes. In this latter case, pipelen[i,j] stores
               the length of the pipe connecting node i with node j. It must
               hold that pipelen[i,j] = pipelen[j,i]
    * mu:      dynamic viscosity (units: mPa x s)

    RETURNS
    * z:      resistance (expressed in mPa x s / mm3). It is a scalar if piperad
              and pipelen are scalar; otherwise, it is a 2D numpy array where
              z[i,j] = z[j,i] stores the resistance of the pipe connecting
              node i with node j, for i = 0, ..., Nnodes - 1 and
              j = 0, ..., Nnodes - 1. 
              
    NOTE
    If the ideal Hagen-Poiseuille law (used, for example, for Newtonian fluids
    like water) is used for blood, then it should be remembered that the 
    actual viscosity of blood in capillaries is about 4x the viscosity of 
    pure plasma. In that case, we advise to set the input parameter mu to 
    4*pure_plasma_viscosity.
    '''

    # Check viscosity
    if(not np.isscalar(mu)):
        raise RuntimeError('The viscosity must be a scalar > 0')
    if(mu<=0):
        raise RuntimeError('The viscosity must be a scalar > 0')

    # Compute resistance for a scalar
    if(np.isscalar(piperad) or np.isscalar(pipelen)):

        if( (np.isscalar(piperad) and np.isscalar(pipelen)) != (np.isscalar(piperad) or np.isscalar(pipelen))  ):
            raise RuntimeError('If you are trying to compute a scalar resistance, then both radius and length must be scalars.')

        # Get geometric characteristics of the current pipe
        pipedim = 2.0*piperad       # pipe diameter in mm

        # Get ideal resistance for a Newtonian fluid
        zmat = (128.0*mu*pipelen)/(pipedim*pipedim*pipedim*pipedim*np.pi)    # Resistance expressed in mPa x s/mm3

    # Compute resistance for a network
    else:

        if(piperad.shape[0]!=piperad.shape[1]):
            raise RuntimeError('ERROR. The radius matrix R must be squared.')
        if(pipelen.shape[0]!=pipelen.shape[1]):
            raise RuntimeError('ERROR. The pipe length matrix L must be squared.')
        if(pipelen.shape[0]!=piperad.shape[0]):
            raise RuntimeError('ERROR. The pipe length matrix L and the radius matrix L must have the same size.')

        # Total number of nodes
        Nnodes = piperad.shape[0]       # Total number of nodes

        # Allocate matrix of flow resistances
        zmat = np.zeros((Nnodes,Nnodes)) + np.nan

        # Calculate resistances
        for ii in range(0,Nnodes):
            for jj in range(ii+1,Nnodes):

                if(piperad[ii,jj]!=piperad[jj,ii]):
                    raise RuntimeError('ERROR. The radius matrix R = [rij] must be s.t. rij = rji')

                if(pipelen[ii,jj]!=pipelen[jj,ii]):
                    raise RuntimeError('ERROR. The pipe length matrix L = [lij] must be s.t. lij = lji')

                # Get geometric characteristics of the current pipe
                radij = piperad[ii,jj]      # pipe radius in mm
                dimij = 2.0*radij           # pipe diameter in mm
                lenij = pipelen[ii,jj]      # pipe length in mm

                # Get ideal resistance for a Newtonian fluid
                zij = (128.0*mu*lenij)/(dimij*dimij*dimij*dimij*np.pi)    # Resistance expressed in mPa x s/mm3

                # Store resistance to flow of the current pipe
                zmat[ii,jj] = zij
                zmat[jj,ii] = zij

    # Return flow resistance matrix/scalar in mPa x s/mm3
    return zmat




def ComputeResModHagPois(piperad,pipelen,mu):
    '''
    Computes flow resistances with the modified Hagen-Poiseuille law as in Blinder et al,
    Nat Neurosci 2013, 16(7): 889-897, doi: 10.1038/nn.3426

    z = ComputeResModHagPois(piperad,pipelen,mu)

    MANDATORY INPUT PARAMETERS
    * piperad: pipe radius in (units: mm). It can be a scalar or a 2D numpy array
               of size Nnodes x Nnodes. In this latter case, piperad[i,j] stores
               the radius of the pipe connecting node i with node j. It must
               hold that piperad[i,j] = piperad[j,i]
    * pipelen: pipe length in (units: mm). It can be a scalar or a 2D numpy array
               of size Nnodes x Nnodes. In this latter case, pipelen[i,j] stores
               the length of the pipe connecting node i with node j. It must
               hold that pipelen[i,j] = pipelen[j,i]
    * mu:      dynamic viscosity of pure blood plasma (units: mPa x s)

    RETURNS
    * z:      resistance (expressed in mPa x s / mm3). It is a scalar if piperad
              and pipelen are scalar; otherwise, it is a 2D numpy array where
              z[i,j] = z[j,i] stores the resistance of the pipe connecting
              node i with node j, for i = 0, ..., Nnodes - 1 and
              j = 0, ..., Nnodes - 1. Note that z depends on the effective
              viscosity of blood in a microvessel, which is approximately 4x the
              viscosity of pure blood plasma
    '''

    # Check viscosity
    if(not np.isscalar(mu)):
        raise RuntimeError('The viscosity must be a scalar > 0')
    if(mu<=0):
        raise RuntimeError('The viscosity must be a scalar > 0')

    # Compute resistance for a scalar
    if(np.isscalar(piperad) or np.isscalar(pipelen)):

        if( (np.isscalar(piperad) and np.isscalar(pipelen)) != (np.isscalar(piperad) or np.isscalar(pipelen))  ):
            raise RuntimeError('If you are trying to compute a scalar resistance, then both radius and length must be scalars.')

        # Get geometric characteristics of the current pipe
        pipedim = 2.0*piperad       # pipe diameter in mm
        rad_um = 1000.0*piperad     # pipe radius in um

        # Get resistance with the modified Hagen-Poiseuille equation as in Blinder P et al, Nat Neurosci 2013, doi: 10.1038/nn.3426
        zmat = 4.0*(1.0 - 0.863*np.exp(-rad_um/14.3) + 27.5*np.exp(-rad_um/0.351))*( (128.0*mu*pipelen)/(pipedim*pipedim*pipedim*pipedim*np.pi) )    # Resistance expressed in mPa x s/mm3

    # Compute resistance for a network
    else:

        if(piperad.shape[0]!=piperad.shape[1]):
            raise RuntimeError('ERROR. The radius matrix R must be squared.')
        if(pipelen.shape[0]!=pipelen.shape[1]):
            raise RuntimeError('ERROR. The pipe length matrix L must be squared.')
        if(pipelen.shape[0]!=piperad.shape[0]):
            raise RuntimeError('ERROR. The pipe length matrix L and the radius matrix L must have the same size.')

        # Total number of nodes
        Nnodes = piperad.shape[0]       # Total number of nodes

        # Allocate matrix of flow resistances
        zmat = np.zeros((Nnodes,Nnodes)) + np.nan

        # Calculate resistances
        for ii in range(0,Nnodes):
            for jj in range(ii+1,Nnodes):

                # Get geometric characteristics of the current pipe
                radij = piperad[ii,jj]      # pipe radius in mm
                dimij = 2.0*radij           # pipe diameter in mm
                radij_um = 1000.0*radij     # pipe radius in um
                lenij = pipelen[ii,jj]      # pipe length in mm

                # Get resistance with the modified Hagen-Poiseuille equation as in Blinder P et al, Nat Neurosci 2013, doi: 10.1038/nn.3426
                zij = 4.0*(1.0 - 0.863*np.exp(-radij_um/14.3) + 27.5*np.exp(-radij_um/0.351))*( (128.0*mu*lenij)/(dimij*dimij*dimij*dimij*np.pi) )   # Resistance expressed in mPa x s/mm3

                # Store resistance to flow of the current pipe
                zmat[ii,jj] = zij
                zmat[jj,ii] = zij

    # Return flow resistance matrix/scalar in mPa x s/mm3
    return zmat




def ComputeFlow(radmat,resmat,nodein,nodeout,qin,solver='numerical'):
    '''
    Computes the volumetric flow rate between each pair of 2 nodes in a vascular network

    Q = ComputeFlow(radmat,resmat,nodein,nodeout,qin,solver)

    MANDATORY INPUT PARAMETERS
    * radmat:  2D numpy matrix storing the radius (units: mm) of the pipes connecting
               the nodes of the vascular network. radmat[i,j] stores the radius of
               the pipe connecting node i with node j (it must hold that
               radmat[i,j] = radmat[j,i]; radmat[i,j] = 0 if the two nodes
               are not connected). Given Nnodes, i = 0, ..., Nnodes - 1 and
               j = 0, ..., Nnodes - 1
    * resmat:  2D numpy matrix storing the resistance (units: mPa x s/mm3) of the
               pipes connecting the nodes of the vascular network. resmat[i,j]
               stores the resistance of the pipe connecting node i with node j
               (it must hold that resmat[i,j] = resmat[j,i]). Given Nnodes,
               i = 0, ..., Nnodes - 1 and j = 0, ..., Nnodes - 1
    * nodein:  index of the input node (must be one of 0, ..., Nnodes - 1)
    * nodeout: index of the output node (must be one of 0, ..., Nnodes - 1;
               nodeout cannot be the same as nodein)

    OPTIONAL INPUT PARAMETERS
    * solver: computation library used to solver the vacular network. Under 40 segments,
              set to 'symbolic' for more accurate calulcations. Default: 'numerical'.

    RETURNS
    * Q        volumetric flow rate matrix Q = [qij] (units: mm3/s) of size Nnodes x Nnodes;
               element Q[i][j] = qij stores the volumetric flow rate between node i
               and node j. NaNs are stored along the diagonal. Moreover:
               qij > 0 --> flow from i to j;
               qij = < 0 --> flow from j to i;
               qij = -qji
    '''

    # Check whether inputs make sense
    if( not np.isscalar(qin) ):
        raise RuntimeError('ERROR. The input flow must be a scalar > 0')
    if( qin <= 0.0 ):
        raise RuntimeError('ERROR. The input flow must be a scalar > 0')

    if( (radmat.ndim!=2) or (resmat.ndim!=2) ):
        raise RuntimeError('ERROR. The radius and resistance matrices must be 2D, and must be square (minimum size 2x2)')
    if( radmat.shape[0]!=radmat.shape[1]):
        raise RuntimeError('ERROR. The radius and resistance matrices must be 2D, and must be square (minimum size 2x2)')
    if( resmat.shape[0]!=resmat.shape[1]):
        raise RuntimeError('ERROR. The radius and resistance matrices must be 2D, and must be square (minimum size 2x2)')

    if ( not np.isscalar(nodein)):
        raise RuntimeError('ERROR. The index of the input node must be a scalar integer')
    if ( not np.isscalar(nodeout)):
        raise RuntimeError('ERROR. The index of the output node must be a scalar integer')
    if( (nodein<0) or (nodeout<0) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')
    if( not isinstance(nodein,int) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')
    if( not isinstance(nodeout,int) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')
    if( (nodein>=radmat.shape[0]) or (nodeout>=radmat.shape[0]) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')
    if(nodein==nodeout):
        raise RuntimeError('ERROR. The output node cannot be the same as the input node.')
    if (solver=='symbolic'):
        # Initialise a graph of resitors with the lcapy python package
        cnet = lcapy.Circuit()

        # Add resistors according to the given vascular network structure
        Nnodes = radmat.shape[0]
        for ii in range(0,Nnodes):
            for jj in range(ii+1,Nnodes):
                radij = radmat[ii,jj]
                if( (~np.isnan(radij)) and (radij!=0) ):
                    rij = resmat[ii,jj]
                    cnet.add('R{}{} {} {} {}'.format(ii,jj,ii,jj,rij))

        # Add a current generator            
        cnet.add('Iq {} {} {}'.format(nodein,nodeout,qin))

        # Get electric currents through each link connecting a pair of nodes and build a volumetric flow matrix
        qmat = np.zeros((Nnodes,Nnodes))
        for ii in range(0,Nnodes):
            for jj in range(ii+1,Nnodes):
                radij = radmat[ii,jj]
                if( (~np.isnan(radij)) and (radij!=0) ):
                    currij = cnet.get_I('R{}{}'.format(ii,jj))   # Get electric current between node i and j
                    mystrij = '{}'.format(currij.evalf(16))      # Convert to numeric (precision: 16 decimal points)
                    try:
                        myvalij = float( mystrij[mystrij.find(' ')+1:len(mystrij)-2] )
                    except:
                        myvalij = 0.0
                    qmat[ii,jj] = myvalij
                    qmat[jj,ii] = -myvalij
        
        # Put NaNs along the diagonal
        for ii in range(0,Nnodes):
            qmat[ii,ii] = np.nan
            
        # Return solution
        return qmat
    else:
        # # Initialise a circuit with the PySpice python package
        cnet = Circuit('')

        # # Add a current generator
        cnet.I('q',f'{nodein}', f'{nodeout}', qin)

        # # Add resistors according to the given vascular network structure
        Nnodes = radmat.shape[0]
        for ii in range(0,Nnodes):
            for jj in range(ii+1,Nnodes):
                radij = radmat[ii,jj]
                if( (~np.isnan(radij)) and (radij!=0) ):
                    rij = resmat[ii,jj]
                    cnet.R(f'{ii}{jj}',f'{ii}', f'{jj}', rij).minus.add_current_probe(cnet)        

        # # Create analysis of circuit 
        simulator = cnet.simulator(temperature=25, nominal_temperature=25)
        analysis = simulator.operating_point()

        # # Get electric currents through each link connecting a pair of nodes and build a volumetric flow matrix
        qmat = np.zeros((Nnodes,Nnodes))
        for key in analysis.branches.keys():
            # kkk = str(key)
            if len(key) == 10:
                ii, jj = int(key[2]), int(key[3])
            elif len(key) == 11:
                ii, jj = int(key[2]), int(key[3:5])
            elif len(key) == 12:
                ii, jj = int(key[2:4]), int(key[4:6])
            else:
                continue
            
            myvalij = float(analysis.branches[key])
            qmat[ii, jj] = myvalij
            qmat[jj, ii] = -myvalij

        # Put NaNs along the diagonal
        for ii in range(0,Nnodes):
            qmat[ii,ii] = np.nan
        return qmat




def ComputeAllPaths(C, idxin, idxout):
    '''
    Computes all paths connecting the input node to the output node given a connecivity matrix.
    Returns a directed graph given the connectivity matrix.

    pathlist = ComputeAllPaths(C, idxin, idxout)

    MANDATORY INPUT PARAMETERS
    * C:   connectivity matrix C = [cij] of size Nnodes x Nnodes; element C[i,j] is the connectivity cij
           between node i and node j, should be is +1/-1 if two nodes are connected, 0 otherwise.
           NaNs should be stored along the diagonal;
           cij = 1 --> flow from i to j; cij = -1 --> flow from j to i;
           cij = -cji
    * idxin:  index of the input node (must be one of 0, ..., Nnodes - 1)
    * idxout: index of the output node (must be one of 0, ..., Nnodes - 1;
               nodeout cannot be the same as nodein)

    RETURNS
    * pathlist: list of P elements containing all possible paths from the input node to the output
                node. Each element pathlist[p] contains a description for the p-th path, with
                p = 0, ..., P - 1, in the form of a list storing the sequence of nodes that
                are encountered by a molecule following path p-th. For example, assume pathlist
                stores something like [[1,2], [1,0,2]]; this would imply that there are two possible
                paths between the input node 1 and the output node 2, and these would be [1,2]
                (path number 0) and [1,0,2] (path number 1)
    * gnet:     a directed Graph object populated with vertices ('node_list') and edges ('edge_list') 
                from the connectivity matrix C.

    '''

    # Check whether inputs make sense
    if( C.ndim!=2 ):
        raise RuntimeError('ERROR. The connectivity matrix must be 2D, and must be square (minimum size 2x2)')
    if( C.shape[0]!=C.shape[1]):
        raise RuntimeError('ERROR. The connectivity matrix must be 2D, and must be square (minimum size 2x2)')

    if ( not np.isscalar(idxin)):
        raise RuntimeError('ERROR. The index of the input node must be a scalar integer')
    if ( not np.isscalar(idxout)):
        raise RuntimeError('ERROR. The index of the output node must be a scalar integer')
    if( (idxin<0) or (idxout<0) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')
    if( (idxin>=C.shape[0]) or (idxout>=C.shape[0]) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')
    if( not isinstance(idxin,int) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')
    if( not isinstance(idxout,int) ):
        raise RuntimeError('ERROR. The indices of the input/output nodes must be integers ranging in 0, ..., number_of_nodes - 1')

    # Build a directed graph given the connectivity matrix
    gnet = gt.Graph()
    node_list = []
    edge_list = []
    Nnodes = C.shape[0]

    # Add all nodes
    for ii in range(0,Nnodes):
        myedge = gnet.add_vertex()
        node_list.append(myedge)

    # Add all between-node connections (known as pipes or segments)
    for ii in range(0,Nnodes):
        for jj in range(ii+1,Nnodes):
            cij = C[ii,jj]
            cji = C[jj,ii]

            # There is a connection
            if(np.abs(cij)!=0):
                if(cij!=-cji):
                    raise RuntimeError('ERROR. The connectivity matrix must satisfy cji = -cij, with cij = {-1, 0, 1}')

                if(cij>0):
                    # Add directed link from node ii to node jj if cij > 0
                    myedge = gnet.add_edge(node_list[ii], node_list[jj])
                else:
                    # Add directed link from node jj to node ii if if cij < 0
                    myedge = gnet.add_edge(node_list[jj], node_list[ii])

                # Store new connection in the graph
                edge_list.append(myedge)

            # There is not a connection
            else:
                if(cji!=0.0):
                    raise RuntimeError('ERROR. If the connectivity matrix features cij = 0, then cji should also be 0')

    # Get all paths from input to output given the current directed graph
    allpaths = gtop.all_paths(gnet,node_list[idxin],node_list[idxout])
    allpath_list = []
    for path in allpaths:
        allpath_list.append(path)

    # Return a list of arrays, where each array stores the nodes to be followed along each path going from input to output,
    # Returns directed graph populated with vertices and edges 
    return allpath_list, gnet



def Traj2Signal(spinpos,gradx,grady,gradz,deltatime):
    '''
    Synthesises a diffusion-weighted MRI signal from a set of spin trajectories

    mag, phase, spinphase = Traj2Signal(spinpos,gradx,grady,gradz,deltatime)

    MANDATORY INPUT PARAMETERS
    * spinpos:      a 3D numpy matrix storing the trajectories of the spins (units: mm). It has size
                    3 x Ntsteps x Nmol, where Nmol is the number of spins, and Ntsetps is the number of simulation
                    time points at a temporal resolution equal to "deltatime" (see input deltatime below).
                    - spinpos[0,k,s] stores the x-component of the trajectory of the s-th spin, at time step k-th
                    - spinpos[1,k,s] stores the y-component of the trajectory of the s-th spin, at time step k-th
                    - spinpos[2,k,s] stores the z-component of the trajectory of the s-th spin, at time step k-th
    * gradx:       a 1D array storing the x component of the diffusion gradient (units: T/mm) over time.
                   It must have length Ntsteps, so that gradx[k] is the gradient at time step k-th
    * grady:       a 1D array storing the y component of the diffusion gradient (units: T/mm) over time.
                   It must have length Ntsteps, so that gradx[k] is the gradient at time step k-th
    * gradz:       a 1D array storing the z component of the diffusion gradient (units: T/mm) over time.
                   It must have length Ntsteps, so that gradx[k] is the gradient at time step k-th
    * deltatime:   temporal resolution used to sample the spin trajectories and the diffusion gradients
                  (units: s)

    RETURNS
    * mag:         a scalar storing the magnitude of the diffusion-weighted signal (mag = 1.0 when there
                   is no diffusion-weighting)
    * phase:       a scalar storing the phase of the diffusion-weighted signal (range: [-pi; pi])
    * spinphase:   1D array storing the phase accrual of all individual spins. The array has size Nmol,
                   so that phasespin[s] stores the phase accrual of the s-th spin, for s = 0, ..., Nmol - 1.
                   It ranges between ]-inf; +inf[ so that one can study the number of phase wraps caused
                   by diffusion-weighting
    '''

    # Check that everything makes sense
    if(spinpos.ndim!=3):
        raise RuntimeError('ERROR. The spin trajectories must be stored as a 3D numpy array of size 3 x number_of_time_points x number_of_spins')
    if(spinpos.shape[0]!=3):
        raise RuntimeError('ERROR. The 0-th dimension of the spin trajectories spinpos must be 3, i.e., spinpos.shape[0] must be 3')
    if( (gradx.ndim!=1) or (grady.ndim!=1) or (gradz.ndim!=1) ):
        raise RuntimeError('ERROR. Each of gx(t), gy(t), gz(t) must be a 1D array')
    if( (gradx.size!=grady.size) or (gradx.size!=gradz.size)  ):
        raise RuntimeError('ERROR. gx, gy and gz do not agree in size')
    if( gradx.size!=spinpos.shape[1] ):
        raise RuntimeError('ERROR. The duration of the gradient waveform does not match the temporal duration of the spin trajectories')
    if( not np.isscalar(deltatime) ):
        raise RuntimeError('ERROR. The sampling period deltatime must be a scalar > 0 expressed in [s]')
    if(deltatime<=0.0):
        raise RuntimeError('ERROR. The sampling period deltatime must be a scalar > 0 expressed in [s]')

    # Get phase of all molecules
    Nspins = spinpos.shape[2]
    spinphase = np.zeros(Nspins)                            # Array that will store the phase accrued by each spin
    spinsig = np.zeros(Nspins) + 1j*np.zeros(Nspins)        # Array that will store the complex-valued signal of each spin
    gammaproton = 267.522187*1e6                            # Proton gyromagnetic ratio in 1/(s x T)
    for kk in range(0,Nspins):

        # Get x(t), y(t), z(t) for current spin
        xt = spinpos[0,:,kk]
        yt = spinpos[1,:,kk]
        zt = spinpos[2,:,kk]

        # Get phase accrual for current spin -->  phase = - gamma * integral( dot_product(G(t),trajectory(t))) dt
        phaseacc = -1.0*gammaproton*deltatime*np.sum( gradx*xt + grady*yt + gradz*zt )   # Phase accrual of the kk-th spin; gradients in T/mm, trajectories in mm, time in s
        spinsig[kk] = np.exp(1j*phaseacc)          # Complex-valued signal of the kk-th spin
        spinphase[kk] = phaseacc                   # Phase accrual of the kk-th spin - can take any value in ]-inf; +inf[

    # Compute total complex-valued signal as well as its magnitude and phase
    sig = np.mean(spinsig)
    mag = np.sqrt( np.real(sig)*np.real(sig) + np.imag(sig)*np.imag(sig) )
    phase = np.angle(sig)    # Phase of the total MRI signal - this is defined in the range [-pi; pi]

    # Return total signal magnitude, total signal phase, and phase of each individual spin
    return mag, phase, spinphase



def getGradPGSE(bval, Gdur, Gsep, Gdir, dt, Nsample=None):
    '''
    Generates a waveform for an ideal pulsed-gradient spin echo (PGSE) sequence

    grx, gry, grz, tm = getGradPGSE(bval, Gdur, Gsep, Gdir, dt, Nsample=None)

    MANDATORY INPUT PARAMETERS
    * bval: b-value for a PGSE sequence (scalar; units: s/mm2);
    * Gdur: gradient duration (small delta) for a PGSE sequence (scalar; units: s)
    * Gsep: gradient separation (large Delta) for a PGSE sequence (scalar; units: s)
    * Gdir: gradient direction for a PGSE sequence; it must be a 1D numpy array of size 3, storing the
            X, y and z components of the gradient direction
    * dt: temporal resolution to be used to sample the diffusion gradient wavform (units: s)

    OPTIONAL INPUT PARAMETERS
    * Nsample: number of temporal samples to be used to generate the gradient waveform. If Nsample
               is larger than the minimum number of samples required to construct a
               waveform with gradien duration Gdur and gradient separation Gsep given the temporal
               resolution dt, then a tail of zeros will be added to the generated gradient waveform.
               If Nsample is smaller than the minimum number of samples required to build the
               waveform, then the latter will be used. Default: None.

    RETURNS
    * grx:       a 1D array storing the x component of the diffusion gradient (units: T/mm) over time
    * gry:       a 1D array storing the y component of the diffusion gradient (units: T/mm) over time
    * grz:       a 1D array storing the z component of the diffusion gradient (units: T/mm) over time
    * tm:        a 1D array storing the time array, s.t. tm[0] = 0, tm[1] = dt, tm[2] = 2*dt, etc

    '''

    # Make sure inputs make sense
    if( not np.isscalar(bval) ):
        raise RuntimeError('ERROR. The b-value must be a scalar (units: s/mm2), but you passed {} s/mm2'.format(bval))
    if( not np.isscalar(Gdur) ):
        raise RuntimeError('ERROR. The gradient duration must be a scalar (units: s), but you passed {} s'.format(Gdur))
    if( not np.isscalar(Gsep) ):
        raise RuntimeError('ERROR. The gradient separation must be a scalar (units: s), but you passed {} s'.format(Gsep))
    if(bval<0):
        raise RuntimeError('ERROR. A negative b-value does not make sense, but you passed {} s/mm2'.format(bval))
    if(Gdur<0):
        raise RuntimeError('ERROR. A negative gradient duration does not make sense, but you passed {} s'.format(bval))
    if(Gsep<0):
        raise RuntimeError('ERROR. A negative gradient separation does not make sense, but you passed {} s'.format(bval))
    if(Gsep<Gdur):
        raise RuntimeError('ERROR. The gradient separation cannot be shorter than the duration - you passed a duration of {} and a separation of {} s'.format(Gdur,Gsep))
    if( (Gdir.ndim!=1) or (Gdir.size!=3) ):
        raise RuntimeError('ERROR. The gradient direction must be a 1D numpy array of size 3, but you passed {}'.format(Gdir))
    if( not np.isscalar(dt) ):
        raise RuntimeError('ERROR. The sampling period must be a scalar expressed in s')
    if(dt<0.0):
        raise RuntimeError('ERROR. The sampling period must be a number > 0 expressed in s')


    # Gradient strength in T/mm
    gammar = 267.522187*1e6                      # Proton gyromagnetic ratio in 1/(s x T)

    # Deal with a b = 0.0 measurement
    if( (bval==0.0) or (Gdur==0.0) or (Gsep==0.0) or (np.linalg.norm(Gdir)==0.0) ):

        if(Nsample is not None):
            grx = np.zeros(Nsample)
            gry = np.zeros(Nsample)
            grz = np.zeros(Nsample)
            timearr = dt*np.linspace(0,float(Nsample)-1,Nsample)    # Time array (units: s)
        else:
            raise RuntimeError('ERROR. For a b = 0 measurement, the number of time points of the gradient waveform Nsample cannot be None')

    # Deal with a DW measurement
    else:

        # Get the required gradient strength
        Gval = np.sqrt( bval / ( gammar*gammar * Gdur*Gdur * (Gsep - Gdur/3.0) ) )	 # Gradient strength in [T/mm] (e.g., 65 mT/m --> 65*1e-6 T/mm)

        # Make sure the gradient direction has unit norm
        gn = Gdir/np.linalg.norm(Gdir)

        # Generate the waveform of the gradient magnitude for PGSE (single linear tensor encoding) - generate the minimum length required given the requested diffusion times and temporal resolution
        Nt_dur = int(np.round(Gdur/dt) + 1)    # Number of time points required to cover the first gradient lobe (gradient duration)
        Nt_sep = int(np.round(Gsep/dt) + 1)    # Number of time points required to cover the gradient up to the rise of the second lobe (gradient separation)
        if(Nt_sep<=Nt_dur):
            Nt_sep = Nt_dur + 1
        Gpos = Gval*np.ones(Nt_dur)
        Gbuffer = np.zeros(Nt_sep)
        Gbuffer[0:Nt_dur] = Gpos
        Gtime = np.concatenate((Gbuffer,-1.0*Gpos))               # Gtime: ideal discrete diffusion gradient waveform G(t) s.t. sum_t(G(t)) = 0.0
        Ntpoints_initial = Gtime.size                             # Number of time points required for the simulation

        # Zero-padd the gradient waveform at the tail if a longer time series was required
        if(Nsample is not None):
            if(Ntpoints_initial>=Nsample):
                Ntpoints_final = Ntpoints_initial     # If the generated waveform is already longer than the requested number of t-samples, just keep this longer one - it is the minimum you need for this PGSE measurement
            else:
                Nmissing = Nsample - Ntpoints_initial
                Gtime = np.concatenate((Gtime,np.zeros(Nmissing)))    # Zero-padd at the end if the requested number of t-samples is larger than the duration of the gradient waveform obtained so far
                Ntpoints_final = Gtime.size
        else:
            Ntpoints_final = Ntpoints_initial

        # Get final time array
        timearr = dt*np.linspace(0,float(Ntpoints_final)-1,Ntpoints_final)    # Time array (units: s)

        # Get time-dependent components Gx(t), Gy(t), Gz(t)
        grx = Gtime*gn[0]
        gry = Gtime*gn[1]
        grz = Gtime*gn[2]

    # Return Gx(t), Gy(t), Gz(t) and time array
    return grx, gry, grz, timearr


