The [syn.py](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/syn.py) module contains the _pipenet_ class. Find below a detailed description of the class, of its methods, and of its attributes.

```
Help on class pipenet in module syn:

class pipenet(builtins.object)
 |  pipenet(nodes, radii, qin, idxin, idxout, radiusin=None, visc=1.2, flowmodel='blinder', solver='numerical')
 |  
 |  Class to define a vascular network and use it to simulate IVIM diffusion MRI signals
 |  
 |  obj = pipenet(<parameters>) initialises a vascular network made of connected pipes
 |  
 |  MANDATORY INPUT PARAMETERS
 |  * nodes:     matrix of node positions (in mm), of size 3 x Nnodes (rows: xpos, ypos, zpos in mm);
 |               a node is defined as the input (or as the output) of a pipe
 |  * radii:     matrix of pipe radii (in mm), of size Nnodes x Nnodes (element (i,j) stores the radius
 |               of the pipe connecting node i with node j)
 |  * qin:       input volumetric flow rate in mm3/s
 |  * idxin:     index of the input node in the node matrix (nodes[:,idxin] provides the x,y,z coordinates
 |               of such a input node)
 |  * idxout:    index of the output node in the node matrix (nodes[:,idxout] provides the x,y,z coordinates
 |               of such an input node)
 |  
 |  OPTIONAL INPUT PARAMETERS
 |  * radiusin:  radius of the pipe bringing fluid to the input node in mm (default: None).
 |               If None, the this will be calculated as the mean of all the
 |               radii for all pipes that emanate from the input node
 |  * visc:      ideal dinaymic viscosity of the fluid in mPa x s (default: 1.20 mPa x s).
 |               If blood is modelled, this value should be set to:
 |              * 4 x the viscosity of pure blood plasma (when no RBCs are present), for flowmodel
 |                "ideal" (the viscosity of pure blood plasma is approx. 1.20 mPa x s, so
 |                4 x 1.20 = 4.80 mPa x s)
 |              * exactly the viscosity of pure blood plasma (when no RBCs are present;
 |                approx. 1.20 mPa x s) for model "blinder"
 |  * flowmodel: string indicating the model that should be used to calculate flow resistance
 |               in each pipe/segment of the vascular network. Values could be:
 |               - "ideal":   to use the Hagen-Poiseuille law for ideal Netwonian fluids;
 |               - "blinder": to use the modified Hagen-Poiseuille law as in Blinder et al,
 |                            Nat Neurosci 2013, 16(7): 889-897, doi: 10.1038/nn.3426. In this model,
 |                            the viscosity of the blood is a function of the vessel radius; for very
 |                            large vessels, the viscosity approaches 4 x the viscosity of pure blood
 |                            plasma when no red blood cells (RBCs) are present;
 |               Default value is "blinder"
 |  * solver:   string indicating the library package that should be used to calculate the flow matrix
 |              of the vascular network. Values could be: 
 |               - "numerical": to use the PySpice package;
 |               - "symbolic": to use the lcapy symbolic computation package;
 |               Default value is "numerical"
 |  
 |  RETURNS
 |  * obj:       instance of class pipenet, storing an initialised vascular network.
 |               It features the following attributes:
 |               - flowmodel: a string indicating the hemodynamics model (e.g., "ideal", "blinder").
 |                            It is a copy of the input parameter "flowmodel"
 |               - nodepos: a 3 x Nnodes matrix storing the positions of the nodes (in mm); a node is defined
 |                          as the input (or as the output) of a pipe. It is a copy of the
 |                          input parameter "nodes"
 |                          nodepos[0,nn] stores the x-component of the nn-th node (units: mm),
 |                          nodepos[1,nn] stores the y-component of the nn-th node (units: mm),
 |                          nodepos[2,nn] stores the z-component of the nn-th node (units: mm).
 |               - totnodes: total number of nodes (Nnodes)
 |               - nodein:  index of the input node (indices start from 0). It stores the value passed
 |                          with input parameter "idxin"
 |               - nodeout: index of the output node (indices start from 0). It stores the value passed
 |                          with input parameter "idxout"
 |               - connmat: connectivity matrix C = [cij] of size Nnodes x Nnodes; element C[i,j] is +1/-1
 |                          if two nodes are connected, 0 otherwise, with NaNs being stored along the diagonal;
 |                          cij = 1 --> flow from i to j; cij = -1 --> flow from j to i;
 |                          cij = -cji
 |               - flowmat: volumetric flow rate matrix Q = [qij] of size Nnodes x Nnodes; element Q[i,j] stores
 |                          the volumetric flow rate qij between node i and node j (units: mm3/s);
 |                          NaNs are stored along the diagonal;
 |                          qij > 0 --> flow from i to j; qij = < 0 --> flow from j to i;
 |                          qij = -qji
 |               - velmat:  flow velocity matrix V = [vij] of size Nnodes x Nnodes; element V[i,j] stores
 |                          the velocity of the flow vij between node i and node j (units: mm/s);
 |                          NaNs are stored along the diagonal;
 |                          vij > 0 --> flow from i to j; vij = < 0 --> flow from j to i;
 |                          vij = -vji
 |               - qin:     input volumetric flow rate (units: mm3/s) arriving at the input node (see nodein)
 |               - radiusin: radius of the pipe bringing fluid to the input node (units: mm). It is a copy
 |                          of input parameter "radiusin" if it was provided; if not, radiusin is set as the
 |                          mean of the radii of all branches emanating from the input node
 |               - radmat:  matrix of pipe radii R = [rij] of size Nnodes x Nnodes. rij stores the radius
 |                          of the pipe connecting node i with node j, s.t. rij = rji (units: mm).
 |                          NaNs are stored along the diagonal. It is essentially a copy of the input
 |                          parameter "radii", with NaNs along the diagonal
 |              - lengthmat: matrix of pipe lengths L = [lij] of size Nnodes x Nnodes. lij stores the length
 |                           radius of the pipe connecting node i with node j, s.t. lij = lji (units: mm).
 |                           NaNs are stored along the diagonal
 |              - resmat:   matrix of pipe resistances Z = [zij] of size Nnodes x Nnodes. zij stores the
 |                          resistance experienced by the flow between nodes i and j, s.t. zij = zji
 |                          (units: mPa x s/mm3). The way resmat is computed depends on the adopted
 |                          hemodynamics model (see flowmodel above)
 |              - iopaths:  list of P elements containing all possible paths from the input node to the output
 |                          node. Each element iopaths[p] contains a description for the p-th path, with
 |                          p = 0, ..., P - 1, in the form of a list storing the sequence of nodes that
 |                          are encountered by a molecule following path p-th. For example, assume iopaths
 |                          stores something like [[1,2], [1,0,2]]; this would imply that there are two possible
 |                          paths between the input node 1 and the output node 2, and these would be [1,2]
 |                          (path number 0) and [1,0,2] (path number 1)
 |              - gnet:     a directed Graph object from the graph_tool module populated with vertices ('node_list')
 |                          and edges ('edge_list') from the connectivity matrix C.
 |              - visco:    dynamic viscosity of the fluid (units: mPa x s). It is a copy of input parameter "visc"
 |              - solver:   a string indicating the computation package (e.g., "numerical", "symbolic").
 |                            It is a copy of the input parameter "solver"
 |  
 |  DEPENDENCIES
 |  To use pipenet you need:
 |  * numpy (developed with numpy 1.23.5; https://numpy.org)
 |  * lcapy (developed with lcapy 1.10; https://lcapy.readthedocs.io/en/latest)
 |  * PySpice (developed with PySpice 1.5; https://github.com/FabriceSalvaire/PySpice)
 |  * graph-tool (developed with graphtool 2.45, commit b1a649d8; https://graph-tool.skewed.de)
 |  * time (python standard library)
 |  pipenet was developed with python 3.10.8
 |  
 |  Methods defined here:
 |  
 |  GetPathTraj(self, dt=2e-05)
 |      Calculates position, velocity and flow experienced by a spin travelling along all possible paths
 |      that connect the input node of the vascular network to the output node
 |      
 |      rs, vs, qs = obj.getPathTraj()
 |      rs, vs, qs = obj.getPathTraj(dt)
 |      
 |      OPTIONAL INPUT PARAMETER
 |      * dt:   temporal resolution used to evaluate trajectories over time
 |              (scalar; units: s; default: 2e-5, i.e., 20 us)
 |      
 |      RETURNS
 |      * rs:   list storing the position over time (units: mm) of a spin along all possible P
 |              paths that connect the input node to the output node. The p-th element rs[p],
 |              for p = 0, ..., P-1 (where P = len(obj.iopaths)), is a 2D numpy array that stores
 |              the position of a spin over time along the p-th paths connecting the input node
 |              to the output node. It has size 3 x Ntpoints, where Ntpoins is the number of time
 |              points of duration dt required for all the spins injected to the network
 |              at time t = 0 to reach the output node. When a molecule reaches the output node,
 |              NaNs are stored in all following time points (the spin has left the network).
 |              - rs[p][0,k] stores the x-position of a molecule following the p-th path
 |                at the k-th time step
 |              - rs[p][1,k] stores the y-position of a molecule following the p-th path
 |                at the k-th time step
 |              - rs[p][2,k] stores the z-position of a molecule following the p-th path
 |                at the k-th time step
 |      * vs:   list storing the velocity over time (units: mm/s) of a spin along all possible P
 |              paths that connect the input node to the output node. The p-th element vs[p],
 |              for p = 0, ..., P-1 (where P = len(obj.iopaths)), is a 2D numpy array that stores
 |              the velocity of a spin along the p-th paths connecting the input node to the
 |              output node. It has size 3 x Ntpoints, where Ntpoins is the number of time points of
 |              duration dt required for all the spins injected to the network
 |              at time t = 0 to reach the output node. When a molecule reaches the output node,
 |              NaNs are stored in all following time points (the spin has left the network).
 |              - vs[p][0,k] stores the x-component of the velocity of a molecule following the
 |                p-th path at the k-th time step
 |              - vs[p][1,k] stores the y-component of the velocity of a molecule following the
 |                p-th path at the k-th time step
 |              - vs[p][2,k] stores the z-component of the velocity of a molecule following the
 |                p-th path at the k-th time step
 |      * qs:   list storing the volumetric flow rate over time (units: mm3/s) experienced by a
 |              spin while it travels along all possible P paths that connect the input node to
 |              the output node. The p-th element qs[p], for p = 0, ..., P-1
 |              (where P = len(obj.iopaths)), is a 1D numpy array that stores
 |              the flow experienced by a spin along the p-th paths connecting the input node to the
 |              output node. It has a length equal to Ntpoints elements, where Ntpoins is the
 |              number of time points of duration dt required for all spins injected to the network
 |              at time t = 0 to reach the output node. When a molecule reaches the output node,
 |              NaNs are stored in all following time points (the spin has left the network).
 |              - qs[p][k] stores the instantaneous flow rate experienced by a molecule following the
 |                p-th path at the k-th time step
 |  
 |  GetTrajUniformSeed(self, Ntpoints, Nspins=2000, dt=2e-05, seednumber=20181019, boundary='periodic')
 |      Calculates spin trajectories for a vascular network and uniformly-seeded spins
 |      
 |      rmat = obj.GetTrajUniformSeed(<parameters>)
 |      
 |      MANDATORY INPUT PARAMETERS
 |      * Ntpoints: number of time points required (number of discrete time steps at which
 |                  spin trajectories will be evaluated). The total duration of the simulation
 |                  will be (Ntpoints-1)*dt, with dt the temporal resolution expressed in s.
 |                  Note that the default dt is 20 us (see below).
 |      
 |      OPTIONAL INPUT PARAMETERS
 |      * Nspins:   number of spins to simulate (default: 2000).
 |      * dt:       temporal resolution of the simulation (scalar; units: s; default: 2e-5,
 |                  i.e., 20 us)
 |      *seednumber: Random seed used for spin position initialisation in the network. Spin positions
 |                   are initialised uniformly within each segment (pipe) of the network. The number
 |                   of spins placed within each segment its proportional to the volume of the
 |                   segment, with respecto to the total volume of the network. Default seed: 20181019
 |      * boundary: string specifying the boundary condition for particles that reach the output node.
 |                  It can be set to "periodic" or "feedback". Default: "periodic".
 |                  - "periodic" -> a particle that reaches the output node starts again
 |                    a new journey through a "copy" of the network that is "attached" to the output node.
 |                    This effectively makes the network periodic.
 |                  - "feedback" -> a spin that reaches the output node is fed back to the input node,
 |                    and starts a new journey through the vascular network. Note that
 |                    in this case the particle positions and velocities over time will exhibit a
 |                    discontinuity when the "jump" from output node to input node occurs
 |      
 |      RETURNS
 |      * rmat:     matrix storing the trajectories of the spins (units: mm). It has size
 |                  3 x Ntsteps x Nspins , where Nspins is the number of requested spins,
 |                  and Ntsetps is the number of simulation time points at a temporal resolution equal to dt.
 |                  Note that for conservation of total mass the vascular network a boundary condition is
 |                  applied depending on the value of the input "boundary". This means that spins that reach
 |                  the output node before the simulation is finshed, start a new journey through the network
 |                  i) in a "copy" of the network that is physically attached to the output node
 |                     (if boundary is set to "periodic"), or
 |                  ii) after having "jumped" back to the input node
 |                  (if boundary is set to "feedback). 
 |                  A detailed description of rmat is:
 |                  - rmat[0,k,s] -> x-component of the trajectory of the s-th spin, at time step k-th
 |                  - rmat[1,k,s] -> y-component of the trajectory of the s-th spin, at time step k-th
 |                  - rmat[2,k,s] -> z-component of the trajectory of the s-th spin, at time step k-th
 |  
 |  __init__(self, nodes, radii, qin, idxin, idxout, radiusin=None, visc=1.2, flowmodel='blinder', solver='numerical')
 |      pipenet object constructor -- type help(syn.pipenet) for the help manual
 |  
 |  dMRISynMea(self, bval, Gdur, Gsep, Gdir, Nspins=2000, dt=2e-05, Gx=None, Gy=None, Gz=None, rndseed=20181019, bcon='periodic')
 |      Synthesises a diffusion-weighted MRI measurement from protons flowing within a vascular network
 |      
 |      magtot, phasetot, phasespin, rmat, tarray, grx, gry, grz = obj.dMRISynMea(<parameters>)
 |      
 |      MANDATORY INPUT PARAMETERS
 |      * bval: b-value for a PGSE sequence (scalar; units: s/mm2);
 |      * Gdur: gradient duration (small delta) for a PGSE sequence (scalar; units: s); if optional 
 |              Gx, Gy, Gz gradient waveforms are provided, this parameter is ignored;
 |      * Gsep: gradient separation (large Delta) for a PGSE sequence (scalar; units: s); if optional
 |              Gx, Gy, Gz gradient waveforms are provided, this parameter is ignored
 |      * Gdir: gradient direction for a PGSE sequence; it must be a 1D numpy array of size 3, storing the
 |              x, y and z components of the gradient direction; ; if optional
 |              Gx, Gy, Gz gradient waveforms are provided, this parameter is ignored
 |              
 |      **** NOTE THAT ALL MANDATORY INPUT PARAMETERS ABOUT b-value, GRADIENT DURATION/SEPARATION/DIRECTION
 |      ARE IGNORED IF CUSTOM GRADIENT WAVEFORMS ARE PROVIDED VIA OPTIONAL INPUT PARAMETERS Gx, Gy, Gz                
 |      
 |      OPTIONAL INPUT PARAMETERS
 |      * Nspins: number of spins to simulate (default: 2000). Note that in practice a slighter smaller
 |                number of particles will be simulated, as some particles will leave the vascular network
 |                during the simulation and will not be used to synthesise the MRI signal, as in
 |                Phi Van V et al, NMR in Biomedicine 2021, 34(7):e4528, doi: 10.1002/nbm.4528
 |      * dt:     temporal resolution of the simulation (scalar; units: s; default: 2e-5, i.e., 20 us)
 |      * Gx:     waveform for x-component of the diffusion encoding gradient, expressed in T/mm.
 |                If specified, previous parameters bval, Gdur, Gsep, Gdir will be ignored.
 |                Gx must be a 1D array storing the values of Gx over time (Gx(t=0), Gx(t=dt), Gx(t=2*dt), ....).
 |                If Gx is specified, Gy and Gz must be specified too. Gx, Gy and Gz must have the same length.
 |                Default value for Gx: None
 |      * Gy:     waveform for y-component of the diffusion encoding gradient, expressed in T/mm.
 |                If specified, previous parameters bval, Gdur, Gsep, Gdir will be ignored.
 |                Gy must be a 1D array storing the values of Gy over time (Gy(t=0), Gy(t=dt), Gy(t=2*dt), ....).
 |                If Gy is specified, Gx and Gz must be specified too. Gx, Gy and Gz must have the same length.
 |                Default value for Gy: None
 |      * Gz:     waveform for z-component of the diffusion encoding gradient, expressed in T/mm.
 |                If specified, previous parameters bval, Gdur, Gsep, Gdir will be ignored.
 |                Gz must be a 1D array storing the values of Gz over time (Gz(t=0), Gz(t=dt), Gz(t=2*dt), ....).
 |                If Gz is specified, Gx and Gy must be specified too. Gx, Gy and Gz must have the same length.
 |                Default value for Gz: None
 |      *rndseed: Random seed used for spin position initialisation in the network. Spin initialisation is
 |                performed by sampling at uniform spatial locations consecutive waves of spins that flow
 |                from the input node to the output node, making sure that the number of spins in each pipe at
 |                t = 0 is proportional to the flow in that pipe (conservation of mass).
 |                Default: 20181019
 |      * bcon:   string specifying the boundary condition for particles that reach the output node.
 |                It can be set to "periodic" or "feedback". Default: "periodic".
 |                - "periodic" -> a particle that reaches the output node starts again
 |                  its trajectory, but starting from the output node, rather than from the input
 |                  node. In practice this makes the vascular network periodic, and a "copy" of the
 |                  network is "attached" to the output node. A particle reaching the output node
 |                  will continue travelling through this "copy" of the network.
 |                - "feedback" -> a spin that reaches the output node is fed back to the input node,
 |                  and it will follow again the same trajectory it had taken before. Note that
 |                  in this case the particle positions and velocities over time will exhibit a
 |                  discontinuity when the "jump" from output node to input node occurs. This could cause
 |                  some unexpected signal behaviours: for example, the signal from a
 |                  straight pipe (where all spins travel in parallel at the same velocity) will show
 |                  magnitude MRI signal attenuation, while in theory should show no signal
 |                  attenuation but only a phase shift. Also, the signal decay could be overstimated when
 |                  the diffusion encoding gradient has a strong component along the direction of
 |                  such a "jump", and a considerable number of particles "jump" during the simulation.
 |      
 |      RETURNS
 |      * magtot:    magnitude of the synthesised diffusion MRI signal (note that in absence of diffusion-weighting,
 |                   a magnitude of 1.0 will be observed)
 |      * phasetot:  phase of the synthesised diffusion MRI signal, defined in the range [-pi; pi]
 |      * phasespin: 1D array storing the phase accrual of all individual spins. The array has size Nspins,
 |                   so that phasespin[s] stores the phase accrual of the s-th spin, for s = 0, ..., Nspins - 1.
 |                   It ranges between ]-inf; +inf[ so that one can study the number of phase wraps caused
 |                   by diffusion-weighting
 |      * rmat:     matrix storing the trajectories of the spins (units: mm). It has size
 |                  3 x Ntsteps x Nspins , where Nspins is the number of requested spins,
 |                  and Ntsetps is the number of simulation time points at a temporal resolution equal to dt.
 |                  Note that for conservation of total mass the vascular network a boundary condition is
 |                  applied depending on the value of the input "bcon". This means that spins that reach
 |                  the output node before the simulation is finshed, can either i) start again their
 |                  trajectory in a "copy" of the network attached to the output node
 |                  (if bcon is set to "periodic"), or ii) be fed back to the input node
 |                  (if bcon is set to "feedback). A detailed description of rmat is
 |                  - rmat[0,k,s] -> x-component of the trajectory of the s-th spin, at time step k-th
 |                  - rmat[1,k,s] -> y-component of the trajectory of the s-th spin, at time step k-th
 |                  - rmat[2,k,s] -> z-component of the trajectory of the s-th spin, at time step k-th
 |      * tarray:   time array (units: s) of length Ntsteps. tarray[k] = k*dt stores the time corresponding to
 |                  the k-th simulation step
 |      * grx:      waveform for the x-component of the diffusion encoding gradient (Gx), expressed in T/mm. It has size
 |                  Ntsteps, and grx[k] is the value of Gx at time instant tarray[k] = k*dt. If optional input Gx was
 |                  provided, than grx is just a copy of Gx
 |      * gry:      waveform for the y-component of the diffusion encoding gradient (Gy), expressed in T/mm. It has size
 |                  Ntsteps, and gry[k] is the value of Gy at time instant tarray[k] = k*dt. If optional input Gy was
 |                  provided, than gry is just a copy of Gy
 |      * grz:      waveform for the z-component of the diffusion encoding gradient (Gz), expressed in T/mm. It has size
 |                  Ntsteps, and grz[k] is the value of Gz at time instant tarray[k] = k*dt. If optional input Gz was
 |                  provided, than grz is just a copy of Gz
 |  
 |  dMRISynProt(self, bvarr, Gdurarr, Gseparr, Gdirarr, Nwater=2000, deltat=2e-05, Gradx=None, Grady=None, Gradz=None, Nrep=1, seed=20181019, boundcon='periodic', vb=True)
 |      Synthesises a diffusion-weighted MRI protocol made of multiple measurements from
 |      protons flowing within a vascular network
 |      
 |      mag, phase = obj.dMRISynProt(<parameters>)
 |      
 |      MANDATORY INPUT PARAMETERS
 |      * bvarr: array of b-values for a PGSE sequence (units: s/mm2; it must be a 1D array of size Nmeas)
 |      * Gdurarr: array of gradient durations (small delta) for a PGSE sequence (1D array of size Nmeas; units: s)
 |      * Gseparr: array of gradient separations (large Delta) for a PGSE sequence (1D array of size Nmeas; units: s)
 |      * Gdirarr: array of gradient directions for a PGSE sequence; it must be a 2D numpy array of size Nmeasx3,
 |                 where Gdirarr[m,0] stores the x-component of the gradient direction corresponding to the m-th
 |                 measurement; Gdirarr[m,1] the y-component; Gdirarr[m,2] the z-component
 |      
 |      ***** NOTE THAT ALL INPUT PARAMETERS ABOUT b-value, GRADIENT DIRECTION/DURATION/SEPARATION ARE IGNORED
 |            IF CUSTOM GRADIENT WAVEFORMS ARE PROVIDED WITH OPTION INPUT PARAMETERS Gradx, Grady, Gradz
 |            (SEE BELOW) 
 |      
 |      OPTIONAL INPUT PARAMETERS
 |      * Nwater: number of spins to simulate (default: 2000). Note that in practice a slighter smaller
 |                number of particles will be simulated, as some particles will leave the vascular network
 |                during the simulation and will not be used to synthesise the MRI signal, as in
 |                Phi Van V et al, NMR in Biomedicine 2021, 34(7):e4528, doi: 10.1002/nbm.4528
 |      * deltat: temporal resolution of the simulation (scalar; units: s; default: 2e-5, i.e., 20 us)
 |      * Gradx:  waveforms for x-components (Gx) of the diffusion encoding gradients, expressed in T/mm.
 |                If specified, previous parameters bvarr, Gdurarr, Gseparr, Gdirarr will be ignored.
 |                Gradx must be a 2D array storing the values of Gx over time (Gx(t=0), Gx(t=dt), Gx(t=2*dt), ....).
 |                for all measurements m = 0, ..., Nmeas - 1.
 |                Gx(m,k) stores Gx for the m-th protocol measurement at the k-th time step.
 |                If Gradx is specified, Grady and Gradz must be specified too.
 |                Gradx, Grady and Gradz must have the same size. Default value for Gradx: None
 |      * Grady:  waveforms for y-components (Gy) of the diffusion encoding gradients, expressed in T/mm.
 |                If specified, previous parameters bvarr, Gdurarr, Gseparr, Gdirarr will be ignored.
 |                Grady must be a 2D array storing the values of Gy over time (Gy(t=0), Gy(t=dt), Gy(t=2*dt), ....).
 |                for all measurements m = 0, ..., Nmeas - 1.
 |                Gy(m,k) stores Gy for the m-th protocol measurement at the k-th time step.
 |                If Grady is specified, Gradx and Gradz must be specified too.
 |                Gradx, Grady and Gradz must have the same size. Default value for Grady: None
 |      * Gradz:  waveforms for z-components (Gz) of the diffusion encoding gradients, expressed in T/mm.
 |                If specified, previous parameters bvarr, Gdurarr, Gseparr, Gdirarr will be ignored.
 |                Gradz must be a 2D array storing the values of Gz over time (Gz(t=0), Gz(t=dt), Gz(t=2*dt), ....).
 |                for all measurements m = 0, ..., Nmeas - 1.
 |                Gz(m,k) stores Gz for the m-th protocol measurement at the k-th time step.
 |                If Gradz is specified, Gradx and Grady must be specified too.
 |                Gradx, Grady and Gradz must have the same size. Default value for Gradz: None
 |      * Nrep:   number of repetitions of the simulation. For each repetition, a different seed number is used,
 |                resulting in slightly different initial spin positions, and hence MRI signals.
 |                Note that the calculation of the spin trajectories (the most costly computational step) is
 |                done once per each simulation repetition, and then the calculated trajectories are used
 |                to synthesise signals for all the diffusion MRI measurements of the protocol.
 |      * seed:   Random seed used for spin position initialisation in the network. Spin initialisation is
 |                performed by sampling at uniform spatial locations consecutive waves of spins that flow
 |                from the input node to the output node, making sure that the number of spins in each pipe at
 |                t = 0 is proportional to the flow in that pipe (conservation of mass).
 |                Default: 20181019
 |      * boundcon: string specifying the boundary condition for particles that reach the output node.
 |                It can be set to "periodic" or "feedback". Default: "periodic".
 |                - "periodic" -> a particle that reaches the output node starts again
 |                  its trajectory, but starting from the output node, rather than from the input
 |                  node. In practice this makes the vascular network periodic, and a "copy" of the
 |                  network is "attached" to the output node. A particle reaching the output node
 |                  will continue travelling through this "copy" of the network.
 |                - "feedback" -> a spin that reaches the output node is fed back to the input node,
 |                  and it will follow again the same trajectory it had taken before. Note that
 |                  in this case the particle positions and velocities over time will exhibit a
 |                  discontinuity when the "jump" from output node to input node occurs. This could cause
 |                  some unexpected signal behaviours: for example, the signal from a
 |                  straight pipe (where all spins travel in parallel at the same velocity) will show
 |                  magnitude MRI signal attenuation, while in theory should show no signal
 |                  attenuation but only a phase shift. Also, the signal decay could be overstimated when
 |                  the diffusion encoding gradient has a strong component along the direction of
 |                  such a "jump", and a considerable number of particles "jump" during the simulation.
 |      * vb:       verbose. Set this parameter to True (or 1) if you want feedback on the simulation
 |                  duration to be printed on the standard output. Any other value will be trated as
 |                  a False and no feedback will be printed. Default: True
 |      
 |      RETURNS
 |      * mag:    1D or 2D array storing the magnitude of the synthesised MRI signals given the input protocol.
 |                If Nrep = 1, then mag will be a 1D array of length Nmeas, with mag[m] storing the
 |                magnitude for the m-th protocol measurement. If Nrep > 1, then mag will be a 2D numpy array
 |                of size Nrep x Nmeas. In that case, mag[r,m] will store the magnitude for the m-th protocol
 |                measurement obtained during the r-th simulation repetition (r = 0, ..., Nrep - 1).
 |      * phase:  1D or 2D array storing the phase of the synthesised MRI signals given the input protocol.
 |                If Nrep = 1, then phase will be a 1D array of length Nmeas, with phase[m] storing the
 |                phase for the m-th protocol measurement. If Nrep > 1, then phase will be a 2D numpy array
 |                of size Nrep x Nmeas. In that case, phase[r,m] will store the phase for the m-th protocol
 |                measurement obtained during the r-th simulation repetition (r = 0, ..., Nrep - 1).
 |                phase is defined in [-pi; pi].
 |  
 |  ----------------------------------------------------------------------
 |  Data descriptors defined here:
 |  
 |  __dict__
 |      dictionary for instance variables (if defined)
 |  
 |  __weakref__
 |      list of weak references to the object (if defined)
```  
