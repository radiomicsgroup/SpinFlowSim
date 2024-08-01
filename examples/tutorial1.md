# Tutorial 1: resolving vascular networks for dMRI signal synthesis

This tutorial shows you how to create vascular networks, and how to use them to synthesise diffusion MRI (dMRI) signals. To go through our tutorials, remember that you need:
- [numpy](https://numpy.org) (developed with version 1.24.2)
- [pandas](https://pandas.pydata.org) (developed with version 1.5.3)
- [Lcapy](https://lcapy.readthedocs.io) (developed with version 1.10)
- [PySpice](https://github.com/FabriceSalvaire/PySpice) (developed with version 1.5)
- [graph-tool](https://graph-tool.skewed.de) (developed with version 2.45, commit b1a649d8)
- [nibabel](https://nipy.org/nibabel) (developed with version '5.1.0')
- [matplotlib](https://matplotlib.org) (examples tested with version '3.7.1')
- [mrtrix](https://www.mrtrix.org) (developed with version '3.0.4').

The tutorial contains the following sections:
* [The _pipenet_ class](#pipenet-class)
* [Generating a new network from histology](#drawing-net)
* [Loading our pre-computed vascular networks](#pre-computed-network)
* [Plotting properties of a vascular network](#plot-network)
* [Synthesising vascular diffusion MRI signals](#dMRI-signals)
* [Final remarks](#remarks)


## The _pipenet_ class <a name="pipenet-class"></a>

We have defined the _pipenet_ class to work with vascular networks. In our framework, a vascular network is nothing else but a collection _pipes_, representing capillaries, which connect a set of _nodes_.  

To instantiate a new object from the _pipenet_ class you need the following mandatory input parameters:
* **a matrix of node positions (in mm)**, of size 3 x Nnodes (rows: xpos, ypos, zpos in mm; a node is the input (or the output) of a pipe);
* **a matrix of pipe radii (in mm)**, of size Nnodes x Nnodes (element (i,j) stores the radius of the pipe connecting node i with node j);
* **the input volumetric flow rate** in mm<sup>3</sup>/s;
* **the index of the input node** in the node matrix (nodes[:,idxin] provides the x,y,z coordinates of such a input node);
* **the index of the output node** in the node matrix (nodes[:,idxout] provides the x,y,z coordinates of such an input node).

Additional optional parameters are the **fluid viscosity**, the **radius at the inlet**, the **fluid dynamics model** and the **type of solver**. The full manual of the _pipenet_ class can be read [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/pipenet_manual.md).

For example, this code creates an initialises a _pipenet_ object describing a simple 2D 3-capillary network, made of 3 nodes, namely nodes 0, 1 and 2, as illustrated below. Connections are between nodes 0 and 1, between 0 and 2, and between 1 and 2; the input flow, of 0.0055 mm<sup>3</sup>/s, comes from node 0, and the output node is 2.

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/networkexample.png" width="700" height="auto">
</div>

```
import numpy as np
import sys
sys.path.insert(0, '../code' )      # Add the SpinFlowSim "code" folder where the syn.py and visu.py files are stored
import syn                          # import syn.py, where the pipenet class is defined

nodepos = np.transpose( np.array([ [0.0, 0.0, 0.0], [0.010, 0.020, 0.0], [0.030, 0.0, 0.0] ] ) )    # Node positions in mm

r = np.zeros((3,3))           # Node positions in mm
r[0][0] = r[1][1] = r[2][2] = np.nan
r[0][1] = r[1][0] = 0.005
r[0][2] = r[2][0] = 0.007
r[1][2] = r[2][1] = 0.008

qin = 0.0055                 # Input volumetric flow rate in mm3/s
inputnode = 0                # Inlet of the whole network
outputnode = 2               # Outlet of the whole network

muval = 1.2                  # Define the viscosity of pure plasma in mPa x s
bloodmodel = 'blinder'       # Define the flow model
solvtechnique = 'numerical'  # Solve the pipe circuit numerically using PySpice
net = syn.pipenet(nodepos,r,qin,inputnode,outputnode,visc=muval,flowmodel=bloodmodel,solver=solvtechnique)

```
For the network initialisiation, we ensure that the fluid dynamics follows the model used in Blinder et al (Blinder P et al, Nature Neuroscience volume 16, pages 889–897 (2013), doi: [10.1038/nn.3426](https://doi.org/10.1038/nn.3426); optional parameter `flowmodel=bloodmodel`) with a cell-free plasma viscosity of 1.2 mPa s (`visc=muval`). The network is solved numerically (`solver=solvtechnique`), with default inlet radius computation (mean of the radii emanating from the inlet).

Upon declaration of the _pipenet_ object, the vascular network is solved, meaning that the volumetric flow rate (VFR) is calculated in each capillary segment. This calculation is done by default numerically (`solver='numerical'` parameter in the _pipenet_ constructor), using [PySpice](https://github.com/PySpice-org/PySpice). Exact analytical VFR estimation can also be performed with [Lcapy](https://lcapy.readthedocs.io), setting `solver='symbolic'`. Be aware though that in this case the computation time is much longer, and practically unfeasible when the number of capillaries is higher than approximately 20-25. 

Once we have initialised a _pipenet_ object, we have solved the network, and we have hence reconstructed the VFR (in mm<sup>3</sup>/s) in each capillary. Such VFR values are stored in the `flowmat` attribute,
```
>>> print(net.flowmat)
[[        nan  0.00138352  0.00411648]
 [-0.00138352         nan  0.00138352]
 [-0.00411648 -0.00138352         nan]]
```
while the corresponding blood velocity (BV) values (in mm/s) are stored in the `velmat` attribute: 
```
>>> print(net.velmat)
[[         nan  17.61550894  26.74115613]
 [-17.61550894          nan   6.88105818]
 [-26.74115613  -6.88105818          nan]]
```
We can also check, for example, the set of paths connecting inlet/outlet, represented as lists of consecutive nodes:
```
>>> print(net.iopaths)
[array([0, 1, 2], dtype=uint64), array([0, 2], dtype=uint64)]
```
These have been evaluated with [graph-tools](https://graph-tool.skewed.de/) using graph theory, since our networks are essentially directed graphs. In our simple 3-capillary network example, we can see that there are only two possible paths from inlet to outlet.

A final note: elements `flowmat[i][j]` and `velmat[i][j]` respectively store the VFR and BV of the capillary **going from node _i_ to node _j_. Their sign is > 0 if the flow goes from node _i_ to node _j_, while it is < 0 if the flow goes from node _j_ to node _i_**. By definition then, it follows that `flowmat[i][j]` and `flowmat[j][i]` (and, similarly, `velmat[i][j]` and `velmat[j][i]`) have opposite signs.  

[Here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/pipenet_manual.md) you can find a complete description of all methods and attributes of the _pipenet_ class.


## Generating a new network from histology <a name="drawing-net"></a>
With SpinFlowSim, one can simulate flow in realistic capillary networks that have been carefully reconstructed from histology. Here we show an example of a Hematoxylin and Eosin (HE) image of a mouse kidney slice, where a vascular network has been segmented by tracing visible capillaries. The network is made of 26 nodes, connected among each other through 37 straight capillaries. 

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/kidneynet.png" width="750" height="auto">
</div>

We have stored the information realted to the network as a CSV spreadsheet. Each straight capillary is described by a row in [Network_stats.csv](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/Network_stats.csv). Each row is described by the following variables:
- segmented ID (SegmentNumber)
- mean radius (RadiiMean_um) in µm; computed by averaging three radii measured over the length of a segment
- starting coordinates (Xstart,Ystart,Zstart) in µm
- ending coordinates (Xend,Yend,Zend) in µm
- starting node ID (Node_at_Xstart_Ystart)
- ending node ID (Node_at_Xend_Yend)

The table below illustrates the first few lines of the CSV file.


| SegmentNumber | RadiiMean_um      | Xstart      | Ystart     | Zstart | Xend        | Yend       | Zend | Node_at_Xstart_Ystart | Node_at_Xend_Yend |
|---------------|-------------------|-------------|------------|--------|-------------|------------|------|-----------------------|-------------------|
| 1             | 3.3177            | 10056.36    | 10111.85   | 0      | 10083.96    | 10073.92   | 0    | 0                     | 1                 |
| 2             | 2.4370833333333333| 10083.96    | 10073.92   | 0      | 10182.86    | 10164.91   | 0    | 1                     | 2                 |
| 3             | 2.2855166666666666| 10182.86    | 10164.91   | 0      | 10088.94    | 10136.81   | 0    | 2                     | 3                 |
| 4             | 3.306783333333333 | 10056.36    | 10111.85   | 0      | 10088.94    | 10136.81   | 0    | 0                     | 3                 |
....


We are including this CSV file within the directory [**examples**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples) (file [_Network_stats.csv_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/Network_stats.csv)). The script [_script01_initnet.py_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/script01_initnet.py), also included in the folder [**examples**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples), shows how to load the network information contained in the CSV with [_pandas_](https://pandas.pydata.org), so that it can be used to initialise a _pipenet_ object. In the script, we assume, as an example, that the inlet/outlet of the network are respectively nodes 0/16.


## Loading our pre-computed vascular networks <a name="pre-computed-network"></a>

In the folder [**networks**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks) we include different instantiations of 15 vascular networks that we have drawn on histological images of human liver biopsies. We are including 100 instantiations per network, by varying the inlet/outlet and input VFR, for a total of 1500 unique nets. [Here](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks/README.md) you can find a detailed description of these vascular networks and of their corresponding flow properties. 

All these 1500 nets, already resolved for you, have been stored as binary files. Loading one of such files into python gives you access to an initialised _pipenet_ object, which you can use to synthesise vascular dMRI signals. 

In this code below we show, for example, how to load one of such networks, and how to plot its VFR matrix (attribute `flowmat`). We assume that the code is run from inside the folder [**examples**](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples) where this tutorial is stored, and that matplotlibt is installed in your _python_ environment.

```
import numpy as np
import sys
from matplotlib import pyplot as plt
import pickle as pk
sys.path.insert(0, '../code' )      # Add the SpinFlowSim "code" folder where the syn.py and visu.py files are stored
import syn

# Load network Net10, input node 0, output node 28, input VFR 0.0007 mm3/s
h = open('../networks/Net10/net10_Nin0_Nout28_Qin0.0007.bin','rb')
net = pk.load(h)
h.close()

# Plot the Volumetric Flow Rate matrix
plt.imshow(net.flowmat,aspect='auto')
plt.title('VFR matrix in mm$^3$/s')
plt.xlabel('Node index')
plt.ylabel('Node index')
plt.colorbar()
plt.show()

```

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/connectivityexample.png" width="450" height="auto">
</div>


## Plotting properties of a vascular network <a name="plot-network"></a>
We now illustrate how to visualise some useful properties of a resolved vascular network. Let's assume, for example, that we want to visualise results from the exemplificatory kidney vascular network that we derived from histology above. 

The code below shows how to load the vascular network, previously saved as a binary file by script [_script01_initnet.py_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/script01_initnet.py). We then generate illustrative trajectories for 100 time steps of 750 spins, using a temporal resolution of 10<sup>-5</sup> s (i.e., 10 µs). We also use a seed number to initialise spin positions (`seednumber=20181019`), for reproducibility purposes, and adopt a _periodic_ boundary condition (option `boundary='periodic`, which is the default). This means that if a spin reaches the network output node, it will continue its trajectory on a "virtual copy" of the network itself, whose input node coincides exactly with the output node.

We can visualise the trajectories with the [_visu.py_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/visu.py) module (whose manual can be found [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/visu_manual.md)): the module contains the `spin_animation()` function, which allows you to show a video of the flowing spins, and to save it as a GIF file.

```
import numpy as np
from matplotlib import pyplot as plt
import sys
import pickle as pk
sys.path.insert(0, '../code' )      # Add the SpinFlowSim folder where the syn.py and visu.py files are stored
import syn
import visu

### Load the kidney vascular networks created with previous script script01_initnet.py
h = open('./script01_initnet_kidneyexample.bin','rb')
net = pk.load(h)
h.close()

### Generate spins trajectories: 750 spins for 100 time steps of dt = 10 us; boundary condition 'periodic'  
rperiodic = net.GetTrajUniformSeed(100, Nspins=750, dt=1e-05, seednumber=20181019, boundary='periodic' )

### Let's save them as a GIF
visu.spin_animation(rperiodic, 'spinvideo_periodic.gif')
```

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/spinvideo_periodic.gif" width="1000" height="auto" >
</div>



Another avilable boundary condition is the _feedback_ one (option `boundary='feedback'`). In this case, spins reaching the outlet are fed back instantaneously to the inlet. Note that this introduces "jumps" on the spin trajectories, which may affect the synthesis of your dMRI signals - so use it with care. Below are the same spin trajectories of above, but with the "feedback" condition.

```
### Generate spins trajectories: again, 750 spins, 100 time steps, dt = 10 us, but boundary condition 'feedback'  
rfeedback = net.GetTrajUniformSeed(100, Nspins=750, dt=1e-05, seednumber=20181019, boundary='feedback' )

### Let's save them as a GIF
visu.spin_animation(rfeedback, 'spinvideo_feedback.gif')
```

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/spinvideo_feedback.gif" width="1000" height="auto" >
</div>


From the GIF animations, it is apparent that spins flowing through different capillaries experience very different velocities. Below we show how to generate a plot in which we colour each capillary of the network according to the velocity of the spins flowing through the capillary itself:

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/kidneynet_velocityplot.png" width="1000" height="auto" >
</div>

```
import pickle as pk
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.colors import Normalize

# Load initialized vascular network
infile = 'script01_initnet_kidneyexample.bin'
h = open(infile,'rb')
net = pk.load(h)
h.close()

# Retrieve network properties
v = net.velmat  # Get velocity matrix
connectivity = net.connmat # Get connectivity matrix
velmax = np.nanmax( np.abs( np.unique(net.velmat) ) )
velmin = 0.0

# Create plot
fig, ax = plt.subplots(figsize=(14, 12))
# Loop over each connection
for i in range(connectivity.shape[0]):
    for j in range(connectivity.shape[1]):
        if connectivity[i, j] > 0:  # Exclude zero or negative values
            velocity = v[i, j]
            if velocity > 0:
                x1, y1, z1 = net.nodepos[:, i]
                x2, y2, z2 = net.nodepos[:, j]
                # Plot each line with color based on velocity
                color = plt.cm.viridis(velocity / velmax)
                ax.plot([x1, x2], [y1, y2], color=color,linewidth=4)


# Add colorbar
norm = Normalize(vmin=velmin, vmax=velmax)
cmap = plt.set_cmap('viridis')
sm = plt.cm.ScalarMappable(cmap=cmap, norm=norm)
sm.set_array(v)
cbar = plt.colorbar(sm, ax=ax, pad=0.015, shrink=0.85)
cbar.set_label('$[mm/s]$', rotation=90, fontsize=11)
cbar.set_ticks([velmin, velmax])
cbar.ax.set_title('$v$', fontsize=15)

plt.axis('equal')
plt.xlabel('x-position [mm]')
plt.ylabel('y-position [mm]')
plt.title('Spin paths coloured by velocity')
plt.show()
```


## Synthesising vascular diffusion MRI signals <a name="dMRI-signals"></a>
Finally, once we have drawn a vascular network on histology, and used such data to initialise a _pipenet_ vascular network obkect, we can use it to generate realistic vascular dMRI signals! Continuing the example above, we now show you how to generate signals from the kidney vascular network initialised by previous script [_script01_initnet.py_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/script01_initnet.py).

Let's assume that we want to generate the complex-valued dMRI signal for a standard Pulsed Gradient Spin Echo (PGSE) sequence. We are interested in a rich set of b-values, ranging in [0; 100] s/mm<sup>2</sup> and three diffusion times, e.g., Δ = {20, 50, 80} ms, with δ fixed to 15 ms. We can use the method `dMRISynProt()` of a _pipenet_ object.
```
import numpy as np
from matplotlib import pyplot as plt
import sys
import pickle as pk
sys.path.insert(0, '../code' )      # Add the SpinFlowSim folder where the syn.py and visu.py files are stored
import syn
import visu

## Load previously resolved vascular network
h = open('./script01_initnet_kidneyexample.bin','rb')
net = pk.load(h)
h.close()

### Simulate two diffusion times for 100 b-values between 0 and 100 s/mm2. 
# We will get magnitude and phase for 3000 spins and a temporal resolution of 30 us 
# We use the 'periodic' boundary condition, and fix the seed number of reproducibility
# We will also use a fixed gradient -- g = [1 0 0] for all b-values, as an example
# We fix the gradient duration to 15 ms, and try two different gradient separations
Nspins=3000 
tempres=3e-05
seednumber=100287
mycond='periodic'
bval = np.linspace(0,100,100)   
gdir = np.zeros((bval.size,3))
gdir[:,0] = 1.0
gdir[:,1] = 0.0
gdir[:,2] = 0.0
gdir[0,0] = 0.0     # Make sure b = 0 as gradient direction [0 0 0]
gdur = 0.015*np.ones(bval.shape)  
gdur[bval==0] = 0

### Get magnitue and phase for a gradient separation of 20 ms
gsep = 0.020*np.ones(bval.shape)      
gsep[bval==0] = 0
smag1,sph1 = net.dMRISynProt(bval, gdur, gsep, gdir, Nwater=Nspins, deltat=tempres, seed=seednumber, boundcon=mycond)

### Get magnitue and phase for a gradient separation of 70 ms
gsep = 0.070*np.ones(bval.shape)      
gsep[bval==0] = 0
smag2,sph2 = net.dMRISynProt(bval, gdur, gsep, gdir, Nwater=Nspins, deltat=tempres, seed=seednumber, boundcon=mycond)

### Let's plot magnitude
plt.subplot(1,2,1)
plt.plot(bval,smag1,label='δ = 15ms; Δ = 20 ms')
plt.plot(bval,smag2,label='δ = 15ms; Δ = 70 ms')
plt.yscale('log')
plt.xlabel('b-values [s/mm$^2$]')
plt.ylabel('Signal magnitude')
plt.legend()

### Let's plot phase in deg
plt.subplot(1,2,2)
plt.plot(bval,sph1*180.0/np.pi,label='δ = 15ms; Δ = 20 ms')
plt.plot(bval,sph2*180.0/np.pi,label='δ = 15ms; Δ = 70 ms')
plt.xlabel('b-values [s/mm$^2$]')
plt.ylabel('Signal phase [deg]')
plt.legend()

plt.show()
```

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/signal_example.png" width="1000" height="auto" >
</div>


The plots reveal several interesting features, e.g., non-mono-exponential decay of the signal magnitude with oscillatory patterns, characteristic of spins whose flow regime is ballistic. We also observe wraps in the signal phase as a function of _b_, seen more clearly at long diffusion time.


## Final remarks <a name="remarks"></a>
This tutorial has shown you how to use the _pipenet_ class defined in module [_syn_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/syn.py) to resolve vascular networks and synthesise dMRI signals from them. We have taken you through some basic usage; however, be aware that SpinFlowSim can do much more! 

* The methods [`dMRISynMea()`](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/pipenet_manual.md) and [`dMRISynProt()`](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/pipenet_manual.md) of a _pipenet_ object can synthesise dMRI signal synthesis for **any, arbitary, input custom gradient waveform**. In this tutorial we have shown signals from standard PGSE, but more advanced diffusion encodings can also be investigated (e.g., b-tensor encoding, Oscillating Gradient Spin Echo (OGSE), flow-compensated diffusion encoding, etc).
* The [`syn`](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/syn.py) module includes functions outside of the _pipenet_ class. For example:
  * `SegmentIn2Out()` calculates the trajectory of a spin flowing at a given velocity along a straight line connecting two points;
  * `ComputeResHagPois()` calculates flow resistances with the Hagen-Poiseuille law;
  * `ComputeResModHagPois()` calculates flow resistances with the modified Hagen-Poiseuille law as in Blinder et al, Nat Neurosci 2013, 16(7): 889-897, doi: [10.1038/nn.3426](https://doi.org/10.1038/nn.3426);
  * `ComputeFlow()` calculates the volumetric flow rate between each pair of 2 nodes in a vascular network;
  * `ComputeAllPaths()` finds all paths connecting the input node to the output node of a network given its connecivity matrix. It returns a directed graph given the connectivity matrix;
  * `Traj2Signal()` synthesises a diffusion-weighted MRI signal from a set of spin trajectories (not necessarily from flowing spins, but even for pure diffusion!);
  * `getGradPGSE()` generates diffusion encoding gradient wave forms for an ideal PGSE sequence.
