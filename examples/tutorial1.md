# Initialising _pipenet_ objects to resolve a vascular network


## The _pipenet_ class

We have defined the _pipenet_ to work with vascular networks, which can later be used to simulated diffusion MRI (dMRI) signals for any given protocol of interest. In our framework, a vascular network is nothing else but a collection _pipes_, representing capillaries, which conned a set of _nodes_.  

To create an object from the _pipenet_ class ones needs the following mandatory input parameters:

    * nodes:     matrix of node positions (in mm), of size 3 x Nnodes (rows: xpos, ypos, zpos in mm);
                 a node is defined as the input (or as the output) of a pipe
    * radii:     matrix of pipe radii (in mm), of size Nnodes x Nnodes (element (i,j) stores the radius
                 of the pipe connecting node i with node j)
    * qin:       input volumetric flow rate in mm3/s
    * idxin:     index of the input node in the node matrix (nodes[:,idxin] provides the x,y,z coordinates
                 of such a input node)
    * idxout:    index of the output node in the node matrix (nodes[:,idxout] provides the x,y,z coordinates
                 of such an input node)

Additional optional parameters are the fluid viscosity, the radius at the inlet, the model used to solve the fluid dynamics and the type of solver. The full manual of the _pipenet_ class constructor is provided here.

For example, this code creates an initialises a _pipenet_ object describing a simple 2D 3-capillary network, depicted below. 

```
import numpy as np
import sys
sys.path.insert(0, '../code' )      # Add the SpinFlowSim "code" folder where the syn.py and visu.py files are stored

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
For the network initialisiation, we ensure that the fluid dynamics follows the model used in Blinder et al (Blinder P et al, Nature Neuroscience volume 16, pages 889–897 (2013), doi: [10.1038/nn.3426](https://doi.org/10.1038/nn.3426); optional parameter `flowmodel=bloodmodel`) with a cell-free plasma viscosity of 1.2 mPa s (`visc=muval`). The circuit is solved numerically (`solver=solvtechnique`), with default inlet radius computation (mean of the radii emanating from the inlet): 

Upon declaration of the _pipenet_ object, the vascular network is solved. We can now check, for example, the volumetric flow rate (VFR) in mm$^3$/s in each capillary:
```
>>> print(net.flowmat)
[[        nan  0.00138352  0.00411648]
 [-0.00138352         nan  0.00138352]
 [-0.00411648 -0.00138352         nan]]
```
or the blood velocity (BV) in mm/s in each capillary:
```
>>> print(net.velmat)
[[         nan  17.61550894  26.74115613]
 [-17.61550894          nan   6.88105818]
 [-26.74115613  -6.88105818          nan]]
```

Note that the VFR and BV are two attributes of the object `net`, namely `flowmat` and `velmat`. Elements `flowmat[i][j]` and `velmat[i][j]` respectively store the VFR and BV of the capillary going from node _i_ to node _j_. These two have the same sign, and is > 0 if the flow goes from node _i_ to node _j_, while it is < 0 if the flow goes from node _j_ to node _i_. By definition then, it follows that `flowmat[i][j]` and `flowmat[j][i]` (as well as `velmat[i][j]` and `velmat[j][i]`) have opposite signs.  


## Drawing a network on histology and using it to initialise a _pipenet_ object
With SpinFlowSim, one can simulate flow in realistic capillary networks that have been carefully reconstructed from histology. Here we show an example of a Hematoxylin and Eosin (HE) image of a mouse kidney slice, where a vascular network has been segmented by tracing visible capillaries.

The network is made of 26 nodes, connected among each other through 37 straight capillaries. The easiest way to start building this object is by loading a spreadsheet in .csv format with the network data as described below.
![labels_githhubreoi](https://github.com/user-attachments/assets/0364164f-4f12-4cf2-9fae-3c7f37770e81)

We have stored the information realted to the network as a CSV spreadsheet. Each straight capillary is described by a row in [Network_stats.csv](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/Network_stats.csv). Each row is described by the following variables:
- segmented ID (SegmentNumber)
- mean radius (RadiiMean_um) in µm; computed by averaging three radii measured over the length of a segment
- starting coordinates (Xstart,Ystart,Zstart) in µm
- ending coordinates (Xend,Yend,Zend) in µm
- starting node ID (Node_at_Xstart_Ystart)
- ending node ID (Node_at_Xend_Yend)

The table below illustrates the content of the CSV file.


| SegmentNumber | RadiiMean_um      | Xstart      | Ystart     | Zstart | Xend        | Yend       | Zend | Node_at_Xstart_Ystart | Node_at_Xend_Yend |
|---------------|-------------------|-------------|------------|--------|-------------|------------|------|-----------------------|-------------------|
| 1             | 3.3177            | 10056.36    | 10111.85   | 0      | 10083.96    | 10073.92   | 0    | 0                     | 1                 |
| 2             | 2.4370833333333333| 10083.96    | 10073.92   | 0      | 10182.86    | 10164.91   | 0    | 1                     | 2                 |
| 3             | 2.2855166666666666| 10182.86    | 10164.91   | 0      | 10088.94    | 10136.81   | 0    | 2                     | 3                 |
| 4             | 3.306783333333333 | 10056.36    | 10111.85   | 0      | 10088.94    | 10136.81   | 0    | 0                     | 3                 |
| 5             | 1.7060666666666666| 10088.94    | 10136.81   | 0      | 10062.82    | 10167.82   | 0    | 3                     | 4                 |
| 6             | 2.5319333333333334| 10062.82    | 10167.82   | 0      | 10093.35    | 10203.21   | 0    | 4                     | 5                 |
| 7             | 3.3291            | 10093.35    | 10203.21   | 0      | 10141.79    | 10180.75   | 0    | 5                     | 6                 |
| 8             | 2.8676999999999997| 10182.86    | 10164.91   | 0      | 10141.79    | 10180.75   | 0    | 2                     | 6                 |
| 9             | 1.6541833333333333| 10062.82    | 10167.82   | 0      | 10141.79    | 10180.75   | 0    | 4                     | 6                 |
| 10            | 2.0759166666666666| 10141.79    | 10180.75   | 0      | 10183.72    | 10200.81   | 0    | 6                     | 7                 |
| 11            | 2.3742666666666667| 10183.72    | 10200.81   | 0      | 10166.24    | 10247.25   | 0    | 7                     | 8                 |
| 12            | 3.064966666666667 | 10093.35    | 10203.21   | 0      | 10166.24    | 10247.25   | 0    | 5                     | 8                 |
| 13            | 2.6477            | 10166.24    | 10247.25   | 0      | 10185.83    | 10291.86   | 0    | 8                     | 9                 |
| 14            | 1.9695            | 10185.83    | 10291.86   | 0      | 10224.19    | 10278.96   | 0    | 9                     | 10                |
| 15            | 7.66195           | 10224.19    | 10278.96   | 0      | 10215.91    | 10251.86   | 0    | 10                    | 11                |
| 16            | 7.133399999999999 | 10215.91    | 10251.86   | 0      | 10199.99    | 10211.34   | 0    | 11                    | 12                |
| 17            | 1.7906166666666667| 10183.72    | 10200.81   | 0      | 10199.99    | 10211.34   | 0    | 7                     | 12                |
| 18            | 1.3943333333333332| 10199.99    | 10211.34   | 0      | 10270.11    | 10199.48   | 0    | 12                    | 13                |
| 19            | 3.5345333333333335| 10270.11    | 10199.48   | 0      | 10287.6     | 10217.21   | 0    | 13                    | 14                |
| 20            | 6.929283333333333 | 10215.91    | 10251.86   | 0      | 10287.6     | 10217.21   | 0    | 11                    | 14                |
| 21            | 4.208966666666666 | 10224.19    | 10278.96   | 0      | 10278.48    | 10278.52   | 0    | 10                    | 15                |
| 22            | 5.100483333333334 | 10278.48    | 10278.52   | 0      | 10287.6     | 10217.21   | 0    | 15                    | 14                |
| 23            | 2.8586833333333335| 10287.6     | 10217.21   | 0      | 10310.79    | 10215.24   | 0    | 14                    | 16                |
| 24            | 4.654883333333333 | 10265.38    | 10165.38   | 0      | 10287.6     | 10217.21   | 0    | 17                    | 14                |
| 25            | 3.3654833333333333| 10182.86    | 10164.91   | 0      | 10210.56    | 10180.14   | 0    | 2                     | 18                |
| 26            | 4.27815           | 10210.56    | 10180.14   | 0      | 10265.38    | 10165.38   | 0    | 18                    | 17                |
| 27            | 2.386616666666667 | 10265.38    | 10165.38   | 0      | 10271.09    | 10139.0    | 0    | 17                    | 19                |
| 28            | 2.9713666666666665| 10271.09    | 10139.0    | 0      | 10249.72    | 10119.15   | 0    | 19                    | 20                |
| 29            | 2.4316            | 10249.72    | 10119.15   | 0      | 10234.06    | 10135.73   | 0    | 20                    | 21                |
| 30            | 3.2838333333333334| 10234.06    | 10135.73   | 0      | 10210.74    | 10131.01   | 0    | 21                    | 22                |
| 31            | 3.657083333333333 | 10182.86    | 10164.91   | 0      | 10210.74    | 10131.01   | 0    | 2                     | 22                |
| 32            | 2.4389            | 10202.13    | 10106.51   | 0      | 10210.74    | 10131.01   | 0    | 23                    | 22                |
| 33            | 1.8349166666666665| 10202.13    | 10106.51   | 0      | 10172.3     | 10101.79   | 0    | 23                    | 24                |
| 34            | 1.6256833333333331| 10172.3     | 10101.79   | 0      | 10163.87    | 10137.31   | 0    | 24                    | 25                |
| 35            | 2.3390666666666666| 10163.87    | 10137.31   | 0      | 10182.86    | 10164.91   | 0    | 25                    | 2                 |
| 36            | 2.1939666666666664| 10172.3     | 10101.79   | 0      | 10134.48    | 10068.24   | 0    | 24                    | 26                |
| 37            | 1.82385           | 10083.96    | 10073.92   | 0      | 10134.48    | 10068.24   | 0    | 1                     | 26                |


Place the .csv file in the same directory as [script01_initnet.py](./script01_initnet.py), and run the latter to initalize the network.  


## Loading the realistic vascular networks that we distribute with SpinFlowSim

TEXT TO GO HERE
