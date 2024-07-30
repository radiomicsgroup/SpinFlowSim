import numpy as np
from matplotlib import pyplot as plt
import syn
import visu
import pandas as pd
import pickle as pk
import time

### Load spreadsheet and set user-defined network parameters
nodein = 0                 # Define an input node
nodeout = 1                # Define an output node
qinput = 0.0055            # Define an input volumetric flow rate in mm3/s
muval = 1.2                # Define the viscosity of the fluid in mPa x s
bloodmodel = 'blinder'     # Define the flow model
outfile = 'script01_initnet.bin'        # Name of output file
dset = pd.read_csv('Network_stats.csv') # Name of spreadsheet with vascular network data

### Get relevant fields from the spreadsheet
rads = dset.RadiiMean_um   # Mean radius in um


ninstr = dset.Node_at_Xstart_Ystart # Input node of segment
noutstr = dset.Node_at_Xend_Yend    # Output node of segment


Xend = dset.Xend
Yend = dset.Yend
Zend = dset.Zend


Xstart = dset.Xstart
Ystart = dset.Ystart
Zstart = dset.Zstart


ninlist = []
for item in ninstr:
    ninlist.append( int( item ) )  
ninlist = np.array(ninlist)


noutlist = []
for item in noutstr:
    noutlist.append( int( item ) ) 
noutlist = np.array(noutlist)


radslist = []
for item in rads:
    radslist.append( float( item ) )
radslist = np.array(radslist)


Xstartlist = []
for item in Xstart:
    Xstartlist.append( float( item ) )
Xstartlist = np.array(Xstartlist)


Ystartlist = []
for item in Ystart:
    Ystartlist.append( float( item ) )
Ystartlist = np.array(Ystartlist)


Zstartlist = []
for item in Zstart:
    Zstartlist.append( float( item ) )
Zstartlist = np.array(Zstartlist)


Xendlist = []
for item in Xend:
    Xendlist.append( float( item ) )
Xendlist = np.array(Xendlist)


Yendlist = []
for item in Yend:
    Yendlist.append( float( item ) )
Yendlist = np.array(Yendlist)


Zendlist = []
for item in Zend:
    Zendlist.append( float( item ) )
Zendlist = np.array(Zendlist)

### Find out how many nodes and segments there are
Nseg = dset.shape[0]      # Number of segments
Nnodes = np.max( np.array([np.max(ninlist),np.max(noutlist)]) )  + 1  # Number of nodes (maximum index + 1)


### Allocate matrices to store the network properties
nodepos = np.zeros((3,Nnodes))      # Coordinates of the nodes
already_checked = np.zeros(Nnodes,dtype='bool')

# Scan both ninlist and noutlist to find all nodes
for ii in range(0,Nseg):
    nodepos[0,ninlist[ii]] = Xstartlist[ii]
    nodepos[1,ninlist[ii]] = Ystartlist[ii]
    nodepos[2,ninlist[ii]] = Zstartlist[ii]
    already_checked[ninlist[ii]] = True

for ii in range(0,Nseg):
    if(not already_checked[noutlist[ii]]):
        nodepos[0,noutlist[ii]] = Xendlist[ii]
        nodepos[1,noutlist[ii]] = Yendlist[ii]
        nodepos[2,noutlist[ii]] = Zendlist[ii]


# Center in 0 and convert to mm
nodepos[0,:] = nodepos[0,:] - nodepos[0,nodein]
nodepos[1,:] = nodepos[1,:] - nodepos[1,nodein]
nodepos[2,:] = nodepos[2,:] - nodepos[2,nodein]
nodepos = nodepos/1000.0

# Create radius matrix in mm
radiusmat = np.zeros((Nnodes,Nnodes))
for ii in range(0,Nseg):
    radiusmat[ninlist[ii],noutlist[ii]] = 0.001*radslist[ii]
    radiusmat[noutlist[ii],ninlist[ii]] = 0.001*radslist[ii]

for nn in range(0,Nnodes):
    radiusmat[nn,nn] = np.nan


### Initialise network and save it
t1 = time.time()
net = syn.pipenet(nodepos,radiusmat,qinput,nodein,nodeout,visc=muval,flowmodel=bloodmodel)
t2 = time.time()
save_file = open(outfile,'wb')
pk.dump(net,save_file,pk.HIGHEST_PROTOCOL)      
save_file.close()

save_file = open(outfile,'wb')
pk.dump(net,save_file,pk.HIGHEST_PROTOCOL)      
save_file.close()

print('')
print('*** It took {} seconds to initialise the network'.format(t2 - t1))
print('')

