import numpy as np
import pickle as pk  
import sys
sys.path.insert(0, '../code' )      # Add the SpinFlowSim folder where the syn.py and visu.py files are stored
import syn
import visu
import os


# List of vascular networks
networkID = [f"Net{i}" for i in range(1, 16)]
Nnets = len(networkID)

# Loop over all folders containing vascular network realizations
for ss in range(Nnets):

    # Intialize lists storing parameters
    Pars_vm = []
    Pars_lm = []
    Pars_rm = []
    Pars_vs = []

    folder_path = os.path.join('Networks', networkID[ss]) # path to folder containg vascular network realizations
    files = os.listdir(folder_path)

    # Loop over all files containing vascular network realizations
    for file in files:

        if file.endswith(".bin"):

            infile_path = os.path.join(folder_path, file)  # Complete path to the infile

            with open(infile_path, 'rb') as net_file:  # Open the '.bin' file to load each network realization

                net = pk.load(net_file)
                connectivity = np.abs(net.connmat)  # Connectivity to know which nodes are connected


                # Access velocities array
                v = np.abs(net.velmat)
                varray = v[connectivity==1]         # Exclude nodes with -1 connectivity due to symmetry
                vm = np.mean(varray)                # Calculate mean velocity  (mm/s)
                vs = np.std(varray)                 # Calculate standard deviation of velocity  (mm/s)

                # Access matrix of segment radii
                r = np.abs(net.radmat)
                marray = r[connectivity==1]         # Exclude nodes with -1 connectivity due to symmetry
                rm = np.mean(marray)                # Calculate mean radius for network (mm)

                # Access matrix of segement lengths 
                l = np.abs(net.lengthmat)
                larray = l[connectivity==1]         # Exclude nodes with -1 connectivity due to symmetry
                lm = np.mean(larray)                # Calculate mean segment length for network (mm)

                Pars_vm.append(vm)
                Pars_vs.append(vs)
                Pars_lm.append(lm)
                Pars_rm.append(rm)

    # Save properties in numpy arrays in folder containg vascular network
    saveas= os.path.basename(folder_path)
    np.save(os.path.join(folder_path, '{}_Pars_vm.npy'.format(saveas)), Pars_vm)
    np.save(os.path.join(folder_path, '{}_Pars_vs.npy'.format(saveas)), Pars_vs)
    np.save(os.path.join(folder_path, '{}_Pars_lm.npy'.format(saveas)), Pars_lm)
    np.save(os.path.join(folder_path, '{}_Pars_rm.npy'.format(saveas)), Pars_rm)


    