import numpy as np
import pickle as pk  
import sys
sys.path.insert(0, '../code' )      # Add the SpinFlowSim folder where the syn.py and visu.py files are stored
import syn
import visu
import os
from parsana import calcMicroPar, anaNetBran


Nmol_sim = 5000 # number of spins in simulation

# List of vascular networks
networkID = [f"Net{i}" for i in range(1, 16)]
Nnets = len(networkID)

# Loop over all folders containing vascular network realizations
for ss in range(Nnets):

    # Intialize lists storing microvascular parameters
    Pars_vm = []
    Pars_lm = []
    Pars_rm = []
    Pars_vs = []
    Pars_anb = []

    folder_path = os.path.join('..', 'networks', networkID[ss]) # path to folder containg vascular network realizations 
    files = os.listdir(folder_path)

    # Loop over all files containing vascular network realizations (files stored with `.bin` extension)
    for file in files:

        if file.endswith(".bin"): 
            infile_path = os.path.join(folder_path, file)  # Complete path to the infile
            with open(infile_path, 'rb') as net_file:  # Open the `.bin` file to load each network realization
                net = pk.load(net_file)

                vm, vs, lm, rm = calcMicroPar(net)  # Calculate the average microvascular parameters per network: 
                                                    # vm: mean velocity, vs: standard deviation of velocity, 
                                                    # lm: mean segment length, rm: mean radius
                                                    
                # Synthesize a diffusion-weighted MRI measurement from protons flowing within a vascular network for a unique b-value
                gdir = np.array((1.0,0.0,0.0))      # use fixed gradient direction
                out = net.dMRISynMea(
                    bval = 100.0, 
                    Gdur = 0.030, 
                    Gsep = 0.070, 
                    Gdir = gdir, 
                    Nspins = Nmol_sim, 
                    dt = 1e-5,  
                    )
                
                rmat = out[3] # Access matrix storing the trajectories of the spins (units: mm)
                
                anb = anaNetBran(rmat, Nmol_sim) # Obtain the average number of changes a spin experiences for the simulated DW-MRI measurement in the network
            
                Pars_vm.append(vm)
                Pars_vs.append(vs)
                Pars_lm.append(lm)
                Pars_rm.append(rm)
                Pars_anb.append(anb) 


    # Save properties in numpy arrays in folder containg vascular network
    saveas= os.path.basename(folder_path)
    target_folder = os.path.join('data_for_plots',networkID[ss])
    os.makedirs(target_folder, exist_ok=True) # create directory if it does not already exist
    
    np.save(os.path.join(target_folder, '{}_Pars_vm.npy'.format(saveas)), Pars_vm)
    np.save(os.path.join(target_folder, '{}_Pars_vs.npy'.format(saveas)), Pars_vs)
    np.save(os.path.join(target_folder, '{}_Pars_lm.npy'.format(saveas)), Pars_lm)
    np.save(os.path.join(target_folder, '{}_Pars_rm.npy'.format(saveas)), Pars_rm)
    np.save(os.path.join(target_folder, '{}_Pars_anb.npy'.format(saveas)), Pars_anb)


    