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

    # Intialize lists storing microvascular parameters
    Pars_vm = []
    Pars_vs = []
    Pars_vw = []
    Pars_qm = []
    Pars_qs = []
    Pars_qw = []
    Pars_lmp = []
    Pars_lm = []
    Pars_npaths = []
    Pars_rm = []
    Pars_rw = []
    Pars_anb = []
    
    folder_path = os.path.join('..', 'networks', networkID[ss]) # path to folder containg vascular network realizations 
    files = os.listdir(folder_path)

    # Loop over all files containing vascular network realizations (files stored with `.bin` extension)
    for file in files:

        if file.endswith(".bin"): 
            infile_path = os.path.join(folder_path, file)  # Complete path to the infile
            with open(infile_path, 'rb') as net_file:  # Open the `.bin` file to load each network realization
                net = pk.load(net_file)

                # Calculate the average microvascular parameters per network
                vm, vs, vm_pw, qm, qs, qm_pw, lmp, lm, Npaths, rm, rm_pw, anb = net.calcMicroPar(Nspins=5000)    # Nspins = number of spins for anb calculation

                Pars_vm.append(vm)
                Pars_vs.append(vs)
                Pars_vw.append(vm_pw)
                Pars_qm.append(qm)
                Pars_qs.append(qs)
                Pars_qw.append(qm_pw)
                Pars_lmp.append(lmp)
                Pars_lm.append(lm)
                Pars_npaths.append(Npaths)
                Pars_rm.append(rm)
                Pars_rw.append(rm_pw)
                Pars_anb.append(anb) 


    # Save properties in numpy arrays in folder containg vascular network
    saveas= os.path.basename(folder_path)
    target_folder = os.path.join('data_for_plots',networkID[ss])
    os.makedirs(target_folder, exist_ok=True) # create directory if it does not already exist
    
    np.save(os.path.join(target_folder, '{}_Pars_vm.npy'.format(saveas)), Pars_vm)
    np.save(os.path.join(target_folder, '{}_Pars_vs.npy'.format(saveas)), Pars_vs)
    np.save(os.path.join(target_folder, '{}_Pars_vw.npy'.format(saveas)), Pars_vw)
    np.save(os.path.join(target_folder, '{}_Pars_qm.npy'.format(saveas)), Pars_qm)
    np.save(os.path.join(target_folder, '{}_Pars_qs.npy'.format(saveas)), Pars_qs)
    np.save(os.path.join(target_folder, '{}_Pars_qw.npy'.format(saveas)), Pars_qw)
    np.save(os.path.join(target_folder, '{}_Pars_lmp.npy'.format(saveas)), Pars_lmp)
    np.save(os.path.join(target_folder, '{}_Pars_lm.npy'.format(saveas)), Pars_lm)
    np.save(os.path.join(target_folder, '{}_Pars_npaths.npy'.format(saveas)), Pars_npaths)
    np.save(os.path.join(target_folder, '{}_Pars_rm.npy'.format(saveas)), Pars_rm)
    np.save(os.path.join(target_folder, '{}_Pars_rw.npy'.format(saveas)), Pars_rw)
    np.save(os.path.join(target_folder, '{}_Pars_anb.npy'.format(saveas)), Pars_anb)




    