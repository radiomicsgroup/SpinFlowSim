import numpy as np
import nibabel as nib
import os
import argparse


np.random.seed(42)

###################################################################################################################################################
###################################################################################################################################################

parser = argparse.ArgumentParser(
    description='Run fitting to estimate microvascular properties.')
parser.add_argument('--protocol', type=str, required=True, help='Input protocol (e.g., subsetPGSE, TRSE, richPGSE, custom)')
parser.add_argument('--snr', type=int, required=True, help='Signal-to-noise ratio')
parser.add_argument('--bval_file', type=str, help='Path to the .bval file if protocol is not one of the predefined (subsetPGSE, TRSE, richPGSE)')
args = parser.parse_args()
protocol = args.protocol
snr = args.snr
bval_file = args.bval_file

print(f"SNR: {snr}")
print(f"Protocol: {protocol}")
if protocol not in ['subsetPGSE', 'TRSE', 'richPGSE']:
    if not bval_file or not os.path.isfile(bval_file):
        raise ValueError("Protocol bval data not found. You must supply a .bval file with the --bval_file argument.")
    print(f'Using custom protocol bvals defined in: {bval_file}')

###################################################################################################################################################
###################################################################################################################################################

# List of vascular networks
networkID = [f"Net{i}" for i in range(1, 16)]
Nnets = len(networkID)

# List of parameters
Pars_list = ["Pars_vm.npy","Pars_vs.npy","Pars_anb.npy"]
Npars = len(Pars_list)

# Load signals for each Net
Sig = [np.load(f"{networkID[i]}/{networkID[i]}_Sigs_{protocol}.npy") for i in range(Nnets)]
# Load Par_vm, Par_vs, and Par_anb for each Net
Par_vm = [np.load(f"{networkID[i]}/{networkID[i]}_Pars_vm.npy") for i in range(Nnets)]
Par_vs = [np.load(f"{networkID[i]}/{networkID[i]}_Pars_vs.npy") for i in range(Nnets)]
Par_anb  = [np.load(f"{networkID[i]}/{networkID[i]}_Pars_anb.npy") for i in range(Nnets)]

# Initialize the sig arrays
Sigsize = Sig[0].size
s_list = [np.zeros(Sigsize) for _ in range(Nnets)]
# Initialize the par arrays
Parsize = Par_vs[0].size  
par_list = [np.zeros((Parsize, Npars)) for _ in range(Nnets)]

# Function to populate parameter 'dictionaries'
def create_Pars_dictionary(vm, vs, anb):
    par = np.zeros((vs.size, Npars))
    par[:, 0] = vm
    par[:, 1] = vs
    par[:, 2] = anb
    return par

# Populate par_list and sig_list with data from each network
for nn in range(15):
    par_list[nn] = create_Pars_dictionary(Par_vm[nn], Par_vs[nn], Par_anb[nn])
    s_list[nn] =Sig[nn]

# Assign ID to each sig and par
s1, s2, s3, s4, s5, par6, s7, s8, s9, s10, s11, s12, s13, s14, s15 = s_list
par1, par2, par3, par4, par5, par6, par7, par8, par9, par10, par11, par12, par13, par14, par15 = par_list



for i in range(len(par_list)):

    dir_name = f'testfit0{i+1}' # Define the test directory name 
    if not os.path.exists(dir_name): # Create directory if non-existing
        os.makedirs(dir_name)
        
    #Create dictionaries in a leave-one-out fashion, storing parameters and signals from N-1 Networks to learn the forward model
    par_train = np.concatenate(par_list[:i] + par_list[i+1:], axis=0)
    sig_train = np.concatenate(s_list[:i] + s_list[i+1:], axis=0)

    np.save(f'par{i+1}_trainset.npy', par_train)
    np.save(f'sig{i+1}_trainset.npy', sig_train)

    print(f'par_trainset shape is: {par_train.shape}')

    ### Create noisy test set
    # snr = 20    # Define SNR
    s_noisy = np.sqrt(   ( s_list[i] + (1.0/snr)*np.random.randn(s_list[i].shape[0],s_list[i].shape[1])  )**2   +   (  (1.0/snr)*np.random.randn(s_list[i].shape[0],s_list[i].shape[1])  )**2  )
    nmeas = s_noisy.shape[1]   # Number of signal measurements 

    ### Convert signals and parameters as a 10x10 matrix of voxels
    sv_mat = np.zeros((10,10,1,nmeas)) # Initalize pure vascular signal
    s_noisy_mat = np.zeros((10,10,1,nmeas))   # Initalize pure vascular signal with Rician noise
    par_mat = np.zeros((10,10,1,Npars)) # Initalize vascular parameters matrix
    
    vv = 0
    for rr in range(0,10):
        for cc in range(0,10):
            sv_mat[rr,cc,0,:] = s_list[i][vv,:] # Populate pure vascular signal matrix
            s_noisy_mat[rr,cc,0,:] = s_noisy[vv,:] # Populate noisy vascular signal matrix
            par_mat[rr,cc,0,:] = par_list[i][vv,:] # Populate parameter matrix
            vv = vv + 1


    ### Create header affine to store signals and properties as NIFTIs
    Amat = np.eye(4)
    Amat[0,0] = 10/1000.0
    Amat[1,1] = 10/1000.0
    Amat[2,2] = 1/1000.0


    img = nib.Nifti1Image(sv_mat, Amat)    
    img.header.set_xyzt_units(2)
    img.set_data_dtype('float64')
    nib.save(img, f'sv{i+1}_mat.nii')      # Save pure vascular signals as NIFTI


    img = nib.Nifti1Image(s_noisy_mat, Amat)
    img.header.set_xyzt_units(2)
    img.set_data_dtype('float64')
    nib.save(img, f's{i+1}_noisy_mat_SNR{snr}.nii')   # Save noisy signals as NIFTI

  
    img = nib.Nifti1Image(par_mat, Amat)
    img.header.set_xyzt_units(2)
    img.set_data_dtype('float64')
    nib.save(img, f'par{i+1}_mat.nii')       # Save vascular parameters as NIFTI


    # Create .bval file with list of b-values, based on input protocol
    if protocol == 'subsetPGSE':
        print(f'The {protocol} bval list is: ')
        print('0.0 50.0 100.0 0.0 50.0 100.0 0.0 50.0 100.0')
        os.system('echo "0.0 50.0 100.0 0.0 50.0 100.0 0.0 50.0 100.0" > s{}_noisy_mat.bval'.format(i+1))
    elif protocol == "RichPGSE":
        print(f'The {protocol} bval list is: ')
        print('[0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0] x 9')
        os.system('echo "0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 0.0 10.0 20.0 30.0 40.0 50.0 60.0 70.0 80.0 90.0 100.0 " > s{}_noisy_mat.bval'.format(i+1))
    elif protocol == 'TRSE':
        print(f'The {protocol} bval list is: ')
        print('0.0 50.0 100.0 0.0 50.0 100.0 0.0 50.0 100.0')
        os.system('echo "0.0 50.0 100.0 0.0 50.0 100.0 0.0 50.0 100.0" > s{}_noisy_mat.bval'.format(i+1))
    else:
        print(f'Using custom protocol with bval: {bval_file}')
        with open(bval_file, 'r') as f:
            bval_content = f.read().strip()
            os.system(f'echo "{bval_content}" > s{i+1}_noisy_mat.bval')


    #### Estimate noise level
    os.system('dwidenoise -extent 3,3,1 -noise s{}_noisy_mat_SNR{}_sigma.nii s{}_noisy_mat_SNR{}.nii --force buffer.nii'.format(i+1,snr,i+1,snr))
    os.system('rm -v buffer.nii')

    #### Run fitting
    os.system('python dictmaxlik.py s{}_noisy_mat_SNR{}.nii sig{}_trainset.npy par{}_trainset.npy testfit0{}/fit --noise s{}_noisy_mat_SNR{}_sigma.nii --ncpu 24 --sldim 0'.format(i+1,snr,i+1,i+1,i+1,i+1,snr))
