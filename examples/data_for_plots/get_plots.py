from scipy import stats
import nibabel as nib
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import argparse


# Set up argument parser
parser = argparse.ArgumentParser(
    description='Generate plots showing ground truth versus prediction microvascular properties.')
parser.add_argument('--save', type=bool, required=True, help='Option to save plots as .png')

args = parser.parse_args()
save = args.save

print(f"Option to save: {save}")



### Function to calculate Spearman and Pearson correlation between prediction and ground truth results
def get_corr_metrics(prediction, groundtruth):
    '''
    Function to calculate Spearman rank-order correlation coefficient 
    and Pearson correlation coefficient, both with associated p-values,
    between a list of prediction and a list of ground truth results. 

    '''

    prediction = np.asarray(prediction).flatten()
    groundtruth = np.asarray(groundtruth).flatten()
    
    spearman_corr = stats.spearmanr(prediction, groundtruth)
    pearson_corr = stats.pearsonr(prediction, groundtruth)
    
    return spearman_corr, pearson_corr



### Function to plot the prediction results, optionally displaying the Spearman and Pearson correlation values
def plot_parameter(prediction, ground_truth, label=None, units=None, spearman_corr=None, pearson_corr=None, save=None):

    '''
    Create plot to visualize the goodness of parameter estimation by means of a scatter plot,
    for Ground truth (x-axis) and Predicted (y-axis) results. 
    Optionally it displays Spearman and Pearson correlation values in text.

    MANDATORY INPUT PARAMETERS
        * prediction: list containing predicted/estimated parameters for all networks (units: parameter units);
        * ground_truth: aist containing ground truth parameters for all networks (units: parameter units);

    OPTIONAL INPUT PARAMETERS
        * spearman_corr: array containing Spearman's rank correlation coefficient `r_s`. If provided, the correlation value will be displayed on the plot.
        * pearson_corr: array containing Pearson's correlation coefficient `r`. If provided, the correlation value will be displayed on the plot.
        * save: (bool) If `True`, the plot will be saved as a `.png` file. If `False` or not provided, the plot will only be displayed on screen.


    '''

    
    # Flatten list arrays and create DataFrame for plotting
    data = np.array([ground_truth, prediction])
    plot_df = pd.DataFrame({'X': data[0].flatten(), 'Y': data[1].flatten()})

    # Create plot
    fig, ax = plt.subplots(figsize=(14, 12))
    ax.scatter(plot_df['X'], plot_df['Y'], s=100, marker="+", linewidths=3)

    # Create and set labels for the axes
    xlabel = f"Ground truth {label} {units}" 
    ylabel = f"Predicted {label} {units}" 
    ax.set_xlabel(xlabel)
    ax.set_ylabel(ylabel)

    # Optional: Display correlation values on the plot
    if spearman_corr is not None and pearson_corr is not None:
        corr_text = f"Spearman's $r_s$:{spearman_corr[0]:.2f}\nPearson's $r$: {pearson_corr[0]:.2f}"
        plt.gca().text(0.05, 0.95, corr_text, transform=ax.transAxes, fontsize=12, verticalalignment='top')


    plt.title(f"{label} estimation")
    plt.tight_layout()
    
    # Save or show the plot
    if save:
        clean_label = label.replace('$', '')
        plt.savefig(f"{clean_label}_estimation_plot.png")
    else:
        plt.show()




# List of vascular networks
networkID = [f"Net{i}" for i in range(1, 16)]
Nnets = len(networkID)

# Define array lists storing parameters for all networks
Par_1_gt= []   # Ground truth for parameter 1
Par_1_pre = [] # Predicted/Estimated values for parameter 1

Par_2_gt = []   # Ground truth for parameter 2
Par_2_pre = [] # Predicted/Estimated values for parameter 2

Par_3_gt = []   # Ground truth for parameter 3
Par_3_pre = [] # Predicted/Estimated values for parameter 3

# Load ground truth and prediction parameters for each network, and store in list
for ss in range(1, Nnets + 1):
    ground_truth_path= f'par{ss}_mat.nii'    # Path to ground truth results 
    img= nib.load(ground_truth_path)          # Load NiBabel image with ground truth results 
    gt_data = img.get_fdata()                 # Access image data as NumPy array
    del img

    gt_1= gt_data[:,:,0,0]                    # Ground truth results for Parameter 1 
    flattened_gt_1 = np.ndarray.flatten(gt_1) # Flatten data
    Par_1_gt.append(flattened_gt_1)           # Add to list

    gt_2= gt_data[:,:,0,1]                    # Access ground truth for Parameter 2   
    flattened_gt_2 = np.ndarray.flatten(gt_2) # Flatten data
    Par_2_gt.append(flattened_gt_2)           # Add to list

    gt_3= gt_data[:,:,0,1]                    # Access ground truth for Parameter 3  
    flattened_gt_3 = np.ndarray.flatten(gt_3) # Flatten data
    Par_3_gt.append(flattened_gt_3)           # Add to list


    prediction_path_1= f'testfit0{ss}/fit_par1.nii'     # Path to prediction results for Parameter 1
    img = nib.load(prediction_path_1)                   # Load NiBabel image with prediction results 
    pre_1_data = img.get_fdata()                        # Access image data as NumPy array
    del img
    flattened_pre_1= np.ndarray.flatten(pre_1_data)
    Par_1_pre.append(flattened_pre_1)                   # Add flattened array to list


    prediction_path_2= f'testfit0{ss}/fit_par2.nii'     # Path to prediction results for Parameter 2
    img = nib.load(prediction_path_2)                   # Load NiBabel image with prediction results 
    pre_2_data = img.get_fdata()                        # Access image data as NumPy array
    del img
    flattened_pre_2= np.ndarray.flatten(pre_2_data)
    Par_2_pre.append(flattened_pre_2)                   # Add flattened array to list


    prediction_path_3= f'testfit0{ss}/fit_par3.nii'     # Path to prediction results for Parameter 3
    img = nib.load(prediction_path_3)                   # Load NiBabel image with prediction results 
    pre_3_data = img.get_fdata()                        # Access image data as NumPy array
    del img
    flattened_pre_3 = np.ndarray.flatten(pre_3_data)
    Par_3_pre.append(flattened_pre_3)                   # Add flattened array to list


# calculate Spearman and Pearson correlation between prediction and ground truth results
sp1, p1 = get_corr_metrics(Par_1_pre, Par_1_gt)  
sp2, p2 = get_corr_metrics( Par_2_pre, Par_2_gt,)  
sp3, p3 = get_corr_metrics( Par_3_pre, Par_3_gt,)  



# Define units and labels 
units = {
    'par1': '[mm/s]',
    'par2': '[mm/s]',
    'par3': '[segments]',
}

label = {
    'par1': r'$v_m$',
    'par2': r'$v_s$', 
    'par3': r'$ANB$',
}


correlations = {
    'par1': {'spearman': sp1, 'pearson': p1},
    'par2': {'spearman': sp2, 'pearson': p2},
    'par3': {'spearman': sp3, 'pearson': p3},
}


prediction = {
    'par1': Par_1_pre,
    'par2': Par_2_pre,
    'par3': Par_3_pre,
}

ground_truth = {
    'par1': Par_1_gt,
    'par2': Par_2_gt,
    'par3': Par_3_gt,
}


# Plot the prediction results for each parameter, optionally saving the plot
for param in ['par1', 'par2','par3']:

    plot_parameter(

        prediction = prediction[param],
        ground_truth = ground_truth[param],
        label = label[param],
        units = units[param],
        spearman_corr=correlations[param]['spearman'],
        pearson_corr=correlations[param]['pearson'],
        save = save,
         
    )

