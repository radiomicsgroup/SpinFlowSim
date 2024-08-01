**WE ARE POLISHING THIS CODE. IT WILL BE HERE IN THE COMING DAYS - APOLOGIES FOR THE INCONVENIENCE!!**

This tutorial will show on our vascular networks data how synthetic signals can be used to inform microvasculature parameter estimation.


## Data for Plots

The data used for generating the plots can be found in the `data_for_plots` folder for each network, grouped by network. 
e.g. Net1 data_for_plots `data_for_plots/Net1/` contains

- [Net1_Pars_vm.npy](data_for_plots/Net1/Net1_Pars_vm.npy)
- [Net1_Pars_vs.npy](data_for_plots/Net1/Net1_Pars_vs.npy)
- [Net1_Pars_anb.npy](data_for_plots/Net1/Net1_Pars_anb.npy)

Each of these .npy arrays has 100 entries, reporting average network properties for each network realization:
- **variable _vm_**: mean value of the blood velocity distribution, in `mm/s`;
- **variable _vs_**: standard deviation of the blood velocity distribution, in `mm/s`;
- **variable _anb_**: apparent network branching (ANB), expressed in `number of segments`.
