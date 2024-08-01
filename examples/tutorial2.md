**WE ARE POLISHING THIS CODE. IT WILL BE HERE IN THE COMING DAYS - APOLOGIES FOR THE INCONVENIENCE!!**

This tutorial will show on our vascular networks data how synthetic signals can be used to inform microvasculature parameter estimation.


## Data for Plots
### Parameters

The parameter data used for generating the plots can be found in the `data_for_plots` folder for each network, grouped by network. 
e.g. Net1 data_for_plots `data_for_plots/Net1/` contains

- [Net1_Pars_vm.npy](data_for_plots/Net1/Net1_Pars_vm.npy)
- [Net1_Pars_vs.npy](data_for_plots/Net1/Net1_Pars_vs.npy)
- [Net1_Pars_anb.npy](data_for_plots/Net1/Net1_Pars_anb.npy)

Each of these .npy arrays has 100 entries, reporting average network properties for each network realization:
- **variable _vm_**: mean value of the blood velocity distribution, in `mm/s`;
- **variable _vs_**: standard deviation of the blood velocity distribution, in `mm/s`;
- **variable _anb_**: apparent network branching (ANB), expressed in `number of segments`.

### Signals

The signals used for the estimation can be found in the `data_for_plots` folder for each network, grouped by network. 
e.g. Net1 data_for_plots `data_for_plots/Net1/` contains

- [Net1_Sigs_TRSE.npy](data_for_plots/Net1/Net1_Sigs_TRSE.npy)
- [Net1_Sigs_subsetPGSE.npy](data_for_plots/Net1/Net1_Sigs_subsetPGSE.npy)
- [Net1_Sigs_richPGSE.npy](data_for_plots/Net1/Net1_Sigs_richPGSE.npy)

- Each of these arrays has 100 entries, storing signals for each simulated protocol:
- **TRSE**: twice-refocused spin echo (TRSE) pulse sequence. The protocol included b-values b = {0, 50, 100} s/mm2 acquired at 3 different diffusion times (δ1 = 8.9 ms, δ2 = 17.6 ms, δ3 = 20.4 ms, δ4 = 6.0 ms, ∆1,2 = 17.4 ms and ∆1,4 = 63.9 ms at short diffusion time; δ1 = 13.2 ms, δ2 = 19.3 ms, δ3 = 24.8 ms, δ4 = 7.7 ms, ∆1,2 = 21.7 ms and ∆1,4 = 74.2 ms at intermediate diffusion time; δ1 = 18.9 ms, δ2 = 21.0 ms, δ3 = 30.5 ms, δ4 = 9.5 ms, ∆1,2 = 27.5 ms and ∆1,4 = 87.5 ms at long diffusion time);
- 
- **subsetPGSE**: pulsed gradient spin echo (PGSE);

- 
- **richPGSE**: ultra-rich, theoretical pulsed gradient spin echo (PGSE).
