
# Tutorial 2: using synthetic dMRI signal to inform microvasculature parameter estimation

This tutorial shows you how the diffusion MRI (dMRI) signals you synthesized in Tutorial 1 can be used to estimate microvascular properties of vascular networks.  Please remember that to go through our tutorials you need python 3 with:
- [numpy](https://numpy.org) (developed with version 1.24.2)
- [pandas](https://pandas.pydata.org) (developed with version 1.5.3)
- [Lcapy](https://lcapy.readthedocs.io) (developed with version 1.10)
- [PySpice](https://github.com/FabriceSalvaire/PySpice) (developed with version 1.5)
- [graph-tool](https://graph-tool.skewed.de) (developed with version 2.45, commit b1a649d8)
- [nibabel](https://nipy.org/nibabel) (developed with version '5.1.0')
- [matplotlib](https://matplotlib.org) (examples tested with version '3.7.1')
- [mrtrix](https://www.mrtrix.org) (developed with version '3.0.4';  note that we use MrTrix command line tool `dwidenoise` from within python with `os.system()`, to estimate the noise level in synthetic vascular dMRI signals).

Tutorial 2 contains the following sections:
* [Data for Plots](#data-4-plots)
* [Parameter estimation from signals](#parameter-est)
* [Scatter plots showing correlation](#corr)
* [Final remarks](#remarks)

## Data for Plots  <a name="data-4-plots"></a>
### Parameters

The parameter data used for generating the plots can be found in the [`data_for_plots`](data_for_plots) folder for each network, grouped by network. 
e.g. Net1 data_for_plots `data_for_plots/Net1/` contains

- [Net1_Pars_vm.npy](data_for_plots/Net1/Net1_Pars_vm.npy)
- [Net1_Pars_vs.npy](data_for_plots/Net1/Net1_Pars_vs.npy)
- [Net1_Pars_anb.npy](data_for_plots/Net1/Net1_Pars_anb.npy)

Each of these .npy arrays has 100 entries, reporting average network properties for each network realization:
- **variable _vm_**: mean value of the blood velocity distribution, in `mm/s`;
- **variable _vs_**: standard deviation of the blood velocity distribution, in `mm/s`;
- **variable _anb_**: apparent network branching (ANB), expressed in `number of segments`.

These parameters are calculated and stored by running the [`calculate_microvascular_properties.py`](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks/calculate_microvascular_properties.py) script. 
### Signals

The signals used for the estimation can be found in the [`data_for_plots`](data_for_plots) folder for each network, grouped by network. 
e.g. Net1 data_for_plots `data_for_plots/Net1/` contains

- [Net1_Sigs_TRSE.npy](data_for_plots/Net1/Net1_Sigs_TRSE.npy)
- [Net1_Sigs_richPGSE.npy](data_for_plots/Net1/Net1_Sigs_richPGSE.npy)
- [Net1_Sigs_subsetPGSE.npy](data_for_plots/Net1/Net1_Sigs_subsetPGSE.npy)

Each of these arrays has 100 x n measurements, storing signals for each simulated protocol:

- **TRSE**: twice-refocused spin echo (TRSE) pulse sequence. The protocol included b-values b = {0, 50, 100} s/mm<sup>2</sup> acquired at 3 different diffusion times (δ1 = 8.9 ms, δ2 = 17.6 ms, δ3 = 20.4 ms, δ4 = 6.0 ms, ∆1,2 = 17.4 ms and ∆1,4 = 63.9 ms at short diffusion time; δ1 = 13.2 ms, δ2 = 19.3 ms, δ3 = 24.8 ms, δ4 = 7.7 ms, ∆1,2 = 21.7 ms and ∆1,4 = 74.2 ms at intermediate diffusion time; δ1 = 18.9 ms, δ2 = 21.0 ms, δ3 = 30.5 ms, δ4 = 9.5 ms, ∆1,2 = 27.5 ms and ∆1,4 = 87.5 ms at long diffusion time).
  
- **richPGSE**: ultra-rich pulsed gradient spin echo (richPGSE) with a total of 99 measurements. Consisting of 9 b = 0 and 10 non-zero b-values b = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100} s/mm<sup>2</sup>,  each acquired for 9 unique diffusion times, corresponding to (δ,∆) = {10, 20, 30} ms × {30, 50, 70} ms.
  
- **subsetPGSE**: a second PGSE protocol. It is a subset of the former, containing 3 b=0 and 6 diffusion-weighted (DW) measurements, namely b = {50, 100} for 3 different diffusion times. The gradient duration δ was fixed to 20 ms, while the 3 diffusion times were achieved by varying ∆ as ∆= {30, 50, 70} ms.

## Parameter estimation from signals <a name="par-est"></a>

The [run_estimation.py](data_for_plots/run_estimation.py) script carries out model fitting to estimate microvascular parameters from synthetic signals. 

The fitting is performed with [mri2micro_dictml.py](data_for_plots/mri2micro_dictml.py) tool, part of bodymritools [mri2micro_dictml.py](https://github.com/fragrussu/bodymritools/blob/main/mrifittools/mri2micro_dictml.py)
Note that mri2micro_dictml.py can be used to fit any equation-free, numerical signal model, given examples of signals and corresponding vascular parameters for any given acquisition protocol.

Begin by cloning this repository, navigate to [data_for_plots](data_for_plots) and execute `run_estimation.py` like in the following command:

```bash
python run_estimation.py --protocol subsetPGSE --snr 20
```
## Scatter plots showing correlation<a name="corr"></a>
<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/nets.png" width="950" height="auto">
</div>

