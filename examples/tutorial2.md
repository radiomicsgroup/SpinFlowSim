
# Tutorial 2: using synthetic dMRI signals to inform microvasculature parameter estimation

This tutorial shows you how the diffusion MRI (dMRI) signals you synthesized as explained in [Tutorial 1](https://github.com/radiomicsgroup/SpinFlowSim/edit/main/examples/tutorial1.md), can be used to inform the estimation of properties of the vascular networks.  Please remember that to go through our tutorials you need python 3 with:
- [numpy](https://numpy.org) (developed with version 1.24.2)
- [pandas](https://pandas.pydata.org) (developed with version 1.5.3)
- [Lcapy](https://lcapy.readthedocs.io) (developed with version 1.10)
- [PySpice](https://github.com/FabriceSalvaire/PySpice) (developed with version 1.5)
- [graph-tool](https://graph-tool.skewed.de) (developed with version 2.45, commit b1a649d8)
- [nibabel](https://nipy.org/nibabel) (developed with version '5.1.0')
- [scipy](https://scipy.org) (developed with version '1.11.3')
- [matplotlib](https://matplotlib.org) (examples tested with version '3.7.1')
- [mrtrix](https://www.mrtrix.org) (developed with version '3.0.4';  note that we use MrTrix command line tool `dwidenoise` from within python with `os.system()`, to estimate the noise level in synthetic vascular dMRI signals).

Tutorial 2 contains the following sections:
* [Description of the data of this tutorial](#data-4-plots)
* [Vascular parameter estimation informed by synthetic signals](#par-est)
* [Scatter plots showing parameter estimation results](#corr)

## Description of the data of this tutorial  <a name="data-4-plots"></a>
To run this tutorial, we need synthetic dMRI signals and corresponding sets of vascular parameters. We will work with the 15 vascular networks used in our [paper](https://doi.org/10.1016/j.media.2025.103531), which we distribute with _SpinFlowSim_ [here](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks).


### Vascular parameters

The vascular parameters required for this tutorial can be found in the [`data_for_plots`](data_for_plots) folder for each network. They are by network: e.g., folder [Net1](data_for_plots/Net1) contains parameters for network 1, [Net2](data_for_plots/Net2) for network 2, etc. Each of this folder contains three arrays: for network 1, for example, these are

- [Net1_Pars_vm.npy](data_for_plots/Net1/Net1_Pars_vm.npy)
- [Net1_Pars_vs.npy](data_for_plots/Net1/Net1_Pars_vs.npy)
- [Net1_Pars_anb.npy](data_for_plots/Net1/Net1_Pars_anb.npy)

Each of these .npy arrays has 100 entries, reporting average network properties for each network realization (remember that we obtain 100 different instantiations of each of our 15 networks, by varying inlet/outlet as well as the input flow - see more about this [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/networks/README.md)). The three array report values of three different vascular properties, namely:
- **variable _vm_**: mean value of the blood velocity distribution, in `mm/s`;
- **variable _vs_**: standard deviation of the blood velocity distribution, in `mm/s`;
- **variable _anb_**: apparent network branching (ANB), expressed in `number of segments`.

These parameters are calculated and stored by running the [`calculate_microPars.py`](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples/calculate_microPars.py) script. The [`parsana.py`](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples/parsana.py)
module contains the functions used for this analysis. You can read more about them in our [paper](https://doi.org/10.1016/j.media.2025.103531), or in this [README file](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks/README.md). 

### Synthetic dMRI signals

Alongside vascular properties, we need synthetic dMRI signals. We have generated signals for different dMRI protocols, which we also include in the [`data_for_plots`](data_for_plots) folder. Signals are again grouped by network, e.g., folder [Net1](data_for_plots/Net1) contains signals for network 1, [Net2](data_for_plots/Net2) for network 2, etc. Each of this folder contains three arrays: for network 1, for example, these are

- [Net1_Sigs_TRSE.npy](data_for_plots/Net1/Net1_Sigs_TRSE.npy)
- [Net1_Sigs_richPGSE.npy](data_for_plots/Net1/Net1_Sigs_richPGSE.npy)
- [Net1_Sigs_subsetPGSE.npy](data_for_plots/Net1/Net1_Sigs_subsetPGSE.npy)

Each of these arrays has 100 x n measurements, storing signals for three different simulated dMRI protocols:

- **TRSE**: a twice-refocused spin echo (TRSE) protocol.
  The protocol included b-values b = {0, 50, 100} s/mm<sup>2</sup> acquired at 3 different diffusion times:
  - δ1 = 8.9 ms, δ2 = 17.6 ms, δ3 = 20.4 ms, δ4 = 6.0 ms, ∆1,2 = 17.4 ms and ∆1,4 = 63.9 ms at short diffusion time;
  - δ1 = 13.2 ms, δ2 = 19.3 ms, δ3 = 24.8 ms, δ4 = 7.7 ms, ∆1,2 = 21.7 ms and ∆1,4 = 74.2 ms at intermediate diffusion time;
  - δ1 = 18.9 ms, δ2 = 21.0 ms, δ3 = 30.5 ms, δ4 = 9.5 ms, ∆1,2 = 27.5 ms and ∆1,4 = 87.5 ms at long diffusion time.
  
- **richPGSE**: an ultra-rich pulsed gradient spin echo (PGSE) protocol, referred to as _richPGSE_.
  This protocol consists of 9 b = 0 and 10 non-zero b-values b = {10, 20, 30, 40, 50, 60, 70, 80, 90, 100} s/mm<sup>2</sup>,  each acquired for 9 unique diffusion times:
  - (δ,∆) = {10, 20, 30} ms × {30, 50, 70} ms.
  
- **subsetPGSE**: a second PGSE protocol that is a subset of the former, referred to as _subsetPGSE_.
  It contains 3 b=0 and 6 diffusion-weighted (DW) measurements, namely b = {50, 100} for 3 different diffusion times:
  - δ was fixed to 20 ms, while the 3 diffusion times were achieved by varying ∆ as ∆= {30, 50, 70} ms.


Signals were synthesised using the same approach described in previous [tutorial 1](https://github.com/radiomicsgroup/SpinFlowSim/edit/main/examples/tutorial1.md). The figure below illustrates the TRSE and PGSE diffusion encodings used in for vascular signal synthesis. 

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/sequences.png" width="650" height="auto">
</div>
<div align="center">
  (a) PGSE pulse sequence (b) TRSE pulse sequence
</div>

## Vascular parameter estimation informed by synthetic signals <a name="par-est"></a>

Now that we have put together rich sets of synthetic dMRI signals and corresponding microvascular parameters from 15 vascular networks, we will illustrate how these can be used to devise a strategy enabling the microvascular parameter estimation from new, unseen signals. 

In short, we use noise-free signals and corresponding vascular parameters from 14 networks to learn a numerical forward signal model that maps vascular parameters to dMRI signals. Afterwards, we plut such a forward signal model into standard maximum-likelihood fitting, through which we estimate vascular parameters from noisy signals from 15th substrate. We perform this experiment in a leave-one-out fashin, obtaining vasuclar parameter estimates for all 15 networks. 

The [`run_estimation.py`](data_for_plots/run_estimation.py) script carries out this experiment. The script relies on the [`mri2micro_dictml.py`](data_for_plots/mri2micro_dictml.py) tool, which is essentially a copy of a script5 included as part of our **bodymritools** repository (available [here](https://github.com/fragrussu/bodymritools)). `mri2micro_dictml.py` can be used to fit any equation-free, numerical signal model, given examples of signals and corresponding vascular parameters for any given acquisition protocol.

To run [`run_estimation.py`](data_for_plots/run_estimation.py), just navigate to [data_for_plots](data_for_plots) folder in your command line, and execute the script, selecting the protocol and a signal-to-noise ratio (SNR) level with options `--protocol` and `--snr`. For example, to estimate vascular parameters for protocol _subsetPGSE_ and an SNR of 20, simply execute:

```
python run_estimation.py --protocol subsetPGSE --snr 20
```

Note that the [`mri2micro_dictml.py`](data_for_plots/mri2micro_dictml.py) tool requires:
- examples of noise-free signals and corresponding vascular parameters, to build a numerical forward signal model;
- the noisy signals on which vascular parameter estimation should be performed as a 4D NIFTI file;
- a text file with the diffusion protocol;
- and an optional noise map, also in NIFTI format, in case you want to account for the noise floor in maximum-likelihood fitting.  

The necessary files are included in [data_for_plots](data_for_plots), and used directly by [`run_estimation.py`](data_for_plots/run_estimation.py). An help manual of [`mri2micro_dictml.py`](data_for_plots/mri2micro_dictml.py) can be found [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/mri2micro_dictml_manual.md).

## Scatter plots showing parameter estimation results<a name="corr"></a>

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/estimation.PNG" width="950" height="auto">
</div>

Visualizing the estimation results as a scatter plot that includes the computation of a correlation coefficient between ground truth and estimated parameters, is as easy as running the script [`get_plots.py`](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/data_for_plots/get_plots.py), which we include in the from [data_for_plots](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/data_for_plots) folder. This can be run easily like this:

```
python get_plots.py --save True 
```
