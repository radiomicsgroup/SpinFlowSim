# Capillary networks for vascular dMRI signal synthesis

SpinFlowSim was showcased on vascular networks reconstructed from 2D histological images. The networks were drawn manually, tracing visible capillaries in 11 stained liver biopsies. We have drawn **15 vascular networks** on human liver biopsies, as detailed in our [paper](https://doi.org/10.1016/j.media.2025.103531):

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/nets.png" width="950" height="auto">
</div>

For each network, we varied the input volumetric flow rate (VFR) in [1.5 × 10⁻⁴ ; 5.5 × 10⁻³] mm³/s and inlet/outlet nodes (10 input VFRs × 10 unique pairs of inlet/outlet), for a total of 100 realisations per network. This lead to **a total of 1500 vascular networks, which we include here in our repository**. You can start running SpinFlowSim by choosing from any of the 1500 network realizations (in binary `.bin` format) stored within this folder, and synthesising corresponding diffusion MRI (dMRI) signals for any protocol of interest. In [this tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial1.md) we show how.

### Directory Description

The [**networks**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks) folder contains 15 sub-folders, each referring to a vascular network.

- **Net{ID}**: folder containing the files, where `{ID}` is the unique identifier for the networks shown in the figure above.

### Sub-Directory Filename Description
Each network folder, contains several files - something like, for example, [_net10_Nin0_Nout28_Qin0.0001.bin_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/networks/Net10/net10_Nin0_Nout28_Qin0.0001.bin). We use the following file naming convention:
- **net{ID}**: Base name of the file, where `{ID}` corresponds to the same unique identifier used in the directory.
- **Nin{input_node}**: The flow input node in the given network realization.
- **Nout{output_node}**: The flow output node in the given network realization.
- **Qin{input_volumetric_flow_rate}**: Represents the input volumetric flow rate in mm³/s.

### _pipenet_ network objects (.bin files)
Resolved vascular networks are stored as binary (.bin) files, containing instantiations of the class _pipenet_ (i.e., _pipenet_ objects; the _pipenet_ class is defined in [_syn.py_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/syn.py)), and is described in detail [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/pipenet_manual.md). 

A file name `net10_Nin0_Nout28_Qin0.0001.bin` indicates that:
- the file is located in the folder corresponding to the `Net10` vascular network;
- the network ID is `10`;
- the flow inlet node is `0`, while the the outlet node is `28`;
- the input volumetric flow rate is `1 × 10⁻⁴ mm³/s`.

### Network properties (.csv files)
Each _pipenet_ object storing a network is accompanied by a small CSV file bearing the same name, and storing summary network properties. For example, the network file  [_net10_Nin0_Nout28_Qin0.0001.bin_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/networks/Net10/net10_Nin0_Nout28_Qin0.0001.bin) has a companion CSV file called  [_net10_Nin0_Nout28_Qin0.0001.csv_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/networks/Net10/net10_Nin0_Nout28_Qin0.0001.csv). Each of these CSV files has 5 columns, reporting 5 network properties:
- **variable _vm_**: mean value of the blood velocity distribution, in `mm/s`;
- **variable _vs_**: standard deviation of the blood velocity distribution, in `mm/s`;
- **variable _anb_**: apparent network branching (ANB), expressed in `number of segments`. This is a count of the number of capillary segments that spins flowing through the network pass through during a given observation time (in this case, 100 ms);
- **variable _rm_**: mean capillary segment radius, in `µm`;
- **variable _lm_**: mean capillary segment length, in `µm`.

As an example, the content of one of these CSV files would look something like this:
```
vm,vs,anb,rm,lm
4.56,3.23,14.76,5.67,19.56
```

### Using vascular networks and their properties
We have written two tutorials to show you how to use these data for your modelling of vascular diffusion MRI signals.  
* This first [tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial1.md) **describes in detail a SpinFlowSim _pipenet_ object, and shows how to resolve a vascular network from capillaries drawn on histology**.
* This second [tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial2.md) shows **how synthetic signals from _pipenet_ objects can be used to inform microvascular parameter estimation**.  
