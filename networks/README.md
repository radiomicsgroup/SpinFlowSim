# Capillary networks for vascular dMRI signal synthesis

SpinFlowSim was showcased on vascular networks reconstructed from 2D histological images. The networks were drawn manually, tracing visible capillaries in 11 stained liver biopsies.

The SpinFlowSim dataset consists of **15 vascular networks**:

![networks](https://github.com/user-attachments/assets/029fdf41-a655-45ef-8aee-d61c3860b416)

The network dataset has been augmented to 1500 realizations by varying input flow in [1.5 × 10⁻⁴ ; 5.5 × 10⁻³] mm³/s, as well as the inlet/outlet nodes.



You can start running SpinFlowSim by choosing from any of the 1500 network realizations (in binary `.bin` format) stored within this folder.

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
Resolved vascular networks are stored as binary (.bin) files, containing instantiations of the class _pipenet_ (_i.e., _pipenet_ objects; the _pipenet_ class is defined in [_syn.py_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/syn.py)). 

A file name `net10_Nin0_Nout28_Qin0.0001.bin` indicates that:
- the file is located in the folder corresponding to the `Net10` vascular network;
- the network ID is `10`;
- the flow inlet node is `0`, while the the outlet node is `28`;
- the input volumetric flow rate is `1 × 10⁻⁴ mm³/s`.

Attributes and methods of objects of the class _pipenet_ are described in detail in this [tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial1.md).

### Network properties (.csv files)
Each _pipenet_ object storing a network is accompanied by a small CSV file bearing the same name, and storing summary network properties. For example, the network file  [_net10_Nin0_Nout28_Qin0.0001.bin_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/networks/Net10/net10_Nin0_Nout28_Qin0.0001.bin) has a companion CSV file called  [_net10_Nin0_Nout28_Qin0.0001.csv_](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/networks/Net10/net10_Nin0_Nout28_Qin0.0001.csv). Each of these CSV files has 5 columns, reporting 5 network properties:
 
