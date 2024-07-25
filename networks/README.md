
SpinFlowSim is showcased on vascular networks reconstructed from 2D histological images.  

The networks were drawn manually, tracing visible capillaries in 11 stained liver biopsies.

The SpinFlowSim dataset consists of **15 vascular networks**:

![networks](https://github.com/user-attachments/assets/029fdf41-a655-45ef-8aee-d61c3860b416)

The network dataset has been augmented to 1500 realizations by varying input flow [1.5 × 10⁻⁴ ; 5.5 × 10⁻³] mm³/s, as well as inlet/outlet nodes.



Start running SpinFlowSim by choosing from any of the 1500 network realizations `.bin` realizations in [networks](https://github.com/annavoronova/SpinFlowSim/tree/main/networks):

### Directory Description:

- **Net{ID}**: Directory containing the files, where `{ID}` is the unique identifier for the network.

### Sub-Directory Filename Description:
- **net{ID}**: Base name of the file, where `{ID}` corresponds to the same unique identifier used in the directory.
- **Nin{input_node}**: The flow input node in the given network realization.
- **Nout{output_node}**: The flow output node in the given network realization.
- **Qin{input_volumetric_flow_rate}**: Represents the input volumetric flow rate in mm³/s.

### Example:

For a file named `net10_Nin0_Nout28_Qin0.0001.bin`, this indicates:

- The file is located in the `Net10` directory.
- The network ID is `10`.
- The flow inlet node is `0`, the outlet node is `28`.
- The input volumetric flow rate is `1 × 10⁻⁴ mm³/s`.



