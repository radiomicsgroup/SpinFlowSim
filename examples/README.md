# _SpinFlowSim_ examples

This folder contains examples that illustrate how to use SpinFlowSim for your diffusion MRI analyses. The examples rely on additional custom-written code, which is included here for convenience. Some third-party packages are also required, namely:
* [numpy](https://numpy.org) (developed with version 1.24.2)
* [pandas](https://pandas.pydata.org) (developed with version 1.5.3)
* [Lcapy](https://lcapy.readthedocs.io) (developed with version 1.10)
* [PySpice](https://github.com/FabriceSalvaire/PySpice) (developed with version 1.5)
* [graph-tool](https://graph-tool.skewed.de) (developed with version 2.45, commit b1a649d8)
* [nibabel](https://nipy.org/nibabel) (developed with version '5.1.0')
* [matplotlib](https://matplotlib.org) (examples tested with version '3.7.1')
* [mrtrix](https://www.mrtrix.org) (developed with version '3.0.4').

In this first [tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial1.md), we show how to create and initialise a SpinFlowSim _pipenet_ object, in order to resolve an exemplificative vascular network drawn on histology.

In this second [tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial2.md), we show on _in silico_ data how synthetic signals can be used to inform microvasculature parameter estimation.
