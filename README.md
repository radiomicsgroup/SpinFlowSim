<h1 align="center">SpinFlowSim: dMRI signal synthesis from spins flowing in capillary networks</h1>

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/overview.png" width="1300" height="auto">
</div>

</div>

**Introducing SpinFlowSim**! 


Here you'll find the code and essential data to simulate diffusion Magnetic Resonance Imaging (dMRI) signals arising from spins flowing in synthetic, realistic microvasculature networks obtained from histology.

The SpinFlowSim code was written by Francesco Grussu (<fgrussu@vhio.net>), Anna Voronova (<annavoronova@vhio.net>) and Athanasios Grigoriou (<agrigoriou@vhio.net>). SpinFlowSim is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/license.txt). 

**The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010"**.

Foster open science by citing our paper if you use SpinFlowSim in your research: 

_Anna Kira Voronova, Athanasios Grigoriou, Kinga Bernatowicz, Sara Simonetti, Garazi Serna, Núria Roson, Manuel Escobar, Maria Vieito, Paolo Nuciforo, Rodrigo Toledo, Elena Garralda, Els Fieremans, Dmitry S. Novikov, Marco Palombo, Raquel Perez-Lopez, Francesco Grussu._  
  **"SpinFlowSim: a blood flow simulation framework for histology-informed diffusion MRI microvasculature mapping in cancer"**   
  Medical Image Analysis, Volume 102, May 2025, 103531, ISSN 1361-8415, doi: [10.1016/j.media.2025.103531](https://doi.org/10.1016/j.media.2025.103531)

<div align="center">
    
<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/Media_screenshot.PNG" width="auto" height="auto">
</div>
    
<a href="#requirements">Requirements</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#installation">Installation</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#repository-description">Repository description</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#usage">Usage</a>&nbsp;&nbsp;&nbsp;
<br/><br/>


</div>


# Requirements

</div>
SpinFlowSim has been developed with python 3.10.8. To use SpinFlowSim and go through all its tutorials, the following third-party packages are required:

- [numpy](https://numpy.org) (developed with version 1.24.2)
- [pandas](https://pandas.pydata.org) (developed with version 1.5.3)
- [Lcapy](https://lcapy.readthedocs.io) (developed with version 1.10)
- [PySpice](https://github.com/FabriceSalvaire/PySpice) (developed with version 1.5)
- [graph-tool](https://graph-tool.skewed.de) (developed with version 2.45, commit b1a649d8)
- [nibabel](https://nipy.org/nibabel) (developed with version '5.1.0')
- [matplotlib](https://matplotlib.org) (examples tested with version '3.7.1')
- [mrtrix](https://www.mrtrix.org) (developed with version '3.0.4').

  
# Installation
</div>
Clone our repository as

 ```sh
 git clone https://github.com/radiomicsgroup/SpinFlowSim
 ```

To perform flow simulations with SpinFlowSim and to go through all the tutorails, you will need all the dependencies listed above. As a bear minimum, to simply synthesise vascular dMRI signals with the [`syn.py`](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/code/syn.py) tools, you will need [graph-tool](https://graph-tool.skewed.de), [Lcapy](https://lcapy.readthedocs.io), [PySpice](https://github.com/FabriceSalvaire/PySpice) and [numpy](https://numpy.org). With Anaconda, run:


 ```sh
conda create --name spinflowsim -c conda-forge graph-tool python=3.10.8
conda activate spinflowsim
conda install -n spinflowsim numpy
conda install -n spinflowsim -c conda-forge pyspice
pip install lcapy
 ```
</div>


# Repository description
The SpinFlowSim repository includes several sub-fodlers:
* The folder [**code**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/code) stores the actual simulation and visualisation tools behind SpinFlowSim:
    * the file [`syn.py`](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/code/syn.py) introduces the _syn_ module, where the _pipenet_ class is defined to work with vascular networks and to synthesise vascular dMRI signals. The manual of the _syn_ module can be found [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/syn_manual.md).
    * the file [`visu.py`](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/code/visu.py) intrdocudes the _visu_ module, where useful tools are defined to visualise spins flowing through vascular networks. The manual of the _visu_ module can be found [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/manuals/visu_manual.md).
* The folder [**examples**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples) contains some tutorials that illustrate how to use SpinFlowSim in practice.
* The folder [**networks**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks) contains the resolved vascular networks generated for our [paper](https://doi.org/10.1016/j.media.2025.103531). These can be used to synthesise vascular dMRI signals for any acquisition protocol of interest (find [here](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks/README.md) how).
    

# Usage

We have included some tutorials in the folder [**examples**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples) to illustrate how to use SpinFlowSim for your diffusion MRI analyses. 

* In [**tutorial 1**](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial1.md), we show **how to create and initialise a SpinFlowSim _pipenet_ object to resolve a vascular network drawn on histology, and to synthesise dMRI signals from the resolved vascular network**.
* In [**tutorial 2**](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial2.md), we show **how synthetic signals can be used to inform microvasculature property inference**, replicating the _in silico_ parameter estimation experiments of our [paper](https://doi.org/10.1016/j.media.2025.103531).



