<h1 align="center">SpinFlowSim: synthesize dMRI signals from spins flowing in capillary networks</h1>

<div align="center">
  <img src="https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/imgs/overview.png" alt="MRHistoillustration" width="900" height="auto">
</div>

</div>

**Introducing SpinFlowSim**! 


Here you'll find the code and essential data to simulate diffusion Magnetic Resonance Imaging (dMRI) signals arising from spins flowing in synthetic, realistic microvasculature networks obtained from histology.

The SpinFlowSim code was written by Francesco Grussu (<fgrussu@vhio.net>), Anna Voronova (<annavoronova@vhio.net>) and Athanasios Grigoriou (<agrigoriou@vhio.net>). SpinFlowSim is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/license.txt). 

**The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010"**.

Foster open science by citing our preprint if you use SpinFlowSim in your research: 

_Anna Voronova, Athanasios Grigoriou, Kinga Bernatowicz, Sara Simonetti, Garazi Serna, Nuria Roson, Manuel Escobar, Maria Vieito, Paolo Nuciforo, Rodrigo Toledo, Elena Garralda, Roser Sala-Llonch, Els Fieremans, Dmitry S. Novikov, Marco Palombo, Raquel Perez-Lopez, Francesco Grussu. **"SpinFlowSim: a blood flow simulation framework for histology-informed diffusion MRI microvasculature mapping in cancer"**. medRxiv 2024.07.15.24310335; doi: [10.1101/2024.07.15.24310335](https://doi.org/10.1101/2024.07.15.24310335)_.

<div align="center">
    
![qr_img](https://github.com/user-attachments/assets/c4c9c69d-48c6-405e-8837-b3afde524312)
    
<a href="#requirements">Requirements</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#installation">Installation</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#repository-description">Repository description</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#usage">Usage</a>&nbsp;&nbsp;&nbsp;
<br/><br/>


</div>


# Requirements

</div>
SpinFlowSim has been developed with python 3.10.8. It requires the following third-party packages:

- [numpy](https://numpy.org) (developed with version 1.24.2)
- [pandas](https://pandas.pydata.org) (developed with version 1.5.3)
- [Lcapy](https://lcapy.readthedocs.io) (developed with Lcapy 1.10)
- [PySpice](https://github.com/FabriceSalvaire/PySpice) (developed with PySpice 1.5)
- [graph-tool](https://graph-tool.skewed.de) (developed with graph-tool 2.45, commit b1a649d8)

# Installation
</div>
Obtain SpinFlowSim by cloning this repository:

 ```sh
 git clone https://github.com/radiomicsgroup/SpinFlowSim
 ```
</div>
Start using SpinFlowSim by creating a new environment:

```
$:- conda create --name pipenet  -c conda-forge graph-tool  python=3.10.8
$:- conda activate pipenet
```
Install the required packages

```
$:- pip install -r requirements.txt
```

# Repository description
The SpinFlowSim repository includes several sub-fodlers:
* The folder [**code**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/code) stores the actual simulation and visualisation tools behind SpinFlowSim:
    * the file [`syn.py`](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/code/syn.py) defines the _pipenet_ class, with which we represent and resolve vascular networks. An initialised _pipenet_ object can be used to synthesise microvascular signals for any dMRI protocol of interest;
    * the file [`visu.py`](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/code/visu.py) defines useful tools to visualise spins flowing within a vascular network as a video.
* The folder [**examples**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples) contains some tutorials that illustrate how to use SpinFlowSim in practice.
* The folder [**networks**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks) contains the resolved vascular networks generated for our [preprint](https://doi.org/10.1101/2024.07.15.24310335). These can be used to synthesise vascular dMRI signals for any acquisition protocol of interest (find [here](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/networks/README.md) how).
    

# Usage

We have included some tutorials in the folder [**examples**](https://github.com/radiomicsgroup/SpinFlowSim/tree/main/examples) to illustrate how to use SpinFlowSim for your diffusion MRI analyses. 

* In this first [tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial1.md), we show **how to create and initialise a SpinFlowSim _pipenet_ object to resolve a vascular network drawn on histology**.
* In this second [tutorial](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/examples/tutorial2.md), we show **how synthetic signals can be used to inform microvasculature property inference**, replicating the _in silico_ parameter estimation experiments of our [preprint](https://doi.org/10.1101/2024.07.15.24310335).
    
**UNDER CONSTRUCTION! Please note that we are still polishing the SpinFlowSim code: the file structure may change slightly and some files may be still missing. Thanks for your patience!**


