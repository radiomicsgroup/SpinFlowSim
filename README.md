<h1 align="center">SpinFlowSim: synthesize dMRI signals from spins flowing in capillary networks</h1>


![Resolvednetworks](https://github.com/user-attachments/assets/e7386d51-27aa-4650-8e75-84545a530daf)

</div>

**Introducing SpinFlowSim**! 


Here you'll find the code and essential data to simulate diffusion Magnetic Resonance Imaging (dMRI) signals arising from spins flowing in synthetic, realistic microvasculature networks obtained from histology.

The SpinFlowSim code was written by Francesco Grussu (<fgrussu@vhio.net>), Anna Voronova (<annavoronova@vhio.net>) and Athanasios Grigoriou (<agrigoriou@vhio.net>). SpinFlowSim is distributed under the **Attribution-NonCommercial-ShareAlike 4.0 International license** ([CC BY-NC-SA 4.0](https://creativecommons.org/licenses/by-nc-sa/4.0)). Copyright (c) 2024, Fundació Privada Institut d’Investigació Oncològica de Vall d’Hebron (Vall d'Hebron Institute of Oncology (VHIO), Barcelona, Spain). All rights reserved. Link to license [here](https://github.com/radiomicsgroup/SpinFlowSim/blob/main/license). 

**The project that gave rise to these results received the support of a fellowship from ”la Caixa” Foundation (ID 100010434). The fellowship code is "LCF/BQ/PR22/11920010"**.

Foster open science by citing our preprint if you use SpinFlowSim in your research: 

Anna Voronova, Athanasios Grigoriou, Kinga Bernatowicz, Sara Simonetti, Garazi Serna, Nuria Roson, Manuel Escobar, Maria Vieito, Paolo Nuciforo, Rodrigo Toledo, Elena Garralda, Roser Sala-Llonch, Els Fieremans, Dmitry S. Novikov, Marco Palombo, Raquel Perez-Lopez, Francesco Grussu. **"SpinFlowSim: a blood flow simulation framework for histology-informed diffusion MRI microvasculature mapping in cancer"**. medRxiv 2024.07.15.24310335; doi: https://doi.org/10.1101/2024.07.15.24310335.

<div align="center">
    
![qr_img](https://github.com/user-attachments/assets/c4c9c69d-48c6-405e-8837-b3afde524312)
    
<a href="#requirements">Requirements</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#installation">Installation</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#usage">Usage</a>&nbsp;&nbsp;&nbsp;|&nbsp;&nbsp;&nbsp;<a href="#other-examples">Other examples</a>&nbsp;&nbsp;&nbsp;
<br/><br/>


</div>


# Requirements

</div>
SpinFlowSim has been developed with python 3.10.8.

It requires the following third-party packages:

- [pandas](https://pandas.pydata.org/) (developed with version 1.5.3)
- [lcapy](https://lcapy.readthedocs.io/en/latest) (developed with lcapy 1.10)
- [PySpice](https://github.com/FabriceSalvaire/PySpice) (developed with PySpice 1.5)
- [graph-tool](https://graph-tool.skewed.de) (developed with graph-tool 2.45, commit b1a649d8)

# Installation
</div>
Obtain SpinFlowSim by cloning this repository:

 ```sh
 git clone https://github.com/radiomicsgroup/SpinFlowSim
 ```
</div>
Start using SpinFlowSim by activating the pipenet environment

```
$:- conda create --name pipenet  -c conda-forge graph-tool  python=3.10.8
$:- conda activate pipenet
```
Install the required packages

```
$:- pip install -r requirements.txt
```
# Usage

Under construction - several tutorials will be available soon!




