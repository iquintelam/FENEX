FENEX
==============================
[//]: # (Badges)
[![GitHub Actions Build Status](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/FENEX/workflows/CI/badge.svg)](https://github.com/REPLACE_WITH_OWNER_ACCOUNT/FENEX/actions?query=workflow%3ACI)
[![codecov](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/FENEX/branch/master/graph/badge.svg)](https://codecov.io/gh/REPLACE_WITH_OWNER_ACCOUNT/FENEX/branch/master)
[![DOI](https://doi.org/10.1063/1.4866764.svg)](https://doi.org/10.1063/1.4866764)


Python implementation of the Free energy extrapolation method (FENEX) for estimating free energies at phase coexistence from polynomial models 

### Instalation

The development version can be installed from this repository:
```bash
git clone https://github.com/iquintelam/FENEX
cd FENEX
pip install
```
### Usage

Import "FENEX" and use simulation or experimental data to estimate the next coexistence point in the thermodynamic integration and calculate free energies. Coexistence properties can also be refined from the “near coexistence” simulated data. According to the [reference paper](https://doi.org/10.1063/1.5006047), one can choose one of two types of integration.

Suppose we want to integrate on a property that can be linearly decoupled from the Halmitonion (type 2).
Read data from the test system, which corresponds to an integration on the pair potential energy parameter($\epsilon$) to find the coexistence pressure(p) between a crystalline and an isotropic phase.

```python
>>> from FENEX import read_test_system
>>> Npoints,f1new,f,free_energy,z,cov,stats = read_test_system() 
```
We have the number of integration points, the next $\epsilon$ in the coexistence line `f1new`, the previous integration points ($\epsilon$,p) for both phases concatenated in an array of shape (2,Npoints,iphase). The free energy of the first point in the integration is in an array of shape (Npoints,iphase) `free_energy[:,0]`, the ensemble average of conjugate variables `z` (u,v) concatenated in an array of shape (2,Npoints,iphase), covariances (cov(u,U),cov(v,V),cov(u,V)] `cov` of previous points for both phases concatenated in an array of shape (3,npoints,iphase). Simulation statistics (acceptance probability of perturbation, potential energy) `stats` for both phases concatenated in an array of shape (2,npoints,iphase). This array is a non-default argument that can be passed to refine the coexistence values. 
We also need to define the integration type that can be 'coupled' or 'decoupled'.
To calculate the next point in the integration and the free energies of the previous point, we call `estimate_coexistence` from the class Integration. For that, we need to initialize the class Integration with all the data from the simulation as attributes:
```python
>>> int_1 = FENEX.Integrate(f1new,f,free_energy,z,cov,'coupled',stats)
>>> FENEX.Integrate.estimate_coexistence(int_1) 
>>> print(int_1.free_energy,int_1.f2new)
```
The coexistence pressure correspondent to the next `f2new` can be used as an input in a new simulation to continue the integration. The free energy of both phases from the previous coexistence points `free_energy` are also calculated.
To refine the coexistence properties, we call `refine_coex`. This function calculates the refined properties: `zsat`,`free_energy_sat`,`enesat`,`f2sat`).
```python
>>> results_sat = FENEX.Integrate.refine_coexistence(int_1)
>>> print(int_1.free_energy_sat)
```
### References

Please cite the original FENEX papers: 
* Escobedo, F. A. J. Chem. Phys. 140, 094102 (2014). [DOI](https://doi.org/10.1063/1.4866764)
* Escobedo, F. A. J. Chem. Phys. 147, 214501 (2017). [DOI](https://doi.org/10.1063/1.5006047)

### Copyright

Copyright (c) 2022, Isabela Quintela Matos


#### Acknowledgements
 
Project based on the 
[Computational Molecular Science Python Cookiecutter](https://github.com/molssi/cookiecutter-cms) version 1.6.
