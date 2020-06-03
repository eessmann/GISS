# GISS
GISS is a powerful and efficient 3D Direct Numerical Simulation (DNS) flow and solid solver to simulate solid-fuild flows at unprecedented detail and accuracy.


This solver has been developed by Erich Essmann (IMT, University of Edinburgh),Pei Shui (IMT, University of Edinburgh), Prashant Valluri (IMT, University of Edinburgh), Rama Govindarajan (ICTS, TIRF) and  ..... This project would not have been possible without finincal support from Thermosmart

The flow subsolver is based on Gerris which has been developed by Stéphane Popinet. 


Features
========
* Highly accurate DNS solver using Gerris 
* Open Dynamics Engine used for solid modelling
* Immerised Boundary Method for complex solid geometries
* High Particle Reynolds numbers
* Solid/Fluid density ratios
* Fluid mixing caused by solid motion


Sample Projects
===============

![Chaotic Ellipsoid](archive/images/Symmertic_Far_field.png)


Installation
============

Unix users (including Windows users under Cygwin):

```shell
./configure
make
make install
```

See the file 'INSTALL' for generic configure instructions and the tutorial
in doc/tutorial for an introduction on how GISS works.


