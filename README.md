# SIMPLE POINT AND PANEL BEM IN PETSc #

This repository contains the code used to produce the results in "Work/Precision Tradeoffs in Continuum Models of Biomolecular Electrostatics", Matthew G. Knepley and Jaydeep P. Bardhan, ASME 2015 International Mechanical Engineering Congress & Exposition.

### Creating Figures ###

```
#!bash

make figures
```

### Building Code ###

You must first [Install PETSc](http://www.mcs.anl.gov/petsc/documentation/installation.html). Then the executable (testSrfOnSurfacePoints) can be built using
```
#!bash

make
```

### Running Regression ###

```
#!bash

make regression
```