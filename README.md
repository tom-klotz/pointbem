# SIMPLE POINT AND PANEL BEM IN PETSc #

This repository contains the code used to produce the results in "Work/Precision Tradeoffs in Continuum Models of Biomolecular Electrostatics", Matthew G. Knepley and Jaydeep P. Bardhan, ASME 2015 International Mechanical Engineering Congress & Exposition.

### Creating Figures ###

```
make figures
```

### Building Code ###

You must first [Install PETSc](http://www.mcs.anl.gov/petsc/documentation/installation.html).
Next, install MSMS (http://mgltools.scripps.edu/packages/MSMS).
Then, set your msms path to the msms executable via:
```
export MSMS_PATH=/your_msms_path/msms_executable
```

To make the srf files used in the examples, run (in the root folder)
```
make meshmaker
make srf
```

Then the executable (testSrfOnSurfacePoints) can be built using
```
cd src
make
```

### Running Regression ###

To run the regression, first set the relative path from your petsc directory to the pointbem directory via:
```
export PETSC_TO_MSMS=../path-to-msms/your-msms-folder
```

Then the regression can be run using
```
cd src
make regression
```