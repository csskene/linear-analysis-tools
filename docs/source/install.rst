Installation
============

Probably the hardest part of using this code is the installation. Here are some instructions but if there are any problems please contact me and I can try to help.

Good luck! :four_leaf_clover:

## PETSc and SLEPc requirements
```
MPI
PETSc built with complex scalars
SLEPc built with complex scalars
```
## Python packages
All using Python 3.7 :snake:
```
numpy
scipy
h5py
mpi4py
petsc4py
slepc4py
pandas (optional)
```

### PETSc
PETSc can be downloaded from the gitlab repository by running
```
git clone -b release https://gitlab.com/petsc/petsc.git petsc
```

For detailed installations instructions please see *https://www.mcs.anl.gov/petsc/documentation/installation.html*

Here are some basic installation hints that worked for me to install a complex version of PETSc with debugging turned on (default). First make sure that the MPI is the same MPI that mpi4py uses is set to the path.
From the directory where PETSc is download run
```
./configure PETSC_ARCH=arch-complex-dbg --with-scalar-type=complex
```
and follow the terminal instructions (this may take a while :tea:). Note, PETSC_ARCH is your choice and multiple versions of PETSc can be simultaneously installed by specifying different versions. I usually also add the argument
```
--prefix=/path/to/folder
```
to specify the directory (a different one for each PETSC_ARCH) where I want PETSc to install the files.
#### Other external libraries
PETSc is able to download and install extra external libraries for you. These can be specified at configure time. For example, to install MUMPS add the following arguments to the configure
```
--download-mumps --download-scalapack --download-parmetis --download-metis --download-ptscotch
```
A full list of solvers is provided at this link *https://www.mcs.anl.gov/petsc/documentation/linearsolvertable.html*.
### SLEPc
SLEPc can be downloaded in a similar manner to PETSc by running
```
git clone -b release https://gitlab.com/slepc/slepc slepc
```
From this directory set the following environmental variables
```
export PETSC_DIR={where you put petsc}
export PETSC_ARCH=arch-complex-dbg (or whatever value you need)
export SLEPC_DIR={where you put slepc}
```
Now cd to where you downloaded SLEPc (``cd $SLEPC_DIR``) and run
```
./configure
```
and follow the terminal instructions. This will install a version of SLEPc that works with your version of PETSc that you installed with PETSC_ARCH.
## Python packages
### Anaconda
If you are using anaconda it is probably best to make a new environment. Make sure that the python version is set to 3.7 (3.8 does not currently work). You can then install numpy and scipy using conda. However, install mpi4py using pip as the mpi that PETSc and SLEPc are built with must match the mpi that mpi4py uses (using conda to install mpi4py installs its own mpi which will cause issues).
### Non Anaconda
Install numpy, scipy and mpi4py using pip.

### petsc4py
petsc4py should also be installed using pip, see *https://www.mcs.anl.gov/petsc/petsc4py-current/docs/usrman/install.html*

Make sure that as when installing slepc you set
```
export PETSC_DIR={where you put petsc}
export PETSC_ARCH=arch-complex-dbg (or whatever value you need)
```
Note, replace the export PETSC_ARCH with
```
export PETSC_ARCH=arch-complex-dbg:arch-complex-dbg:arch-whatever:...
```
to install petsc4py with multiple different versions of petsc concurrently.
Now run
```
pip install petsc4py
```
### slepc4py
After petsc4py has installed make sure that SLEPC_DIR is set via
```
export SLEPC_DIR={where you put slepc}
```
and install slepc4py using pip
```
pip install slepc4py
```
Note, if you need to add versions of petsc to petsc4py I've found that uninstalling the packages using
```
pip uninstall petsc4py
pip uninstall slepc4py
```
and then reinstalling by running
```
pip install petsc4py --no-cache-dir
pip install slepc4py --no-cache-dir
```
works well.
Hopefully with some perseverance :computer: :wrench: you should have everything installed to use the code!
