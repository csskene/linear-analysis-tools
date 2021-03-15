# Linear stability analysis code specifically for Taira Lab

Scalable code for conducting linear analyses using the python bindings of PETSc and SLEPc. Inspired by some code written by [Miguel Fosas de Pando](https://github.com/miguelfp/nthRootsOfUnity).

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

## Installation instructions
Probably the hardest part of using this code is the installation. Here are some instructions but if there are any problems please contact me and I can try to help.

Good luck! :four_leaf_clover:

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
SLEPc is best downloaded as a tarball from *https://slepc.upv.es/download/* (the gitlab repository is the development version and not the release version). From this directory set the following environmental variables
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
# Usage
To keep things as general as possible the workflow is separated into three pieces.

Firstly, there is preprocessing code that writes the linear operators commonly used by our lab to PETSc binary files that can be easily read by the main piece of code.

Secondly, there is the main linear analysis code (which so far can only do resolvent analysis). This will be expanded upon in future updates.

Lastly, there is postprocessing which transfers the PETSc binary output of the main code to a format that can be read by MATLAB/python.

## Preprocessing
The linear operator (and weight matrices) must first be written to PETSc binary files so that the code and read these in as distributed sparse matrices. To do this the python script *preprocess.py* in the folder *labTools* can be used.

Example: To obtain binary files from the linear operators in the file *cavity.mat* using the names *Base*, *OptL* and *OptL_Beta* run
```
python preprocess.py -linop path/to/cavity.mat
```
If the names are instead appended to *Base_NP01*, *OptL_NP01* and *OptL_Beta_NP01* for example instead run
```
python preprocess.py -linop path/to/cavity.mat -baseName _NP01
```
to append *'_NP01'* to the names. Lastly, to specify an output directory for the binary files use the argument
```
-outputdir outdir
```

## The main code
### Resolvent analysis
Given the linear operators, conducting a resolvent analysis requires two main ingredients; a linear systems solver and an svd solver. Before discussing how to set each of these some specific arguments for choosing which resolvent modes to obtain are now given.

The linear operator file is set using
```
-linop filestr
```
The input weight matrix is set using
```
-wifile filestr
```
Likewise, the output weight matrix is set using
```
-wofile filestr
```
The angular frequency at which to compute the resolvent mode is set via the argument
```
-omega float
```
More than one frequency can be set using a string of frequencies is given with each one separated by a space.
```
-omegas 'float float float ...'
```
An equispaced range of frequencies can be set using
```
-omegaRange 'startFreq endFreq numOfFreqs'
```

Discounting parameters can be set in the same way as the frequencies using any of
```
-disc float
-discs 'float float float ...'
-discRange 'startDisc endDisc numOfDiscs'
```
If a spanwise wavenumber is to be used this is set using any of
```
-beta float
-betas 'float float float ...'
-betaRange 'startBeta endBeta numOfBetas'
```
Note, if a spanwise wavenumber is set you must also specify the spanwise linear operator using
```
-linopbeta filestr
```
To save the leading singular vectors use
```
-outputdir 'outdir'
```
to specify an output directory. Then also use
```
-saveLeading
```
to save only the leading singular vector, or
```
-saveAll
```
to save all the singular vectors.

After each solve a file is created in the output directory called 'singularvalues.txt'. This contains all the singular values and can easily be read
and sorted in python using the *pandas* library. For an example see the *plotGains.py*
file provided.

#### Linear systems solves
Secondly, multiplying a vector by the resolvent requires a matrix inverse (which should ***never ever be explicitly computed***). Instead, these resolvent-vector multiplications are carried out by solving a linear system. The code is setup so that the default is to compute the LU decomposition of the inverse resolvent (which is a sparse matrix). This is done using a petsc KSP object. The method for carrying out this linear system solve can be overwritten using command-line arguments at runtime. Some examples of how this behaviour can be set are now provided.

Running the code using
```
python resl.py
```
will compute the resolvent modes using the default settings in series. This will take the LU decomposition using the built in petsc method. If we instead wish to use mumps to take the LU decomposition this can be achieved (as long as PETSc was installed with mumps) by instead running
```
python resl.py -pc_factor_mat_solver_type mumps
```
To run the code in parallel on four processors run
```
mpiexec -n 4 python resl.py -pc_factor_mat_solver_type mumps
```
Note, the petsc LU decomposition cannot run in parallel so mumps or another parallel LU package must be installed (we can also omit the -pc_factor_mat_solver_type argument in this case).

Iterative solvers may also be specified using the argument `-ksp_type KSPType` and the preconditioner can be set using `-pc_type PCType` (see petsc manual for full details and for the options that KSPType and PCType can take). The code will automatically print some information about the linear system solver, however more information can be obtained by using the argument `-ksp_view`. For sanity, if you want to check that the solving the linear system and solving the corresponding linear hermitian system are truly the hermitian transpose of each other use the argument `-test_adj`.

#### Singular value decomposition
Now that the resolvent matrix is setup we can compute the resolvent modes by taking the singular value decomposition of the resolvent matrix. This is done using the slepc svd object. As for the linear system solver there are numerous choices over svd solver and some useful ones are now provided (see the slepc manual for full details)

``-svd_nsv k
``
: number of singular triplets to compute\
``-svd_ncv m``
: maximum dimension of the subspace for the solver\
``-svd_type SVDType``
: sets the SVD solver\
``-svd_tol tol``
: sets the convergence tolerance

To obtain full information about the SVD solver use the argument `-svd_view`.

## Postprocessing
PETSc outputs binary files. For further analysis/visualisation in python/MATLAB the python script *postprocess.py* can be used to convert the data from the binary files to .dat files.

Running
```
python postprocess.py -modesdir modesdir
```
will write all the singular vectors contained in *modesdir* to a MATLAB file. To specify python output use the argument
```
-python True
```
An output directory can be specified using
```
-outputdir outdir
```

## Plotting the gains
The gains can be plotted using the script *plotGains.py*. This script uses the gradient information outputed in the singularvalues.txt to interpolate between consecutive gains with cubic-Hermite splines. The code needs the python library `pandas`. Install this using conda or pip.

To use this code run
```
python plotGains.py -gainsFile /path/to/singularvalues.txt
```
To specify a mode to plot use the argument
```
-mode n
```
Note, that n=0 is the leading singular value. Specify n=-1 to plot all modes. The figure will be saved as a .png file.

## Tips
There could be many options to be set to run this code. Luckily, these can all be written in a file and can be used by running the code with the option `-options_file`. For example, if the code options are in *opts.txt* the code can be used with these arguments by running
```
mpiexec python resl.py -options_file opts.txt
```
For a comprehensive breakdown of the code use the argument `-log_view`.
## Concluding remarks
I hope this code provides valuable to you all. Any comments, suggestions, or things you want implemented please let me know.

:seedling:
