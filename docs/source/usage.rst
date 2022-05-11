Usage
=====

To keep things as general as possible the workflow is separated into three pieces.

Firstly, there is preprocessing code that writes linear operators to PETSc binary files that can be easily read by the main piece of code.

Secondly, there are two main codes resl.py and eigs.py for resolvent analysis and eigenvalue computations, respetively. This will be expanded upon in future updates.

Lastly, there is postprocessing code which transfers the PETSc binary output of the main code to a format that can be read by MATLAB/python.

Preprocessing
-------------

All required matrices must first be written to PETSc binary files so that the code and read these in as distributed sparse matrices. To do this the python script *preprocess.py* in the folder *tools* can be used.

Example: To obtain binary files from the linear operators in the file *cavity.mat* using the names *Base*, *OptL* and *OptL_Beta* run::

  python preprocess.py -linop path/to/cavity.mat

If the names are instead appended to *Base_NP01*, *OptL_NP01* and *OptL_Beta_NP01* for example instead run::

  python preprocess.py -linop path/to/cavity.mat -baseName _NP01

to append *'_NP01'* to the names. Lastly, to specify an output directory for the binary files use the argument::

  -outputdir outdir

Eigenvalue solver
-----------------
Documentation coming soon

Resolvent analysis
------------------

Given the linear operators, conducting a resolvent analysis requires two main ingredients; a linear systems solver, and an svd solver. Before discussing how to set each of these some specific arguments for choosing which resolvent modes to obtain are now given.

The linear operator file is set using::

  -linop filestr

The input weight matrix is set using::

  -wifile filestr

Likewise, the output weight matrix is set using::

  -wofile filestr

The angular frequency at which to compute the resolvent mode is set via the argument::

  -omega float

More than one frequency can be set using a string of frequencies is given with each one separated by a space::

  -omegas 'float float float ...'

An equispaced range of frequencies can be set using::

  -omegaRange 'startFreq endFreq numOfFreqs'

Discounting parameters can be set in the same way as the frequencies using any of::

  -disc float
  -discs 'float float float ...'
  -discRange 'startDisc endDisc numOfDiscs'

If a spanwise wavenumber is to be used this is set using any of::

  -beta float
  -betas 'float float float ...'
  -betaRange 'startBeta endBeta numOfBetas'

Note, if a spanwise wavenumber is set you must also specify the spanwise linear operator using::

  -linopbeta filestr

To save the leading singular vectors use::

  -outputdir 'outdir'

to specify an output directory. Then also use::

  -saveLeading

to save only the leading singular vector, or::

  -saveAll

to save all the singular vectors.

After each solve a file is created in the output directory called 'singularvalues.txt'. This contains all the singular values and can easily be read and sorted in python using the *pandas* library. For an example see the *plotGains.py* file provided.

Linear systems solves
^^^^^^^^^^^^^^^^^^^^^
Secondly, multiplying a vector by the resolvent requires a matrix inverse (which should ***never ever be explicitly computed***). Instead, these resolvent-vector multiplications are carried out by solving a linear system. The code is setup so that the default is to compute the LU decomposition of the inverse resolvent (which is a sparse matrix). This is done using a petsc KSP object. The method for carrying out this linear system solve can be overwritten using command-line arguments at runtime. Some examples of how this behaviour can be set are now provided.

Running the code using::

  python resl.py ...options

will compute the resolvent modes using the default settings in series. This will take the LU decomposition using the built in petsc method. If we instead wish to use mumps to take the LU decomposition this can be achieved (as long as PETSc was installed with mumps) by instead running::

  python resl.py -pc_factor_mat_solver_type mumps

To run the code in parallel, for example on four processors, run::

  mpiexec -n 4 python resl.py -pc_factor_mat_solver_type mumps

Note, the petsc LU decomposition cannot run in parallel so mumps or another parallel LU package must be installed (we can also omit the -pc_factor_mat_solver_type argument in this case).

Iterative solvers may also be specified using the argument `-ksp_type KSPType` and the preconditioner can be set using `-pc_type PCType` (see petsc manual for full details and for the options that KSPType and PCType can take). The code will automatically print some information about the linear system solver, however more information can be obtained by using the argument `-ksp_view`. For sanity, if you want to check that the solving the linear system and solving the corresponding linear hermitian system are truly the hermitian transpose of each other use the argument `-test_adj`.

Singular value decomposition
^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Now that the resolvent matrix is setup we can compute the resolvent modes by taking the singular value decomposition of the resolvent matrix. This is done using the slepc svd object. As for the linear system solver there are numerous choices over svd solver and some useful ones are now provided (see the slepc manual for full details)

Number of singular triplets to compute::

  -svd_nsv k

Maximum dimension of the subspace for the solver::

  -svd_ncv m

Set the SVD solver::

  -svd_type SVDType

Set the convergence tolerance::

  -svd_tol tol

To obtain full information about the SVD solver use the argument `-svd_view`.

Randomised SVD
^^^^^^^^^^^^^^

Documentation coming soon

Postprocessing
--------------
PETSc outputs binary files. For further analysis/visualisation in python/MATLAB the python script *postprocess.py* can be used to convert the data from the binary files to .dat files.

Running::

  python postprocess.py -modesdir modesdir

will write all the singular vectors contained in *modesdir* to a MATLAB file. To specify python output use the argument::

  -python True

An output directory can be specified using::

  -outputdir outdir


Plotting the gains
^^^^^^^^^^^^^^^^^^
The gains can be plotted using the script *plotGains.py*. This script uses the gradient information outputed in the singularvalues.txt to interpolate between consecutive gains with cubic-Hermite splines. The code needs the python library `pandas` in order to sort through the saved data. Install this using conda or pip.

To use this code run::

  python plotGains.py -gainsFile /path/to/singularvalues.txt

To specify a mode to plot use the argument::

  -mode n

Note, that n=0 is the leading singular value. Specify n=-1 to plot all modes. The figure will be saved as a .png file.
