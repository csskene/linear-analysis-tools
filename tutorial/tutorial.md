# Tutorial for the code (in progress)

This tutorial will walk you through using the code for the operator `OptL.mat` kindly provided by Chi-An Yeh.

## Preprocess
First of all we need to read `OptL.mat` and create binary files which can be used by PETSc. To do this run
```
python ../labTools/preprocess.py -linop OptL.mat -outputdir binFiles -baseName _NP01
```
This will create the binary files in the folder binFiles. The argument `-baseName _NP01` is used to tell the code that the files needed in `OptL.mat` have the name *base_NP01* and not just *base*.

## The main code
Now we have created the binary files we are ready to do a resolvent analysis. In series we can achieve this by running
```
python ../resl.py -linop binFiles/L.dat -wifile binFiles/WI.dat -wofile binFiles/WO.dat -wifile binFiles/WI.dat -omega 1 -outputdir out -saveLeading
```
This will run the resolvent for unit frequency and put the output in the folder */out*. The leading singular values will be saved. We can see that a lot of arguments were used. In order to make the code more readable these options can be written in a file. The file *opts.txt* contains the same arguments except that omega is changed to 2, and we specify that we want mumps to do the LU decomposition by adding `-pc_factor_mat_solver_type mumps`. Now run the resolvent using
```
python ../resl.py -options_file opts.txt
```
We see that the file `out/singularvalues.txt` has been appended to with the results for the new frequency.

## Postprocess

In order to change our binary files for the singular vectors into a MATLAB readable format we run
```
python ../labTools/postprocess.py -modesdir out -outputdir mfiles
```
This creates the file *singularVectors.mat* in the folder *mfiles* which contains MATLAB readable versions of the binary vectors contained in the folder *out*.
