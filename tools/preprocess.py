# Preprocess to get the required matrices
from scipy.io import loadmat
import os,sys
sys.path.append(os.environ["PETSC_DIR"]+'/lib/petsc/bin/')
import PetscBinaryIO
import argparse
from scipy import sparse
import numpy as np
from pathlib import Path

if __name__ == '__main__':
    """
    Reads a scipy sparse matrix and outputs a petsc binary file
    """
    parser = argparse.ArgumentParser()
    parser.add_argument('-linop',type=str,help='Location of the linear operator')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary file')
    args = parser.parse_args()

    if(args.outputdir):
        outdir = args.outputdir
        print('Writing files to specified folder, ',  outdir)
    else:
        print('No output folder specified')
        outdir = './'

    if(args.linop):
        print('Reading linop from',args.linop)
        linop = args.linop
    else:
        print('You must specify the path to the linear operator')

    opName = Path(linop).stem

    op = sparse.load_npz(linop)
    io = PetscBinaryIO.PetscBinaryIO()
    io.writeMatSciPy(open(outdir+'/'+opName+'.dat','w'), op)

    f = open(outdir+'/info.txt', 'w')
    f.write('Binary files created from the linop %s \n' % linop)
    f.write('The PETSc arch is % s' % os.environ["PETSC_ARCH"])
    f.close()
