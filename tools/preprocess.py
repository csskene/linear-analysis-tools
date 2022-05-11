# Preprocess to get the required matrices
from weightMatrices import *
from scipy.io import loadmat
import os,sys,getopt
sys.path.append(os.environ["PETSC_DIR"]+'/lib/petsc/bin/')
import PetscBinaryIO
import argparse
import h5py
from scipy import sparse
import numpy as np

if __name__ == '__main__':

    parser = argparse.ArgumentParser()
    parser.add_argument('-linop',type=str,help='Location of the linear operator')
    parser.add_argument('-baseName',type=str,help='Name of the base')
    parser.add_argument('-outputdir',type=str,help='Output folder for the binary files')
    args = parser.parse_args()

    if(args.outputdir):
        outdir = args.outputdir
        print('Writing files to specified folder, ',  outdir)
    else:
        print('No output folder specified')
        outdir = './'

    if(args.baseName):
        baseName = args.baseName
        print('Reading the linop using the extension Base', baseName)
    else:
        print('Will read the linop using the default basename: Base')
        baseName = ''

    if(args.linop):
        print('Reading linop from',args.linop)
        linop = args.linop
    else:
        print('You must specify the path to the linear operator')

    # Load the linear operator and create the required matices
    try:
        matDict = loadmat(linop)
        Base = matDict['Base'+baseName]
        OptL_R = matDict['OptL'+baseName]
        OptL_B = matDict['OptL_Beta'+baseName]
    except:
        # For matlab 7.3 files
        print('Loading using h5py')
        matDict =  h5py.File(linop, 'r')
        Base = np.array(matDict['Base'+baseName]).T
        OptL_R = sparse.csc_matrix((matDict['OptL'+baseName]["data"], matDict['OptL'+baseName]["ir"], matDict['OptL'+baseName]["jc"]))
        data = matDict['OptL_Beta'+baseName]["data"]['real']+1j*matDict['OptL_Beta'+baseName]["data"]['imag']
        OptL_B = sparse.csc_matrix((data, matDict['OptL_Beta'+baseName]["ir"], matDict['OptL_Beta'+baseName]["jc"]))


    sponge = np.array([1, -1, 1, -1], dtype = float)
    GO = cons2Prim(Base)
    GI = prim2Cons(Base)
    FI,FO = EnormQuad(Base)
    WO = FO@GO
    WI = GI@FI
    SP = linopSponge(Base,sponge)
    phi = physVec(Base)

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/L.dat','w'), OptL_R)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/WO.dat','w'), WO)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/WI.dat','w'), WI)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/LB.dat','w'), OptL_B)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/SP.dat','w'), SP)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/PHI.dat','w'), phi)

    f = open(outdir+'/info.txt', 'w')
    f.write('Binary files created from the linop %s \n' % linop)
    f.write('The PETSc arch is % s' % os.environ["PETSC_ARCH"])
    f.close()
