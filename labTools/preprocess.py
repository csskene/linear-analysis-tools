# Preprocess to get the required matrices
from weightMatrices import *
from scipy.io import loadmat
import os,sys,getopt
sys.path.append(os.environ["PETSC_DIR"]+'/lib/petsc/bin/')
import PetscBinaryIO
import argparse

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
    matDict = loadmat(linop)

    Base = matDict['Base'+baseName]
    OptL_R = matDict['OptL'+baseName]
    OptL_B = matDict['OptL_Beta'+baseName]

    GO = cons2Prim(Base)
    GI = prim2Cons(Base)
    FI,FO = EnormQuad(Base)
    WO = FO@GO
    WI = GI@FI

    if not os.path.exists(outdir):
        os.mkdir(outdir)

    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/L.dat','w'), OptL_R)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/WO.dat','w'), WO)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/WI.dat','w'), WI)
    PetscBinaryIO.PetscBinaryIO().writeMatSciPy(open(outdir+'/LB.dat','w'), OptL_B)

    f = open(outdir+'/info.txt', 'w')
    f.write('Binary files created from the linop %s \n' % linop)
    f.write('The PETSc arch is % s' % os.environ["PETSC_ARCH"])
    f.close()
