# Convert petsc binary vecs to numpy arrays
import sys,os
import numpy as np
from scipy.io import savemat
sys.path.append(os.environ["PETSC_DIR"]+'/lib/petsc/bin/')
import PetscBinaryIO
import argparse

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument('-modesdir',type=str,help='Folder where the singular vector binary files are.')
    parser.add_argument('-outputdir',type=str,help='Folder for output files')
    parser.add_argument('-python',type=bool,help='To get python output rather than MATLAB')
    args = parser.parse_args()

    if(args.python):
        flg_py = True
        print('Writing python files')
    else:
        flg_py=False
        print('Writing MATLAB files')
    if(args.outputdir):
        outdir = args.outputdir
        print('writing to folder ',outdir)
        if not os.path.exists(outdir):
            os.mkdir(outdir)
    else:
        outdir = './'
        print('writing to default folder ',outdir)

    args = parser.parse_args()
    modesdir = args.modesdir+'/'
    files = os.listdir(modesdir)
    modeFiles = [file for file in files if file.endswith('.dat')]
    mdict   = {}
    for modeFile in modeFiles:
        vec = PetscBinaryIO.PetscBinaryIO().readBinaryFile(modesdir+modeFile)
        if(flg_py):
            vec = np.squeeze(np.array(vec))
            np.save(outdir+'/'+modeFile[:-4],vec)
        else:
            modeFile = modeFile.replace('.','p')
            mdict[modeFile[:-4]] = vec

    if(not flg_py):
        savemat(outdir+'/singularVectors.mat',mdict,appendmat=True)
