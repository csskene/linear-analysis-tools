import sys,os,datetime
import numpy as np

import os.path
from os import path

# Append petsc and slepc to the path
PETSC_DIR  = os.getenv('PETSC_DIR')
PETSC_ARCH = os.getenv('PETSC_ARCH')
SLEPC_DIR = os.getenv('SLEPC_DIR')
petscPath = PETSC_DIR + '/' + PETSC_ARCH + '/lib'
slepcPath = SLEPC_DIR + '/' + PETSC_ARCH + '/lib'
sys.path.append(petscPath)
sys.path.append(slepcPath)

import petsc4py
import slepc4py

petsc_arch = os.environ["PETSC_ARCH"]
slepc4py.init(sys.argv,arch=petsc_arch)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print


if __name__ == '__main__':
    # Handles all options and running the relevant codes
    opts = PETSc.Options()
    now = datetime.datetime.now()
    comm = PETSc.COMM_WORLD
    size = comm.getSize()
    rank = comm.getRank()

    Print (now.strftime("%Y-%m-%d %H:%M:%S"))
    if(size==1):
        Print("Running in series",comm=comm)
    else:
        Print("Running on",size,"processors",comm=comm)
    Print("Petsc arch = %s" % petsc_arch)
    Print("Using scalar type = ",PETSc.ScalarType)

    # Read in options (second arg is the default)
    linop      = opts.getString('linop', False)
    massop     = opts.getString('massop', False)
    linopbeta  = opts.getString('linopbeta', False)
    beta       = opts.getReal('beta',False)
    betas      = opts.getString('betas',False)
    betaRange  = opts.getString('betaRange',False)

    # Options for saving
    outputdir   = opts.getString('outputdir',False)
    flg_leading = opts.getBool('saveLeading',False)
    flg_all     = opts.getBool('saveAll',False)

    if not outputdir:
        outputdir='./'
        Print('No output directory specified. Will save output in this current directory. \n')
    else:
        outputdir += '/'
        Print('Writing to specified output directory %s \n' % outputdir)
        if(rank==0):
            if not os.path.exists(outputdir):
                os.mkdir(outputdir)

    t1 = PETSc.Log.getTime()
    # Read in the file
    Print('### Reading in files and setting options ###')

    if(linop!=False):
        L = PETSc.Mat().load(PETSc.Viewer().createBinary(linop, 'r'))
        Print('Reading in L from %s' % linop)
    else:
        raise RuntimeError('Must set a linear operator')

    # L.setOption(PETSc.Mat.Option.NEW_NONZERO_LOCATIONS,True)
    # L.shift(1e-3)
    # L.shift(-1e-3)
    # L.shift(10)
    # L.shift(-10)
    # L.shift(-1)
    # print(massop)
    if(massop!=False):
        B = PETSc.Mat().load(PETSc.Viewer().createBinary(massop, 'r'))
        Print('Reading in B from %s' % massop)

    # Implement different betas later

    # if(beta!=0 and linopbeta==False):
    #      raise RuntimeError('To use a non-zero spanwise wavenumber you must set the spanwise linear operator using -linopbeta')
    # else:
    #      LB = PETSc.Mat().load(PETSc.Viewer().createBinary(linopbeta, 'r'))
    #      Print('Reading in LB from %s' % linopbeta)

    # Create iterable for beta
    # betaIter = []
    #
    # if(beta!=False):
    #     betaIter.append(beta)
    #
    # if(betas!=False):
    #     betaIter.extend([float(beta) for beta in str.split(betas)])
    #
    # if(betaRange!=False):
    #     betaR = str.split(betaRange)
    #     NR = int(betaR[2])
    #     betaStart = float(betaR[0])
    #     betaEnd   = float(betaR[1])
    #     betaIter.extend(np.linspace(betaStart,betaEnd,NR))
    #
    # if not betaIter:
    #     betaIter.append(0)

    Print('### Files read and options set ### \n')

    E = SLEPc.EPS()
    E.create()
    # E.setOperators(L,B=B)
    if(massop):
        E.setOperators(L,B=B)
    else:
        E.setOperators(L)
    E.setFromOptions()

    E.setWhichEigenpairs(E.Which.TARGET_MAGNITUDE)
    st = E.getST()
    st.setType(st.Type.SINVERT)

    Print()
    Print("******************************")
    Print("*** Eigenvalue setup *********")
    Print("******************************")
    Print()

    eps_type = E.getType()
    Print("Solution method: %s" % eps_type)
    eps_which = E.getWhichEigenpairs()
    Print( "Which eigenvalues to find: %s" % eps_which)

    st_type = st.getType()
    Print("Spectral transformation method: %s" % st_type)

    shift = st.getShift()
    Print( "Shift = : %f + 1i*(%f)" % (np.real(shift),np.imag(shift)))

    target= E.getTarget()
    Print( "Target = : %f + 1i*(%f)" % (np.real(target),np.imag(target)))

    ksp = st.getKSP()
    ksp_type = ksp.getType()
    Print( "Linear systems solved via %s" % ksp_type)
    pc = ksp.getPC()
    pc_type = pc.getType()

    Print( "Preconditioner type %s" % pc_type)

    nev, ncv, mpd = E.getDimensions()
    Print("Number of requested eigenvalues: %d" % nev)
    nev, ncv, mpd = E.getDimensions()
    Print("Max dimension of the subspace (ncv): %d" % ncv)
    tol, maxit = E.getTolerances()
    Print("Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit))

    E.solve()
    Print()
    Print("******************************")
    Print("*** SLEPc Solution Results ***")
    Print("******************************")
    Print()

    its = E.getIterationNumber()
    Print("Number of iterations of the method: %d" % its)

    nconv = E.getConverged()
    Print("Number of converged eigenpairs %d" % nconv)

    if(not path.exists(outputdir+'singularvalues.txt')):
        Print('Creating the file %s ' % outputdir+'eigenvalues.txt')
        data = PETSc.Viewer().createASCII(outputdir+'eigenvalues.txt', 'w')
        data.printfASCII(" evalue     residual norm \n")
    else:
        Print('Appending data to the previous file %s' % outputdir+'eigenvalues.txt')
        data = PETSc.Viewer().create()
        data.setType('ascii')
        data.setFileMode(PETSc.Viewer().Mode.APPEND)
        data.setFileName(outputdir+'eigenvalues.txt')

    if nconv > 0:
        # Create the results vectors
        vr, wr = L.getVecs()
        vi, wi = L.getVecs()

        Print()
        Print("         k            ||Ax-kx||/||kx|| ")
        Print("-------------------- ------------------")
        for i in range(nconv):
            k = E.getEigenpair(i, vr, vi)
            # k = E.getEigenvalue(i)
            error = E.computeError(i,etype=E.ErrorType.RELATIVE)
            if k.imag != 0.0:
                Print(" %9f%+9fi %12g" % (k.real, k.imag, error))
                data.printfASCII("%9f%+9fi %17.14g \n" % (k.real, k.imag, error))
            else:
                Print(" %12f      %12g" % (k.real, error))
                data.printfASCII("%9f %17.14g \n" % (k.real,  error))
            if(flg_leading and i==0):
                PETSc.Viewer().createBinary('%sevector%9f%+9f.dat' % (outputdir,k.real, k.imag),'w')(vr)
                # viewer(vr)
                # viewer = PETSc.Viewer().createBinary('eigavi.dat','w')
                # viewer(vi)
            elif(flg_all):
                PETSc.Viewer().createBinary('%sevector%9f%+9f.dat' % (outputdir,k.real, k.imag),'w')(vr)

        Print()

    Print()
    t1 = PETSc.Log.getTime() -t1
    Print('Time Taken = %f' % t1)
