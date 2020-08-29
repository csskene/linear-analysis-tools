import sys,os,datetime
import numpy as np

import os.path
from os import path

import petsc4py
import slepc4py

petsc_arch = os.environ["PETSC_ARCH"]
slepc4py.init(sys.argv,arch=petsc_arch)

from petsc4py import PETSc
from slepc4py import SLEPc

Print = PETSc.Sys.Print

class resolvent(object):
    def __init__(self, L,WO,WI,omega,alpha,LB=False,beta=False,SP=False):
        # Apply the weight matrices to get the weighted inverse resolvent

        Print( "***********************" )
        Print( "*** Resolvent setup ***" )
        Print( "***********************\n" )
        # Create the inverse resolvent
        self.RI = L.duplicate(copy=True)
        if LB:
            self.LBR = LB.duplicate(copy=True)
            self.LBR.realPart()
            self.LBI = LB.duplicate(copy=True)
            self.LBI.imagPart()
            self.RI += beta**2*self.LBR+1j*beta*self.LBI
        if SP:
            self.RI -= SP

        # Allow to create new non-zeros on the diagonal
        self.RI.setOption(PETSc.Mat.Option.NEW_NONZERO_LOCATIONS,True)
        self.RI.shift(1j*(omega+1j*alpha))
        self.RI.scale(-1)

        self.R = PETSc.KSP().create()
        self.R.setOperators(self.RI)
        # Default Solve linear systems with LU decomposition
        self.R.setType('preonly')
        pc = self.R.getPC()
        pc.setType('lu')
        # Can overwrite from options
        self.R.setFromOptions()
        pc.setFromOptions()
        # Setup now so can print useful information before any solves
        self.R.setUp()
        ksp_type = self.R.getType()
        Print( "Linear system method: %s" % ksp_type )
        pc_type = pc.getType()
        Print( "Preconditioner: %s" % pc_type )
        pc_solver_type = pc.getFactorSolverType()
        Print( "Factor solver type: %s" % pc_solver_type )
        Print()
        self.WO = WO
        self.WI = WI

    def changeRes(self,L,omega,alpha,beta=False):
        L.copy(self.RI)
        if beta:
            self.RI += beta**2*self.LBR+1j*beta*self.LBI
        self.RI.shift(1j*(omega+1j*alpha))
        self.RI.scale(-1)

    def mult(self, A, x, y):
        f, q = self.WO.getVecRight(), self.WO.getVecRight()
        self.WI.mult(x, f)
        self.R.solve(f, q)
        self.WO.mult(q, y)

    def multHermitian(self, A, x, y):
        f, q = self.WI.getVecRight(), self.WI.getVecRight()
        self.WO.multTranspose(x, f)
        f.conjugate()
        self.R.solveTranspose(f, q)
        q.conjugate()
        self.WI.multTranspose(q, y)

class randomisedResolvent(object):
    def __init__(self, gains,leftVecs,rightVecs,WR):
        # Class to hold randomised resolvent solution so can use the same commands
        # as for the traditional SVD
        self.gains     = gains
        self.leftVecs  = leftVecs
        self.rightVecs = rightVecs
        self.WR = WR

    def getSingularTriplet(self,i, u, v):
        self.leftVecs.getColumnVector(i, u)
        self.rightVecs.getColumnVector(i, v)
        return self.gains[i]

    def computeError(self,i):
        u, v = self.WR.getVecs()
        self.leftVecs.getColumnVector(i, u)
        self.rightVecs.getColumnVector(i, v)
        Av, AHu = self.WR.getVecs()
        Av = self.WR*v
        self.WR.multHermitian(u,AHu)
        u1 = Av-self.gains[i]*u
        v1 = AHu-self.gains[i]*v

        norm1 = u1.norm()
        norm2 = v1.norm()

        return np.sqrt(norm1**2+norm2**2)/self.gains[i]

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
    linopbeta  = opts.getString('linopbeta', False)
    wifile     = opts.getString('wifile', False)
    wofile     = opts.getString('wofile',False)
    testadj    = opts.getBool('test_adj', False)
    omega      = opts.getReal('omega',False)
    alpha      = opts.getReal('disc',False)
    beta       = opts.getReal('beta',False)
    omegas     = opts.getString('omegas',False)
    betas      = opts.getString('betas',False)
    betaRange  = opts.getString('betaRange',False)
    omegaRange = opts.getString('omegaRange',False)
    alphas     = opts.getString('discs',False)
    alphaRange = opts.getString('discRange',False)
    sponge     = opts.getString('sponge',False)
    physvec    = opts.getString('physvec',False)

    # Options for saving
    outputdir   = opts.getString('outputdir',False)
    flg_leading = opts.getBool('saveLeading',False)
    flg_all     = opts.getBool('saveAll',False)

    # Use randomised resolvent
    flgRandRes = opts.getBool('randomisedResolvent',False)
    k          = opts.getInt('randk',10)

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
    if(wofile!=False):
        WO = PETSc.Mat().load(PETSc.Viewer().createBinary(wofile, 'r'))
        Print('Reading in WO from %s' % wofile)
    else:
        Print('No WO matrix set')

    if(wifile!=False):
        WI = PETSc.Mat().load(PETSc.Viewer().createBinary(wifile, 'r'))
        Print('Reading in WI from %s' % wifile)
    else:
        Print('No WI matrix set')

    if(sponge!=False):
        SP = PETSc.Mat().load(PETSc.Viewer().createBinary(sponge, 'r'))
        Print('Reading in SP from %s' % sponge)
    else:
        Print('No sponge set')

    if(physvec!=False):
        phi = PETSc.Mat().load(PETSc.Viewer().createBinary(physvec, 'r'))
        Print('Reading in Physics-Based Vector from %s' % physvec)
    else:
        Print('No sacling in random test vector')


    if(linop!=False):
        L = PETSc.Mat().load(PETSc.Viewer().createBinary(linop, 'r'))
        Print('Reading in L from %s' % linop)
    else:
        raise RuntimeError('Must set a linear operator')

    # Create an iterable for omega
    omegaIter = []
    if(omega!=False):
        omegaIter.append(omega)

    if(omegas!=False):
        omegaIter.extend([float(om) for om in str.split(omegas)])

    if(omegaRange!=False):
        omegaR = str.split(omegaRange)
        NR = int(omegaR[2])
        omegaStart = float(omegaR[0])
        omegaEnd   = float(omegaR[1])
        omegaIter.extend(np.linspace(omegaStart,omegaEnd,NR))
    if not omegaIter:
        raise RuntimeError('Must set omegas to calculate the resolvent modes for')

    # Create iterable for alpha
    alphaIter = []

    if(alpha!=False):
        alphaIter.append(alpha)

    if(alphas!=False):
        alphaIter.extend([float(alph) for alph in str.split(alphas)])

    if(alphaRange!=False):
        alphaR = str.split(alphaRange)
        NR = int(alphaR[2])
        alphaStart = float(alphaR[0])
        alphaEnd   = float(alphaR[1])
        alphaIter.extend(np.linspace(alphaStart,alphaEnd,NR))
    if not alphaIter:
        alphaIter.append(0)

    if((beta!=0 or betas!=0 or betaRange!=0) and linopbeta==False):
         raise RuntimeError('To use a non-zero spanwise wavenumber you must set the spanwise linear operator using -linopbeta')
    elif((beta!=0 or betas!=0 or betaRange!=0)):
         LB = PETSc.Mat().load(PETSc.Viewer().createBinary(linopbeta, 'r'))
         Print('Reading in LB from %s' % linopbeta)

    # Create iterable for beta
    betaIter = []

    if(beta!=False):
        betaIter.append(beta)

    if(betas!=False):
        betaIter.extend([float(beta) for beta in str.split(betas)])

    if(betaRange!=False):
        betaR = str.split(betaRange)
        NR = int(betaR[2])
        betaStart = float(betaR[0])
        betaEnd   = float(betaR[1])
        betaIter.extend(np.linspace(betaStart,betaEnd,NR))

    if not betaIter:
        betaIter.append(0)

    Print('### Files read and options set ### \n')

    params = [(om,alph,bet) for om in omegaIter for alph in alphaIter for bet in betaIter]
    # To more easily identify cases in the output
    for iter,param in enumerate(params):
        omega = param[0]
        alpha = param[1]
        beta  = param[2]
        # Set up the resolvent linear systems solver
        # If WO or WI are provided this resolvent is weighted
        if(iter==0):
            if linopbeta:
                resCtx = resolvent(L,WO,WI,omega,alpha,LB=LB,beta=beta)
            else:
                resCtx = resolvent(L,WO,WI,omega,alpha)
            n, m = L.getSize()
            WR = PETSc.Mat().createPython((n, n), resCtx)
            WR.setUp()

            if(testadj):
                Print('Test adj')
                x1, x2 = WR.getVecs()
                x1.setRandom()
                x2.setRandom()
                Ax1, AHx2 = WR.getVecs()
                WR.mult(x1,Ax1)
                WR.multHermitian(x2,AHx2)
                norm1 = x2.dot(Ax1)
                norm2 = AHx2.dot(x1)
                err = np.abs(norm1-norm2)/np.abs(norm1)
                if(err < 1e-12):
                    Print('Adjoint test passed with err = %g' % err)
                else:
                    Print('Error is %g Check the code' % err)
        else:
            if linopbeta:
                resCtx.changeRes(L,omega,alpha,beta=beta)
            else:
                resCtx.changeRes(L,omega,alpha)

        if(not flgRandRes):
            S = SLEPc.SVD()
            S.create()
            S.setOperator(WR)
            S.setFromOptions()
            S.setUp()
            if(iter==0):
                Print( "*****************" )
                Print( "*** SVD setup ***" )
                Print( "*****************\n" )

                svd_type = S.getType()
                tol, maxit = S.getTolerances()
                nsv, ncv, mpd = S.getDimensions()

                Print( "Solution method: %s" % svd_type )
                Print( "Stopping condition: tol=%.4g, maxit=%d" % (tol, maxit) )
                Print( "Number of requested singular values: %d" % nsv )
                Print( "Max dimension of the subspace (ncv): %d" % ncv )

                Print()
                Print( "***********************" )
                Print( "*** Running the SVD ***" )
                Print( "***********************\n" )

            S.solve()

            nconv = S.getConverged()
        else:
            Print('Running the SVD using the randomised resolvent method')
            Print('k = ',k)
            test = WR.getVecRight()
            sketch = WR.getVecRight()

            rand = PETSc.Random().create()
            rand.setType(rand.Type.RAND)

            bv = SLEPc.BV()
            bv.create()
            bv.setType('vecs')
            bv.setSizesFromVec(test, k)
            for i in range(k):
                test.setRandom(rand)
                # added conditional for physics random test vector
                if(physvec!=False):
                    phi.mult(test,omega)
                else:
                    omega = test
                WR.mult(omega, sketch)
                bv.insertVec(i,sketch)

            bv.orthogonalize()

            n = test.getSize()
            B = PETSc.Mat().createDense((n,k))
            B.setUp()

            y = WR.getVecRight()
            for i in range(k):
                v = bv.getColumn(i)
                WR.multHermitian(v, y)
                bv.restoreColumn(i,v)
                low,high=y.getOwnershipRange()
                B.setValues(np.arange(low,high,dtype='int32'),i,y)
            bv.destroy()
            B.assemblyBegin(B.AssemblyType.FINAL)
            B.assemblyEnd(B.AssemblyType.FINAL)

            S = SLEPc.SVD()
            S.create()
            S.setOperator(B)
            S.solve()

            nconv = S.getConverged()
            v, u = B.getVecs()

            US = PETSc.Mat().createDense((n,nconv))
            US.setUp()

            V = PETSc.Mat().createDense((n,nconv))
            V.setUp()

            for i in range(nconv):
                S.getSingularTriplet(i, u, v)
                low,high=u.getOwnershipRange()
                V.setValues(np.arange(low,high,dtype='int32'),i,u)
                y = WR*u
                low,high=y.getOwnershipRange()
                US.setValues(np.arange(low,high,dtype='int32'),i,y)

            V.assemblyBegin(V.AssemblyType.FINAL)
            V.assemblyEnd(V.AssemblyType.FINAL)

            US.assemblyBegin(US.AssemblyType.FINAL)
            US.assemblyEnd(US.AssemblyType.FINAL)

            S.destroy()
            S = SLEPc.SVD()
            S.create()

            S.setOperator(US)
            S.solve()
            nconv = S.getConverged()

            Vbar = PETSc.Mat().createDense((nconv,nconv))
            Vbar.setUp()


            leftVecs = PETSc.Mat().createDense((n,nconv))
            leftVecs.setUp()

            gains = []
            for i in range(nconv):
                sigma = S.getSingularTriplet(i, u, v)

                low,high=v.getOwnershipRange()
                Vbar.setValues(i,np.arange(low,high,dtype='int32'),v)
                gains.append(sigma)
                low,high=u.getOwnershipRange()
                leftVecs.setValues(np.arange(low,high,dtype='int32'),i,u)

            leftVecs.assemblyBegin(leftVecs.AssemblyType.FINAL)
            leftVecs.assemblyEnd(leftVecs.AssemblyType.FINAL)

            Vbar.assemblyBegin(Vbar.AssemblyType.FINAL)
            Vbar.assemblyEnd(Vbar.AssemblyType.FINAL)

            rightVecs = V.matMult(Vbar)
            S.destroy()
            S = randomisedResolvent(gains,leftVecs,rightVecs,WR)

            # for i in range(nconv):
            #     rightVecs.getColumnVector(i, u)
            #     PETSc.Viewer().createBinary('%sMf_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(u)

        if nconv > 0:
            v, u = WR.getVecs()
            Miv, Miu = WR.getVecs()
            q, f = WO.getVecRight(), WO.getVecRight()
            if(iter==0):
                if(not path.exists(outputdir+'singularvalues.txt')):
                    Print('Creating the file %s ' % outputdir+'singularvalues.txt')
                    data = PETSc.Viewer().createASCII(outputdir+'singularvalues.txt', 'w')
                    data.printfASCII("     Omega        alpha               beta         i         sigma          dsigma/domega      dsigma/dalpha       dsigma/dbeta       residual norm \n")
                else:
                    Print('Appending data to the previous file %s' % outputdir+'singularvalues.txt')
                    data = PETSc.Viewer().create()
                    data.setType('ascii')
                    data.setFileMode(PETSc.Viewer().Mode.APPEND)
                    data.setFileName(outputdir+'singularvalues.txt')
                    # The next line would be better but doesn't seem to work
                    # data = PETSc.Viewer().createASCII(outputdir+'singularvalues.txt',mode=PETSc.Viewer().Mode.APPEND)
                Print()
                Print("     Omega               alpha               beta         i         sigma          dsigma/domega      dsigma/dalpha       dsigma/dbeta       residual norm ")
                Print("-----------------  -----------------  -----------------  ---  -----------------  -----------------  -----------------  -----------------  -------------------")
            for i in range(nconv):
                sigma = S.getSingularTriplet(i, u, v)
                error = S.computeError(i)

                # Compute the sensitivies
                sens  = -sigma**2*u.dot(v)
                if linopbeta:
                    DA = 2.*beta*resCtx.LBR+1j*resCtx.LBI
                    Miu = WO*DA*WI*u
                    sensB = sigma**2*Miu.dot(v)
                else:
                    sensB = 0.
                Print( "%17.14g  %17.14g  %17.14g  %03d  %17.14g  %17.14g  %17.14g  %17.14g  %17.14g" % (omega,alpha,beta,i, sigma,np.imag(sens),np.real(sens),np.real(sensB), error))
                data.printfASCII("%17.14g  %17.14g  %17.14g  %03d  %17.14g  %17.14g  %17.14g  %17.14g  %17.14g \n" % (omega,alpha,beta,i, sigma,np.imag(sens),np.real(sens),np.real(sensB), error))

                if(flg_leading and i==0):
                    PETSc.Viewer().createBinary('%s_Mf_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(v)
                    PETSc.Viewer().createBinary('%s_Mq_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(u)
                    # Miv = WI*v
                    # Miu = WI*u
                    # PETSc.Viewer().createBinary('%sf_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(Miv)
                    # PETSc.Viewer().createBinary('%sq_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(Miu)
                elif(flg_all):
                    PETSc.Viewer().createBinary('%s_Mf_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(v)
                    PETSc.Viewer().createBinary('%s_Mq_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(u)
                    # Miv = WI*v
                    # Miu = WI*u
                    # PETSc.Viewer().createBinary('%sf_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(Miv)
                    # PETSc.Viewer().createBinary('%sq_k%03d_om%05.2f_alpha%05.2f_beta%05.2f.dat' % (outputdir, i, omega, alpha, beta), 'w')(Miu)
        if(not flgRandRes):
            S.destroy()


    Print()
    t1 = PETSc.Log.getTime() -t1
    Print('Time Taken = %f' % t1)
