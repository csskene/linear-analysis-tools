import scipy.sparse as sp
import numpy as np
from scipy.interpolate import griddata,interp2d
def prim2Cons(Base):
    # Builds a matrix G that maps [rho, u, v, w, T] to
    # [rho, rhou, rhov, rhow, rhoe]
    Gamma = 1.4
    R_Gas = 1./Gamma
    NCV = Base.shape[0]

    RUBar = Base[:,4]
    RVBar = Base[:,5]
    RWBar = Base[:,6]
    REBar = Base[:,7]
    RBar  = Base[:,3]

    UBar = RUBar/RBar
    VBar = RVBar/RBar
    WBar = RWBar/RBar
    PBar = (Gamma - 1.)*(REBar-0.5*RBar*(UBar**2+VBar**2+WBar**2))

    # Blocks for rho (identity matrix)
    RR = sp.identity(NCV)

    # Blocks for rho*u
    RU = sp.spdiags(UBar,0,NCV,NCV)
    UU = sp.spdiags(RBar,0,NCV,NCV)

    # Blocks for rho*v
    RV = sp.spdiags(VBar,0,NCV,NCV)
    VV = sp.spdiags(RBar,0,NCV,NCV)

    # Blocks for rho*w
    RW = sp.spdiags(WBar,0,NCV,NCV)
    WW = sp.spdiags(RBar,0,NCV,NCV)

    # Blocks for rho*e
    D = 0.5*(UBar**2+VBar**2+WBar**2)+PBar/RBar/(Gamma-1)
    RT = sp.spdiags(D,0.,NCV,NCV)
    UT = sp.spdiags(RBar*UBar, 0, NCV, NCV)
    VT = sp.spdiags(RBar*VBar, 0, NCV, NCV)
    WT = sp.spdiags(RBar*WBar, 0, NCV, NCV)
    TT = sp.spdiags(RBar*R_Gas/(Gamma-1), 0, NCV, NCV)

    Z0 = 0*sp.identity(NCV)

    G = sp.bmat([[RR,Z0,Z0,Z0,Z0],[RU, UU, Z0, Z0, Z0],[RV, Z0, VV, Z0, Z0 ],[RW, Z0, Z0, WW, Z0 ],[RT, UT, VT, WT, TT]])

    return G

def cons2Prim(Base):
    # Builds a matrix G that maps [rho, rhou, rhov, rhow, rhoe] to
    # [rho, u, v, w, T]
    Gamma = 1.4
    R_Gas = 1./Gamma
    NCV = Base.shape[0]

    RUBar = Base[:,4]
    RVBar = Base[:,5]
    RWBar = Base[:,6]
    REBar = Base[:,7]
    RBar  = Base[:,3]

    UBar = RUBar/RBar
    VBar = RVBar/RBar
    WBar = RWBar/RBar
    PBar = (Gamma - 1.)*(REBar-0.5*RBar*(UBar**2+VBar**2+WBar**2))

    # Blocks for rho (identity matrix)
    RR = sp.identity(NCV)

    # Blocks for rho*u
    RU = sp.spdiags(-UBar/RBar,0,NCV,NCV)
    UU = sp.spdiags(1/RBar,0,NCV,NCV)

    # Blocks for rho*v
    RV = sp.spdiags(-VBar/RBar,0,NCV,NCV)
    VV = sp.spdiags(1/RBar,0,NCV,NCV)

    # Blocks for rho*w
    RW = sp.spdiags(-WBar/RBar,0,NCV,NCV)
    WW = sp.spdiags(1/RBar,0,NCV,NCV)

    # Blocks for rho*e
    D = 0.5*(Gamma - 1)/R_Gas/RBar*(UBar**2 + VBar**2 + WBar**2) - PBar /(R_Gas*RBar**2)
    RE = sp.spdiags(D,0.,NCV,NCV)
    UE = sp.spdiags(-(Gamma - 1)/R_Gas *UBar/RBar, 0, NCV, NCV)
    VE = sp.spdiags(-(Gamma - 1)/R_Gas *VBar/RBar, 0, NCV, NCV)
    WE = sp.spdiags(-(Gamma - 1)/R_Gas *WBar/RBar, 0, NCV, NCV)
    EE = sp.spdiags((Gamma - 1)/R_Gas/RBar, 0, NCV, NCV)

    Z0 = 0*sp.identity(NCV)

    G = sp.bmat([[RR,Z0,Z0,Z0,Z0],[RU, UU, Z0, Z0, Z0],[RV, Z0, VV, Z0, Z0 ],[RW, Z0, Z0, WW, Z0 ],[RE, UE, VE, WE, EE]])

    return G

def EnormQuad(Base):
    # Builds a matrix G that maps [rho, rhou, rhov, rhow, rhoe] to
    # [rho, u, v, w, T]
    Gamma = 1.4
    R_Gas = 1./Gamma
    NCV5 = 5*Base.shape[0]

    RUBar = Base[:,4]
    RVBar = Base[:,5]
    RWBar = Base[:,6]
    REBar = Base[:,7]
    RBar  = Base[:,3]

    UBar = RUBar/RBar
    VBar = RVBar/RBar
    WBar = RWBar/RBar
    TBar  = (Gamma - 1.0) *(REBar - 0.5*(RUBar**2 + RVBar**2 + RWBar**2)/RBar)/(R_Gas*RBar) ;
    Volm = Base[:,2]
    Volm = Volm/np.linalg.norm(Volm)

    # Compute diagonal elements
    DR = Volm*R_Gas*TBar/RBar
    DU = Volm*RBar
    DT = Volm*R_Gas/(Gamma-1.0)*RBar/TBar


    DAll = np.hstack((DR,DU,DU,DU,DT))

    # Input
    FI = sp.spdiags(1/np.sqrt(DAll),0,NCV5,NCV5)
    # Output
    FO = sp.spdiags(np.sqrt(DAll),0,NCV5,NCV5)

    return FI,FO

def linopSponge(Base,sponge):

    NCV  = Base.shape[0]
    NCV5 = 5*NCV

	#Baseflow
    X     = Base[:, 0]
    Y     = Base[:, 1]
    Volm  = Base[:, 2]

	# Max grid values
    Xmax           = np.amax(X)
    Xmin           = np.amin(X)
    Ymax           = np.amax(Y)
    Ymin           = np.amin(Y)

	# Sponge size limits
    sponge_Xmax  = sponge[0]
    sponge_Xmin  = sponge[1]
    sponge_Ymax  = sponge[2]
    sponge_Ymin  = sponge[3]

    # sponge parameters
    sponge_strenth = 2.0
    sponge_a       = 0.068
    sponge_b       = 0.845
    sponge_n       = 2.0
    sponge_m       = 8.0

	# Sponge matrix
    SP = np.zeros((NCV5))   # Sponge coefficients

    for icv in range(0, NCV):
        xloc = X[icv]
        yloc = Y[icv]
        #vloc = Volm[icv]
        vloc = 1.0

        if ( yloc >= sponge_Ymax ):
            sponge_xp    =  ( yloc - sponge_Ymax ) / np.absolute(Ymax - sponge_Ymax)
            sponge_coeff = 2.0*sponge_strenth*( sponge_a *(sponge_xp**sponge_n) + sponge_b *(sponge_xp**sponge_m) )
            #sponge_coeff = ( yloc - sponge_Ymax )**2
            SP[icv] = SP[icv] + vloc*sponge_coeff

        if ( yloc <= sponge_Ymin ):
            sponge_xp    = -( yloc - sponge_Ymin ) / np.absolute(Ymin - sponge_Ymin)
            sponge_coeff = 2.0*sponge_strenth*( sponge_a *(sponge_xp**sponge_n) + sponge_b *(sponge_xp**sponge_m) )
            #sponge_coeff = ( yloc - sponge_Ymin )**2
            SP[icv] = SP[icv] + vloc*sponge_coeff

        if ( xloc >= sponge_Xmax ):
            sponge_xp    =  ( xloc - sponge_Xmax ) / np.absolute(Xmax - sponge_Xmax)
            sponge_coeff = 1.5*sponge_strenth*( sponge_a *(sponge_xp**sponge_n) + sponge_b *(sponge_xp**sponge_m) )
            #sponge_coeff = ( xloc - sponge_Xmax )**2
            SP[icv] = SP[icv] + vloc*sponge_coeff

        if ( xloc <= sponge_Xmin ):
            sponge_xp    =  -( xloc - sponge_Xmin ) / np.absolute(Xmin - sponge_Xmin)
            sponge_coeff = 2.0*sponge_strenth*( sponge_a *(sponge_xp**sponge_n) + sponge_b *(sponge_xp**sponge_m) )
            #sponge_coeff = ( xloc - sponge_Xmin )**2
            SP[icv] = SP[icv] + vloc*sponge_coeff

    SP = np.concatenate((SP,SP,SP,SP,SP), axis=None)
    SP = sp.spdiags(SP,0,NCV5,NCV5)

    return SP


def physVec(Base):

    Gamma = 1.4
    R_Gas = 1./Gamma
    NCV   = Base.shape[0]
    NCV5  = 5*NCV

    X     = Base[:,0]
    Y     = Base[:,1]
    Volm  = Base[:,2]
    RBar  = Base[:,3]
    RUBar = Base[:,4]
    RVBar = Base[:,5]
    RWBar = Base[:,6]
    REBar = Base[:,7]

    # Max grid values
    xmax  = np.amax(X)
    xmin  = np.amin(X)
    ymax  = np.amax(Y)
    ymin  = np.amin(Y)

    # Convert to primitive variables
    UBar = np.divide(RUBar,RBar)
    VBar = np.divide(RVBar,RBar)
    WBar = np.divide(RWBar,RBar)
    PBar = (Gamma - 1.)*(REBar-0.5*RBar*(UBar**2+VBar**2+WBar**2))

    # Generate structured mesh for interpolation
    dx = 0.1
    dy = 0.1
    xm = np.arange(xmin,xmax,dx)
    ym = np.arange(ymin,ymax,dy)
    XG, YG = np.meshgrid(xm, ym, sparse=False, indexing='ij')

    # Interpolation onto structured grid
    Rbar_Interp = griddata((X,Y), RBar, (XG, YG), method='cubic') ;
    Ubar_Interp = griddata((X,Y), UBar, (XG, YG), method='cubic') ;
    Vbar_Interp = griddata((X,Y), VBar, (XG, YG), method='cubic') ;
    Wbar_Interp = griddata((X,Y), WBar, (XG, YG), method='cubic') ;
    Pbar_Interp = griddata((X,Y), PBar, (XG, YG), method='cubic') ;

    # Computing gradients in structured grid
    drdx,drdy = np.gradient(Rbar_Interp,dx,dy)
    dudx,dudy = np.gradient(Ubar_Interp,dx,dy)
    dvdx,dvdy = np.gradient(Vbar_Interp,dx,dy)
    dwdx,dwdy = np.gradient(Wbar_Interp,dx,dy)
    dpdx,dpdy = np.gradient(Pbar_Interp,dx,dy)

    dr = np.sqrt(np.square(drdx) + np.square(drdy))
    du = np.sqrt(np.square(dudx) + np.square(dudy))
    dv = np.sqrt(np.square(dvdx) + np.square(dvdy))
    dw = np.sqrt(np.square(dwdx) + np.square(dwdy))
    dp = np.sqrt(np.square(dpdx) + np.square(dpdy))


    # Interpolate back to unstructured grid
    f = interp2d(xm, ym, np.transpose(dr))
    DR=[]
    for i, j in zip(X,Y):
        DR.append(f(i,j))
    DU=[]
    for i, j in zip(X,Y):
        DU.append(f(i,j))
    DV=[]
    for i, j in zip(X,Y):
        DV.append(f(i,j))
    DW=[]
    for i, j in zip(X,Y):
        DW.append(f(i,j))
    DP=[]
    for i, j in zip(X,Y):
        DP.append(f(i,j))

    # Stack gradients per variable in sparse diagonal matrix
    phi = np.concatenate((DR,DU,DV,DW,DP), axis=None)
    phi = sp.spdiags(phi,0,NCV5,NCV5)

    return(phi)
