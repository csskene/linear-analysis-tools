import scipy.sparse as sp
import numpy as np
def prim2Cons(Base):
    # Builds a matrix G that maps [rho, u, v, w, T] to
    # [rho, rhou, rhov, rhow, rhoe]
    Gamma = 1.4
    R_Gas = 1./Gamma
    NCV = Base.shape[0]

    RUBar = Base[:,4]
    RVBar = Base[:,5]
    REBar = Base[:,6]
    RBar  = Base[:,3]

    UBar = RUBar/RBar
    VBar = RVBar/RBar
    PBar = (Gamma - 1.)*(REBar-0.5*RBar*(UBar**2+VBar**2))

    # Blocks for rho (identity matrix)
    RR = sp.identity(NCV)

    # Blocks for rho*u
    RU = sp.spdiags(UBar,0,NCV,NCV)
    UU = sp.spdiags(RBar,0,NCV,NCV)

    # Blocks for rho*v
    RV = sp.spdiags(VBar,0,NCV,NCV)
    VV = sp.spdiags(RBar,0,NCV,NCV)

    # Blocks for rho*w
    WW = sp.spdiags(RBar,0,NCV,NCV)

    # Blocks for rho*e
    D = 0.5*(UBar**2+VBar**2)+PBar/RBar/(Gamma-1)
    RT = sp.spdiags(D,0.,NCV,NCV)
    UT = sp.spdiags(RBar*UBar, 0, NCV, NCV)
    VT = sp.spdiags(RBar*VBar, 0, NCV, NCV)
    TT = sp.spdiags(RBar*R_Gas/(Gamma-1), 0, NCV, NCV)

    Z0 = 0*sp.identity(NCV)

    G = sp.bmat([[RR,Z0,Z0,Z0,Z0],[RU, UU, Z0, Z0, Z0],[RV, Z0, VV, Z0, Z0 ],[Z0, Z0, Z0, WW, Z0 ],[RT, UT, VT, Z0, TT]])

    return G

def cons2Prim(Base):
    # Builds a matrix G that maps [rho, rhou, rhov, rhow, rhoe] to
    # [rho, u, v, w, T]
    Gamma = 1.4
    R_Gas = 1./Gamma
    NCV = Base.shape[0]

    RUBar = Base[:,4]
    RVBar = Base[:,5]
    REBar = Base[:,6]
    RBar  = Base[:,3]

    UBar = RUBar/RBar
    VBar = RVBar/RBar
    PBar = (Gamma - 1.)*(REBar-0.5*RBar*(UBar**2+VBar**2))

    # Blocks for rho (identity matrix)
    RR = sp.identity(NCV)

    # Blocks for rho*u
    RU = sp.spdiags(-UBar/RBar,0,NCV,NCV)
    UU = sp.spdiags(1/RBar,0,NCV,NCV)

    # Blocks for rho*v
    RV = sp.spdiags(-VBar/RBar,0,NCV,NCV)
    VV = sp.spdiags(1/RBar,0,NCV,NCV)

    # Blocks for rho*w
    WW = sp.spdiags(1/RBar,0,NCV,NCV)

    # Blocks for rho*e
    D = 0.5*(Gamma - 1)/R_Gas/RBar*(UBar**2 + VBar**2) - PBar /(R_Gas*RBar**2)
    RE = sp.spdiags(D,0.,NCV,NCV)
    UE = sp.spdiags(-(Gamma - 1)/R_Gas *UBar/RBar, 0, NCV, NCV)
    VE = sp.spdiags(-(Gamma - 1)/R_Gas *VBar/RBar, 0, NCV, NCV)
    EE = sp.spdiags((Gamma - 1)/R_Gas/RBar, 0, NCV, NCV)

    Z0 = 0*sp.identity(NCV)

    G = sp.bmat([[RR,Z0,Z0,Z0,Z0],[RU, UU, Z0, Z0, Z0],[RV, Z0, VV, Z0, Z0 ],[Z0, Z0, Z0, WW, Z0 ],[RE, UE, VE, Z0, EE]])

    return G

def EnormQuad(Base):
    # Builds a matrix G that maps [rho, rhou, rhov, rhow, rhoe] to
    # [rho, u, v, w, T]
    Gamma = 1.4
    R_Gas = 1./Gamma
    NCV5 = 5*Base.shape[0]

    RUBar = Base[:,4]
    RVBar = Base[:,5]
    REBar = Base[:,6]
    RBar  = Base[:,3]

    UBar = RUBar/RBar
    VBar = RVBar/RBar
    TBar  = (Gamma - 1.0) *(REBar - 0.5*(RUBar**2 + RVBar**2)/RBar)/(R_Gas*RBar) ;
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
