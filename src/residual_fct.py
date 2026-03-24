
import numpy as np


def residual_XY(tobs, f, T, deltaE, deltaP, a, m, theta, chi):  
    """
    Returns X or -Y according to the formulas in the thesis

    Parameters 
    ----------
    tobs : float, observation time [GeV^-1]
    f : fct, either np.sin-> returns X  or  np.cos-> returns -Y 
    -> see fct: timing_residual for descritpion of the rest
    """
    Theta = lambda t: m*t + theta
    t1, t2, t3, t4 = tobs, tobs+T, tobs-a*chi, tobs-a*chi+T
    XY = deltaE * (f(2*Theta(t1)) - f(2*Theta(t2))) - deltaP * (f(2*Theta(t3)) - f(2*Theta(t4)))
    XY /= 2*m    
    return XY

def timing_residual(G, T, rho, deltaE, deltaP, H, a, tstart, tend, m , k, theta, chi):  
    """
    Returns the timing residual according to the formulas in the thesis.
    Standard errors are not included in the parameters or the returns. 

    Parameters [all float, all in GeV^x]
    ----------
    G : gravitational constant 
    T : pulsar period
    rho : NSN energy density
    deltaE : NSN overdensity close to earth
    deltaE : NSN overdensity close to pulsar
    H : cosmological const
    a : expansion factor
    tstart : measurement startingtime
    tend : measurement endtime
    m : mass
    theta : phase
    k : scale 
    chi : comoving earth-pulsar distnace

    Returns
    -------
    R : float, timing residual [GeV^-1]
    """
    #prefactors
    PsiX = 8*np.pi*G*rho/m**2
    PsiY = H/m * 8*np.pi**(0.5)*G*rho/(12*np.pi*G*rho-k**2) 
    #X and Y
    X = lambda tobs: residual_XY(tobs, np.sin, T, deltaE, deltaP, a, m, theta, chi)
    Y = lambda tobs: -residual_XY(tobs, np.cos, T, deltaE, deltaP, a, m, theta, chi)
    #residual
    R = PsiX/(2*m*T) * (Y(tend)-Y(tstart)) + PsiY/(2*m*T) * (X(tend)-X(tstart))
    return R #GeV^-1




