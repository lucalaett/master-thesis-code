
import numpy as np


def residual_helper(G, chi, T, rho, delta, H, a, deltadot, m , theta, k, time):  
    """
    Returns Psi_X*Y + Psi_Y*X according to the formulas in the thesis
    -> see fct: timing_residual for descritpion
    """
    Theta = lambda t: m*t + theta
    #Prefactors:
    preX = 4*np.pi*G*rho*delta/m**3
    preY = 8*np.pi*G*rho*deltadot/(m**2*(12*np.pi*G*rho-k**2))
    #2mX (f=sin), -2mY (f=cos):
    t1, t2, t3, t4 = time+a*chi, time, time+a*chi+T, time+T
    XY = lambda f : f(2*Theta(t1)) - f(2*Theta(t2)) - f(2*Theta(t3)) + f(2*Theta(t4))
    #X and Y:
    X = XY(np.sin)/(2*m)
    Y = -XY(np.cos)/(2*m)
    return preX*Y + preY*X

def timing_residual(G, chi, T, rho, delta, H, a, deltadot, m , theta, k, tend):  
    """
    Returns the timing residual according to the formulas in the thesis.
    Standard errors are not included in the parameters or the returns. 

    Parameters [all float, all in GeV^x]
    ----------
    G : gravitational constant 
    chi : comoving earth-pulsar distnace
    T : pulsar period
    rho : NSN energy density
    delta : NSN overdensity
    H : cosmological const
    a : expansion factor
    deltadot : (physical) time derivative of delta
    m : mass
    theta : phase
    k : scale 
    tend : timespan of measurement 

    Returns
    -------
    R : float, timing residual 
    """
    tend_part = residual_helper(G, chi, T, rho, delta, H, a, deltadot, m , theta, k, tend)/T
    t0_part = residual_helper(G, chi, T, rho, delta, H, a, deltadot, m , theta, k, 0)/T
    return tend_part - t0_part




