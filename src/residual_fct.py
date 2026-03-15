


def deltaT(G, chi, T, rho, delta, H, a, deltadot, m , theta, k, time):  
    """
    Returns deltaT according to the formulas in the thesis.
    Standard errors are not included in the parameters or the returns. 

    Parameters [all float, all in GeV^x]
    ----------
    -> see fct: timing_residual for descritpion

    Returns
    -------
    deltaT : float, delay in pulsar period 
    """
    #Theta
    Theta = lambda t: m*t + theta
    #Prefactor: Psi_X 
    preX = 4*np.pi*G*rho*delta/m**3
    #Prefactor: Psi_Y
    preY = 8*np.pi*G*rho*deltadot/(m**2*(12*np.pi*G*rho-k**2))
    #Calculate 2mX (f=sin) and -2mY (f=cos)
    t1, t2, t3, t4 = time+a*chi, time, time+a*chi+T, time+T
    XY = lambda f : f(2*Theta(t1)) - f(2*Theta(t2)) - f(2*Theta(t3)) + f(2*Theta(t4))
    #Calculate deltaT
    return preX * XY(np.sin) + preY * XY(np.cos)




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
    innerfct = lambda t: deltaT(G, chi, T, rho, delta, H, a, deltadot, m , theta, k, t)/T
    #integration over full period = 0 -> neglegtion before integration is less expensive
    tend = tend % (np.pi/m) 
    return quad(innerfct, 0.0, tend)[0]




