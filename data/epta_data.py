

def p_name(): 
    """
    Pulsar Jnames of pulsars measured in EPTA DR2
    
    Returns 
    ----------
    p_name: array, Pulsar Jnames
    """
    return np.array(['J0030+0451', 'J0613−0200', 'J0751+1807', 'J0900−3144',
                     'J1012+5307', 'J1022+1001', 'J1024−0719', 'J1455−3330',
                     'J1600−3053', 'J1640+2224', 'J1713+0747', 'J1730−2304',
                     'J1738+0333', 'J1744−1134', 'J1751−2857', 'J1801−1417',
                     'J1804−2717', 'J1843−1113', 'J1857+0943', 'J1909−3744',
                     'J1910+1256', 'J1911+1347', 'J1918−0642', 'J2124−3358',
                     'J2322+2057'])

def p_period():
    """
    Pulsar emission periods of pulsars measured in EPTA DR2. 
    The period is based on the spin-frequency in tables B1-B7 of paper I, standard errors not included
    
    Returns 
    ----------
    T: array, emission period of pulsar [GeV^-1]
    """
    nu = np.array([205.530695938456, 326.6005620234831, 287.457853995106, 90.011841919354, 
               190.2678344415654, 60.7794479566973, 193.715683448548, 125.200243244993, 
               277.9377069896062, 316.123979331869, 218.8118404171605, 123.110287147370,
               170.937369887100, 245.4261196898081, 255.43611088568, 275.85470899694,
               107.031649219533, 541.809745036152, 186.4940783779890, 339.3156872184705,
               200.658802230113, 216.171227371979, 130.789514123371, 202.793893746013,
               207.96816335836]) #Hz
    return 1.52e24/nu #GeV^-1


def p_tend():
    """
    Measurement timespan of each Pulsar in EPTA DR2. 
    The value is copied from table 2 of paper I, standard errors not included
    
    Returns 
    ----------
    tend: array, measurement timespan of pulsar [GeV^-1]
    """
    tend = np.array([22.0, 22.9, 24.2, 13.6, 23.7, 24.5, 23.1, 15.7, 14.3, 24.4, 24.5, 16.1,
                 14.1, 24.0, 14.7, 13.7, 14.7, 16.8, 24.1, 15.7, 15.2, 14.2, 19.7, 16.0, 14.7]) #yr
    return tend*4.8e31

def p_wrms():
    """
    RMS (whitened) of timing residuals of each Pulsar in EPTA DR2. 
    The value is copied from table 2 of paper I, standard errors not included
    
    Returns 
    ----------
    wmrs: array, rms(whitened) of timing residual of pulsar [GeV^-1]
    """
    wrms = np.array([2.30, 1.06, 1.50, 2.60, 1.02, 1.56, 1.02, 2.46, 0.37, 1.10, 0.20, 0.83,
                 2.33, 0.56, 2.34, 2.46, 1.63, 0.81, 1.05, 0.14, 1.77, 0.75, 1.31, 2.17, 4.08]) #microsec
    return wrms*1.52e18 #GeV^-1
    








