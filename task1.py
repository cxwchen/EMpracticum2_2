import numpy as np
import matplotlib.pyplot as plt
import scipy.io as spio

#------------------------------------------------------------------------------
#                                FUNCTIONS
#------------------------------------------------------------------------------

def PolarPhasor(phi, tau):
    '''
    1.1 - Polarization Phasor (eq. 2)
    Parameters
    ----------
    phi : angle of orientation in radians [-pi, pi]
    tau : angle of ellipticity in radians [-pi/4, pi/4]

    Returns
    -------
    p : complex polarization ratio (polarization phasor)

    '''
    p = (np.tan(phi)+1j*np.tan(tau))/(1-1j*np.tan(phi)*np.tan(tau))
    return p

def FresnelCoeff(Er_1, Er_2, th_i):
    '''
    1.2 Complex Fresnel reflection coefficient (eq. 14)
    Parameters
    ----------
    Er_1 : relative dielectric permittivity of medium 1
    Er_2 : Er_2 relative dielectric permittivity of medium 2
    th_i : angle of incidence [radians]

    Returns
    -------
    Gamma_h : Fresnel reflection coefficient for horizontal polarization
    Gamma_v : Fresnel reflection coefficient for vertical polarization
    
    '''
    Gamma_h = (np.cos(th_i)-np.sqrt(Er_2/Er_1-np.square(np.sin(th_i)))) / (np.cos(th_i)+np.sqrt(Er_2/Er_1-np.square(np.sin(th_i))))              
    Gamma_v = -(Er_2*np.cos(th_i)-Er_1*np.sqrt(Er_2/Er_1-np.square(np.sin(th_i)))) / (Er_2*np.cos(th_i)+Er_1*np.sqrt(Er_2/Er_1-np.square(np.sin(th_i))))
    return Gamma_h, Gamma_v