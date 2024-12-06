# -*- coding: utf-8 -*-
"""
Created on Wed Nov 27 12:10:15 2024

Script illustrating the polarization pattern shown in the EE2P1
lab seesion 2

@author: PJA
"""
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
    Gamma_h = (np.cos(th_i)-np.sqrt(Er_2/Er_1-np.square(np.sin(th_i))))/ \
        (np.cos(th_i)+np.sqrt(Er_2/Er_1-np.square(np.sin(th_i))))              
    Gamma_v = -(Er_2*np.cos(th_i)-Er_1*np.sqrt(Er_2/Er_1-np.square(np.sin(th_i))))/ \
        (Er_2*np.cos(th_i)+Er_1*np.sqrt(Er_2/Er_1-np.square(np.sin(th_i))))
    return Gamma_h, Gamma_v


def MultipathRxPolarState(pt, f, R, H, Er):
    '''
    1.3 Polarization state of the received wave in multipath propagation  (eq. 13)
    Parameters
    ----------
    pt : transmit wave polarization phasor
    f : frequency [Hz]
    R: distance between antennas [m]
    H: reflector height [m]
    Er: reflector material relative permittivity
    
    Returns
    -------
    pr: received wave polarization phasor
    
    '''
    k = 2*np.pi*f/3e8;         # wave number 
    th_i = np.arctan(R/2/H);     # angle of incidence
    deltaR = 2*np.sqrt(R**2/4+H**2)-R;
    [G_h, G_v]=FresnelCoeff(1, Er, th_i);

    pr = pt*(1+G_v*np.exp(1j*k*deltaR))/(1+G_h*np.exp(1j*k*deltaR));
    return pr


def PolarPattern(p):
    '''
    Parameters
    ----------
    p : polarization phasor

    Returns
    -------
    x : angles [0 2pi]
    V : polarization pattern

    '''
    x = np.linspace(0, 2*np.pi, 361)
    Et_h = np.exp(1j*x);
    Et_v = Et_h*p;
    V = Et_h*np.cos(x)+Et_v*np.sin(x)
    return x, V


def loadmat(filename):
    '''
    this function should be called instead of direct spio.loadmat
    as it cures the problem of not properly recovering python dictionaries
    from mat files. It calls the function check keys to cure all entries
    which are still mat-objects
    '''
    data = spio.loadmat(filename, struct_as_record=False, squeeze_me=True)
    return _check_keys(data)

def _check_keys(dict):
    '''
    checks if entries in dictionary are mat-objects. If yes
    todict is called to change them to nested dictionaries
    '''
    for key in dict:
        if isinstance(dict[key], spio.matlab.mio5_params.mat_struct):
            dict[key] = _todict(dict[key])
    return dict        

def _todict(matobj):
    '''
    A recursive function which constructs from matobjects nested dictionaries
    '''
    dict = {}
    for strg in matobj._fieldnames:
        elem = matobj.__dict__[strg]
        if isinstance(elem, spio.matlab.mio5_params.mat_struct):
            dict[strg] = _todict(elem)
        else:
            dict[strg] = elem
    return dict


#-------------------------  END OF FUNCTIONS ----------------------------------


# Example plot, transmit circular polarization + dielectric plate multipath
grData = loadmat('O:\Running EDUCATION\EE2P1\Session 2.2 Matlab\Group data\group-01.03.mat')
f = grData['session2']['task2']['frequency']
R = grData['session2']['task2']['antennas_distance']
Er = grData['session2']['task2']['dielectric_prermittivity']
H = grData['session2']['task2']['reflection_height']

pt = PolarPhasor(0, np.pi/4)    # transmit circular pol

fig, ax = plt.subplots(subplot_kw={'projection': 'polar'})
for h in H:
    pr = MultipathRxPolarState(pt, f, R, h, Er)
    th, rho = PolarPattern(pr)
    ax.plot(th, np.abs(rho)/np.max(np.abs(rho)),label=round(h,2))



ax.set_rlabel_position(-22.5)  # Move radial labels away from plotted line
ax.grid(True)

ax.set_title("Polarization Pattern", va='bottom')
plt.legend(title="H",loc='best', fontsize=6, fancybox=True)
plt.show()
#PlotPolarPattern(pt)

## Fresnel / Brewster
## glass example (n2 = 1.5)
# th_i = np.linspace(0, np.pi/2, 91)
# Gh, Gv = FresnelCoeff(1, 1.5**2, th_i)
# plt.figure()
# plt.plot(180*th_i/np.pi,np.abs(Gh))
# plt.plot(180*th_i/np.pi,np.abs(Gv))
# plt.show()