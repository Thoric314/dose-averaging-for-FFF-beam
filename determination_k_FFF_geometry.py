#!/usr/bin/env python
# -*- coding: utf-8 -*-

"""
The k_{FFF} is supposed to correct for the non-homogeneity of the
dose around the detector (here FARMER 2571).
"""

import sys
import getopt
import numpy as np

from uncertainties import ufloat
from uncertainties.umath import *
from math import sqrt
import matplotlib.pyplot as plt

import scipy as sp
import numpy as np

# Around 10 cm depth the PDD is approximated by a line
# slope of this line is :
# a06 for the 06 FFF
# a10 for the 10 FFF

#a06             = -0.0627403       +/- 0.0008565    (1.365%)
#a10             = -0.0521787       +/- 0.0006379    (1.223%)
#k06             = 0.0626403        +/- 0.0003033    (0.4843%)
#k10             = 0.0520527        +/- 0.0003844    (0.7385%)

_debug = False

def usage():
    print()
    print("determination_k_FFF.py")
    print("                        -h, --help : this reading")
    print("                        -d, --debug : to get more details")
    print()
    

def PDD(k, z):
    return exp(k*(10-z))


def Profile(q, x):
    return 1.0 - q*x**2


def dv(n,ra,le,cr):
    return (2/n)**3*ra*le*ra


def plotg(s, mu, sigma, name):
    count, bins, ignored = plt.hist(s, 30, density=True)
    plt.plot(bins, 1/(sigma * np.sqrt(2 * np.pi)) *
             np.exp( - (bins - mu)**2 / (2 * sigma**2) ),
             linewidth=2, color='r')
    plt.title(name)
    #plt.show()
    nameplt=name.replace('$','').replace('\\','').replace('{','').replace('}','')
    plt.savefig(nameplt+'.png')
    plt.close()

def k_FFF(k, q, label):
    n = 10000
    sx = 0.10 # sigma [cm]
    sy = 0.10 # sigma [cm]
    sz = 0.01 # sigma [cm]
    dx = np.random.normal(0.0, sx, n)
    dy = np.random.normal(0.0, sy, n)
    dz = np.random.normal(0.0, sz, n)
    kFFF = np.zeros(n)
   
    for i in range(n):
        vol, kFFF[i] = k_FFF_jitter(k, q, dx[i], dy[i], dz[i])
       
    kFT = kFFF.mean()

    plotg(dx, 0.0, sx, label+"$\Delta_x$")
    plotg(dy, 0.0, sy, label+"$\Delta_y$")
    plotg(dz, 0.0, sz, label+"$\Delta_z$")
    sdev = kFFF.std()
    plotg(kFFF, kFT, sdev, label+"$k_{FFF-S_{w-air}}$")
       
    return vol, kFT
   

def k_FFF_jitter(k, q, dx, dy, dz):
    """
    IC cylinder
    radius = 0.305 cm             0.305 cm
    length = 2.300 cm             1.150 cm
    ceanod = 0.110 cm (diameter)  0.055 cm
    """
    ra = 0.305
    le = 1.150
    cr = 0.055

    decalx = 0.00 
    decaly = 0.09
    decalz = 10.0

    #print()
    #print()
    
    dvol = lambda r, y, th: r
    volume_cylindre = sp.integrate.tplquad(dvol, 0, 2*np.pi, -le+decaly, le-decaly, cr, ra)

    volume_cone = sp.integrate.tplquad(dvol, 0, 2*np.pi,
                                     le-decaly, le+decaly,
                                             0, lambda th, y : (ra*(1-(y-le+decaly)/(2*decaly))) )
    #print("       Volume cylindre : ", volume_cylindre)                                    
    #print("           Volume cone : ", volume_cone)
    vol = volume_cylindre[0] + volume_cone[0]
    #print("Volume par integration : ", vol)
    
    dose = lambda r, y, th: r*PDD(k,r*sin(th)+decalz+dz)*Profile(q,sqrt((r*cos(th)+dx)**2+(y+decaly+dy)**2))

    dose_cylindre = sp.integrate.tplquad(dose, 0, 2*np.pi, -le+decaly, le-decaly, cr, ra)
    dose_cone = sp.integrate.tplquad(dose, 0, 2*np.pi,
                                le-decaly, le+decaly,
                                        0, lambda th, y : (ra*(1-(y-le+decaly)/(2*decaly))) )
                                     
    dosem = (dose_cylindre[0] + dose_cone[0])  / vol
    #print("  Dose par integration : ",dosem)
    #print(" k_FFF par integration : ",1/dosem)
    #print()
    #print()
    return vol, 1/dosem



def calcul():
    print()
    # PDDs
    k06 = ufloat(0.0626403, 0.0003033)
    k10 = ufloat(0.0520527, 0.0003844)
    
    # Profiles X and Y
    # qx06            = 0.00688693       +/- 0.0001704    (2.474%)
    # qy06            = 0.00699005       +/- 0.0003118    (4.46%)
    # qx10            = 0.0125959        +/- 0.0003635    (2.886%)
    # qy10            = 0.0131           +/- 0.000438     (3.343%)

    qx06 = ufloat(0.00688693, 0.0001704)
    qy06 = ufloat(0.00699005, 0.0003118)
    q06 = (qx06+qy06)/2.0

    qx10 = ufloat(0.0125959, 0.0003635)
    qy10 = ufloat(0.0131, 0.000438)

    q10 = (qx10+qy10)/2.0

    print("\n06 FFF")
    print("======")
    print(f"{k06=}")
    print(f"{q06=}")

    volume, kFFF = k_FFF(k06.nominal_value, q06.nominal_value, "06FFF")
    print(f'{volume=} [cm^3]')
    print(f'{kFFF=}')

    print("\n10 FFF")
    print("======")
    print(f"{k10=}")
    print(f"{q10=}")

    volume, kFFF = k_FFF(k10.nominal_value, q10.nominal_value, "10FFF")
    print(f'{volume=} [cm^3]')
    print(f'{kFFF=}')


    
def main(argv):
    try:                                
        opts, args = getopt.getopt(argv, "dhn:", ["debug", "help", "number="])
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()                  
            elif opt == '-d':
                global _debug
                _debug = True
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    calcul()
    print()



if __name__ == "__main__":
    main(sys.argv[1:])
