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
    print("                        -n, --number [num] : to set the number of iterations [DEFAULT=10]")
    print()
    

def PDD(k, z):
    return exp(k*(10-z))


def Profile(q, x):
    return 1.0 - q*x**2


def dv(n,ra,le,cr):
    return (2/n)**3*ra*le*ra


def volume_box_kFFF(ra, le, cr, n, k, q):
    nx=np.linspace(-ra, ra, num=n)
    ny=np.linspace(-le, le, num=n)
    nz=np.linspace(-ra+10, ra+10, num=n)
    v = dv(n,ra,le,cr)
    volume = sum( v
                  for x in nx
                  for y in ny
                  for z in nz )
    ik = sum( v*Profile(q, sqrt(x**2+y**2))*PDD(k, z)
              for x in nx
              for y in ny
              for z in nz )
    ik = ik/volume
    return volume, 1/ik


def volume_cyl_kFFF(ra, le, cr, n, k, q):
    decaly = 0.09
    decalz = 10.0
    nx=np.linspace(-ra, ra, num=n)
    ny=np.linspace(-le+decaly, le+decaly, num=n)
    nz=np.linspace(-ra+decalz, ra+decalz, num=n)
    v = dv(n,ra,le,cr)
    volume = sum( v
                  for x in nx
                  for y in ny
                  for z in nz if x*x + (z-decalz)*(z-decalz) <= ra*ra )
    ik = sum( v*Profile(q, sqrt(x**2+y**2))
               *PDD(k, z)
              for x in nx
              for y in ny
              for z in nz if x*x + (z-decalz)*(z-decalz) <= ra*ra)
    ik = ik/volume
    return volume, 1/ik


def volume_cyl_anode_kFFF(ra, le, cr, n, k, q):
    decaly = 0.09
    decalz = 10.0
    nx=np.linspace(-ra, ra, num=n)
    ny=np.linspace(-le+decaly, le+decaly, num=n)
    nz=np.linspace(-ra+decalz, ra+decalz, num=n)
    v = dv(n,ra,le,cr)
    volume = sum( v
                  for x in nx
                  for y in ny
                  for z in nz
                  if ( x*x + (z-decalz)*(z-decalz) <= ra*ra)
                  and (x*x + (z-decalz)*(z-decalz) > cr*cr) )
    ik = sum( v*Profile(q, sqrt(x**2+y**2))
               *PDD(k, z)
              for x in nx
              for y in ny
              for z in nz
              if (x*x + (z-decalz)*(z-decalz) <= ra*ra)
              and (x*x + (z-decalz)*(z-decalz) > cr*cr) )
    ik = ik/volume
    return volume, 1/ik


def volume_cyl_anode_cone_kFFF(ra, le, cr, n, k, q):
    decaly = 0.09
    decalz = 10.0
    hco = 0.18
    nx = np.linspace(-ra, ra, num=n)
    ny = np.linspace(-le+decaly, le+decaly, num=n)
    nz = np.linspace(-ra+decalz, ra+decalz, num=n)
    v = dv(n,ra,le,cr)

    if _debug:
        for x in nx:
            for y in ny:
                for z in nz:
                    if (  (x**2 + (z-decalz)**2 <= ra*ra) and
                          ( (x**2 + (z-decalz)**2  > cr*cr) or (y >= le-decaly) ) and
                          ( (x**2 + (z-decalz)**2 <= ( ra/hco*(hco+le-decaly-y) )**2 or (y < le-decaly ) ) ) ):
                        print(f'* {x}, {y}, {z}')
    
    volume = sum( v
                  for x in nx
                  for y in ny
                  for z in nz
                  if  (  (x**2 + (z-decalz)**2 <= ra*ra) and
                          ( (x**2 + (z-decalz)**2  > cr*cr) or (y >= le-decaly) ) and
                          ( (x**2 + (z-decalz)**2 <= ( ra/hco*(hco+le-decaly-y) )**2 or (y < le-decaly ) ) ) )
                 )
                  
    ik = sum( v*Profile(q, sqrt(x*x+y*y))
              *PDD(k, z)
              for x in nx
              for y in ny
              for z in nz
              if  (  (x**2 + (z-decalz)**2 <= ra*ra) and
                          ( (x**2 + (z-decalz)**2  > cr*cr) or (y >= le-decaly) ) and
                          ( (x**2 + (z-decalz)**2 <= ( ra/hco*(hco+le-decaly-y) )**2 or (y < le-decaly ) ) ) )
              )

    ik = ik/volume


    return volume, 1/ik


def volume_kFFF(n, k, q):
    """
    IC cylinder
    radius = 0.305 cm             0.305 cm
    length = 2.300 cm             1.150 cm
    ceanod = 0.110 cm (diameter)  0.055 cm
    """
    ra = 0.305
    le = 1.150
    cr = 0.055

    if _debug:
        vo, kF = volume_box_kFFF(ra, le, cr, n, k, q)
        print('\nBox')
        print(f'{vo=} [cm^3]')
        print(f'{kF=}')

    if _debug:
        vo, kF = volume_cyl_kFFF(ra, le, cr, n, k, q)
        print('\nCyl')
        print(f'{vo=} [cm^3]')
        print(f'{kF=}')

    if _debug:
        vo, kF = volume_cyl_anode_kFFF(ra, le, cr, n, k, q)
        print('\nCyl w/o anode')
        print(f'{vo=} [cm^3]')
        print(f'{kF=}')

    vo, kF = volume_cyl_anode_cone_kFFF(ra, le, cr, n, k, q)
    print('\nCyl w/o anode, with cone tip')
    print(f'{vo=} [cm^3]')
    print(f'{kF=}')
    print()
    return vo, kF



def k_FFF(k, q):
    """
    IC cylinder
    radius = 0.305 cm             0.305 cm
    length = 2.300 cm             1.150 cm
    ceanod = 0.110 cm (diameter)  0.055 cm
    """
    ra = 0.305
    le = 1.150
    cr = 0.055
    decaly = 0.09
    decalz = 10.0

    print()
    print()
    
    dvol = lambda r, y, th: r
    volume_cylindre = sp.integrate.tplquad(dvol, 0, 2*np.pi, -le+decaly, le-decaly, cr, ra)

    volume_cone = sp.integrate.tplquad(dvol, 0, 2*np.pi,
                                     le-decaly, le+decaly,
                                             0, lambda th, y : (ra*(1-(y-le+decaly)/(2*decaly))) )
    print("       Volume cylindre : ", volume_cylindre)                                    
    print("           Volume cone : ", volume_cone)
    vol = volume_cylindre[0] + volume_cone[0]
    print("Volume par integration : ", vol)
    
    dose = lambda r, y, th: r*PDD(k,r*sin(th)+decalz)*Profile(q,sqrt((r*cos(th))**2+(y+decaly)**2))

    dose_cylindre = sp.integrate.tplquad(dose, 0, 2*np.pi, -le+decaly, le-decaly, cr, ra)
    dose_cone = sp.integrate.tplquad(dose, 0, 2*np.pi,
                                le-decaly, le+decaly,
                                        0, lambda th, y : (ra*(1-(y-le+decaly)/(2*decaly))) )
                                     
    dosem = (dose_cylindre[0] + dose_cone[0])  / vol
    print("  Dose par integration : ",dosem)
    print(" k_FFF par integration : ",1/dosem)
    
    print()
    print()
    
    return vol, 1/dosem



def calcul(n):
    print()
    print(f'Number of iteration per direction (x,y,z) : {n}')
    # PDDs
    #a06 = ufloat(-0.0627403,0.0008565)
    #b06 = paramb(a06)
    k06 = ufloat(0.0626403, 0.0003033)
    k10 = ufloat(0.0520527, 0.0003844)
        
    #a10 = ufloat(-0.0521787,0.0006379)
    #b10 = paramb(a10)

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

    volume, kFFF = volume_kFFF(n, k06.nominal_value, q06.nominal_value)

    print(f'{volume=} [cm^3]')
    print(f'{kFFF=}')

    volume, kFFF = k_FFF(k06.nominal_value, q06.nominal_value)
    print(f'{volume=} [cm^3]')
    print(f'{kFFF=}')

    
    print("\n10 FFF")
    print("======")
    print(f"{k10=}")
    print(f"{q10=}")

    volume, kFFF = volume_kFFF(n, k10.nominal_value, q10.nominal_value)

    print(f'{volume=} [cm^3]')
    print(f'{kFFF=}')

    volume, kFFF = k_FFF(k10.nominal_value, q10.nominal_value)
    print(f'{volume=} [cm^3]')
    print(f'{kFFF=}')


    
def main(argv):
    n=10
    try:                                
        opts, args = getopt.getopt(argv, "dhn:", ["debug", "help", "number="])
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()                  
            elif opt == '-d':
                global _debug
                _debug = True
            elif opt in ("-n", "--number"):
                n = int(arg)
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    calcul(n)
    print()



if __name__ == "__main__":
    main(sys.argv[1:])

