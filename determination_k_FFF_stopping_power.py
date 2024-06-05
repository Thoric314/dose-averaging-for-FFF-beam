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

_debug = False

def usage():
    print()
    print("determination_k_FFF_stopping_power.py")
    print("                        -h, --help : this reading")
    print("                        -d, --debug : to get more details")
    print("                        -t, --TPR2010 [float] : Quality Index [DEFAULT=0.627]")
    print()
    


def Lrho_water_air_WFF(t):
    """
    t : TPR^20_10
    Fit de Czarnecki et al. 2017, Medical Physics
    Pour les rapport de pouvoir d'arrêt eau-air avec les cones égalisateurs (WFF)
    """
    value = -0.942215 * t**2 + 1.08739 * t + 0.81529
    return ufloat(value, value*0.0017)

def Lrho_water_air_FFF(t):
    """
    t : TPR^20_10
    Fit de Czarnecki et al. 2017, Medical Physics
    Pour les rapport de pouvoir d'arrêt eau-air sans les cones égalisateurs (FFF)
    """
    value = -0.945130 * t**2 + 1.04130 * t + 0.84365
    return ufloat(value, value*0.0015)

    
def calcul(t):
    print()

    print("\nTPR 20,10")
    print("===========")

    print(f"TPR^20_10 = {t:0.4f}")

    L_WFF =  Lrho_water_air_WFF(t)
    L_FFF =  Lrho_water_air_FFF(t)
    
    print(f'\nL_WFF = {L_WFF:0.4f}')
    print(f'L_FFF = {L_FFF:0.4f}')

    ratio = L_FFF/L_WFF

    print(f'\nCorrection for kQ from METAS fit based on WFF to FFF : {ratio:0.4f}')


    
def main(argv):
    t=0.627
    try:                                
        opts, args = getopt.getopt(argv, "dht:", ["debug", "help", "TPR2010="])
        for opt, arg in opts:
            if opt in ("-h", "--help"):
                usage()
                sys.exit()                  
            elif opt == '-d':
                global _debug
                _debug = True
            elif opt in ("-t", "--TPR2010"):
                t = float(arg)
    except getopt.GetoptError:
        usage()
        sys.exit(2)

    calcul(t)
    print()



if __name__ == "__main__":
    main(sys.argv[1:])

