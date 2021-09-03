#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@Description: The main program of Raser induced current simulation      
@Date       : 2021/08/31 12:14:31
@Author     : tanyuhang
@version    : 1.0
'''
import sys
import os
from setting import Setting
from geometry import R3dDetector
from pyfenics import FenicsCal
from g4particles import Particles
from calcurrent import CalCurrent
from elereadout import Amplifier
from drawsave import drawplot   

def main():
    """
    Description:
        The main program of Raser induced current simulation      
    Parameters:
    ---------
    dset : class
        Parameters of simulation
    Function or class:
        R3dDetector -- Define the structure of the detector
        FenicsCal -- Get the electric field and weighting potential 
        Particles -- Electron and hole paris distibution
        CalCurrent -- Drift of e-h pais and induced current
        Amplifier -- Readout electronics simulation
        drawplot -- Draw electric field ,drift path and energy deposition        
    Modify:
    ---------
        2021/09/02
    """
    args = sys.argv[1:]
    dset = Setting(args)
    my_d = R3dDetector(dset.detector) 
    my_f = FenicsCal(my_d, dset.fenics)
    my_g4p = Particles(my_d, my_f, dset)
    my_current = CalCurrent(my_d, my_f, my_g4p, dset)
    ele_current = Amplifier(my_d, dset.amplifer)
    drawplot(my_d,ele_current,my_f,my_g4p,my_current)
    

if __name__ == '__main__':
    main()
    os._exit(0) 

