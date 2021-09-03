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
import pyraser
from setting_test import Setting as Stest
# sys.path.append("..")
# from geometry import R3dDetector
# from pyfenics import FenicsCal
# from g4particles import Particles
# from calcurrent import CalCurrent
# from elereadout import Amplifier
# from drawsave import drawplot   

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
    dset = Stest(args)
    my_d = pyraser.R3dDetector(dset.detector) 
    my_f = pyraser.FenicsCal(my_d, dset.fenics)
    my_g4p = pyraser.Particles(my_d, my_f, dset)
    my_current = pyraser.CalCurrent(my_d, my_f, my_g4p, dset)
    ele_current = pyraser.Amplifier(my_d, dset.amplifer)
    # pyraser.drawsave.draw_unittest(my_d,ele_current,my_f,my_g4p,my_current)
    judge_test(my_current.laudau_t_pairs,dset.det_model)
    os._exit(0)

def judge_test(e_h_pair,det_model):
    """
    Description:
        Judge whether the raser process is normal
    Parameters:
    ---------
    laudau_t_pairs : float
        Genergy e-h pairs in sensor    
    @Returns:
    ---------
        Flase or True
    @Modify:
    ---------
        2021/09/03
    """
    if ( (int(e_h_pair) == 9754 and det_model == "planar3Dscan") 
        or  (int(e_h_pair) == 15954 and det_model == "plugin3Dscan")):
        print("The raser can work normal. Successful !!!")
        # print("Further compare the 'right_%s.pdf' and 'test.pdf' files "%(det_model),end="") 
        # print("in the 'pyraser/unittest/' path. When the two are roughly ",end="")
        # print("the same, unit test will be more reliable.") 
    else:
        print("The raser can't work. Failed !!!")
        print("Your code make the original version of the code unable to run!!!")

if __name__ == '__main__':
    main() 

