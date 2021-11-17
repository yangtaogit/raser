# -*- encoding: utf-8 -*-
'''
@Description: Detector structure definition      
@Date       : 2021/08/31 11:09:40
@Author     : tanyuhang
@version    : 1.0
''' 

import ROOT
import math
import sys
from raser.model import Material

#Detector structure
class R3dDetector:
    def __init__(self,dset):
        """
        Description:
            Differnet types detectors parameters assignment.
        Parameters:
        ---------
        det_dic : dictionary
            Contain all detector parameters 
        meter : int
            mater = 1, SiC
            mater = 0, Si          
        Modify:
        ---------
            2021/09/02
        """ 
        det_dic = dset.detector
        self.l_x = det_dic['lx'] 
        self.l_y = det_dic['ly']  
        self.l_z = det_dic['lz'] 
        self.d_neff = det_dic['doping'] 
        self.v_voltage = det_dic['voltage'] 
        self.temperature = det_dic['temp']
        self.det_model = dset.det_model
        self.mater = 1    
        self.current_define()
        if 'plugin3D' in self.det_model:
            self.set_3D_electrode(det_dic['e_ir'],det_dic['e_gap'])

    def current_define(self):
        """
        @description: 
            Parameter current setting     
        @param:
            matter -- 0 is silicon, 1 is silicon carbide
            positive_cu -- Current from holes move
            negative_cu -- Current from electrons move
            sum_cu -- Current from e-h move
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        self.t_bin = 50e-12
        self.t_end = 6.0e-9
        self.t_start = 0
        self.n_bin = int((self.t_end-self.t_start)/self.t_bin)
        
        self.positive_cu = ROOT.TH1F("charge+", "Positive Current",
                                     self.n_bin, self.t_start, self.t_end)
        self.negative_cu = ROOT.TH1F("charge-", "Negative Current",
                                     self.n_bin, self.t_start, self.t_end)
        self.sum_cu = ROOT.TH1F("charge","Total Current",
                                self.n_bin, self.t_start, self.t_end)

    def set_3D_electrode(self,e_ir,e_gap=0):
        """
        @description: 
            3D plug-in detector electrodes setting     
        @param:
            e_ir -- The radius of electrode
            e_gap -- The spacing between the electrodes  
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        self.e_gap = e_gap
        e_r = e_ir  
        e_int = e_gap 
        e_t_y = self.infor_ele(e_r,e_int)
        self.e_tr=[]
        self.e_t_1 = [self.l_x*0.5          ,self.l_y*0.5      ,e_r,0,self.l_z,"n"]
        self.e_t_2 = [self.l_x*0.5-e_int    ,self.l_y*0.5      ,e_r,0,self.l_z,"p"]
        self.e_t_3 = [self.l_x*0.5+e_int    ,self.l_y*0.5      ,e_r,0,self.l_z,"p"]
        self.e_t_4 = [self.l_x*0.5-e_int*0.5,self.l_y*0.5+e_t_y,e_r,0,self.l_z,"p"]
        self.e_t_5 = [self.l_x*0.5+e_int*0.5,self.l_y*0.5+e_t_y,e_r,0,self.l_z,"p"]
        self.e_t_6 = [self.l_x*0.5-e_int*0.5,self.l_y*0.5-e_t_y,e_r,0,self.l_z,"p"]
        self.e_t_7 = [self.l_x*0.5+e_int*0.5,self.l_y*0.5-e_t_y,e_r,0,self.l_z,"p"]
        for i in range(7):
           n_e = eval('self.e_t_' + str(i+1))
           self.e_tr.append(n_e)

    def infor_ele(self,e_r,e_int):
        """
        @description: 
            3D plug-in detector electrodes spacing    
        @param:
            e_x_gap -- Judge whether the electrode is outer the detector
            e_t_y -- Distance between electrodes at y bottom or to and center
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        e_x_gap = self.l_x - 2*e_r - 2*e_int
        if e_x_gap < 0:
            print("the electrode at x position is large than sensor length")
            sys.exit(0)
        e_t_y = math.sqrt(e_int*e_int*0.75)
        if 2*e_t_y > self.l_y:
            print("the electrode at y position is large than sensor length")
            sys.exit(0)            
        return e_t_y



class R2dDetector:

    def __init__(self,det_dic):

        self.det_model = det_dic["name"]
        self.det_width = det_dic["det_width"] #um X-axis
        self.det_thin = det_dic["det_thin"] #um Y-axis
        self.x_step = det_dic["x_step"]
        self.y_step = det_dic["y_step"]
        self.mat = Material(det_dic["material"])
        self.doping_epr = det_dic["doping_epr"]
        self.bias_voltage = det_dic["bias_voltage"]
        self.temperature = det_dic["temperature"]

        self.mat_name = self.mat.mat_name
        self.par_dict = self.mat.mat_database()
        self.nx = int(self.det_width/self.x_step)
        self.ny = int(self.det_thin/self.y_step)

        self.t_start = 0
        self.t_end = 3e-9
        self.t_bin = 50e-12
        self.n_bin = int((self.t_end-self.t_start)/self.t_bin)
        self.positive_cu = ROOT.TH1F("charge+","Positive Current",self.n_bin,0,self.t_end)
        self.negative_cu = ROOT.TH1F("charge-","Negative Current",self.n_bin,0,self.t_end)

        self.gain_positive_cu = ROOT.TH1F("gain_charge+","Current Contribution",self.n_bin,0,self.t_end)
        self.gain_negative_cu = ROOT.TH1F("gain_charge-","Gain Negative Current",self.n_bin,0,self.t_end)

        self.gain_n_n_cu = ROOT.TH1F("gain_n_n","Current Contribution",self.n_bin,0,self.t_end)
        self.gain_n_p_cu = ROOT.TH1F("gain_n_p","Current Contribution",self.n_bin,0,self.t_end)
        self.gain_p_n_cu = ROOT.TH1F("gain_p_n","Current Contribution",self.n_bin,0,self.t_end)
        self.gain_p_p_cu = ROOT.TH1F("gain_p_p","Current Contribution",self.n_bin,0,self.t_end)

        self.sum_cu = ROOT.TH1F("charge","Total Current",self.n_bin,0,self.t_end)
    
    def mesh(self,x_step,y_step):

        self.x_step = x_step
        self.y_step = y_step

        self.nx = int(self.det_width/x_step)
        self.ny = int(self.det_thin/y_step)

    def set_mat(self,mat):
        self.mat =mat
        self.mat_name = mat.mat_name
        self.par_dict = self.mat.mat_database()

    def set_doping(self,doping_epr):
        self.doping_epr = doping_epr

    def set_bias_voltage(self,bias_voltage):
        self.bias_voltage = bias_voltage

    def set_temperature(self,temperature):
        self.temperature = temperature
