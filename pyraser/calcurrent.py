#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Description:  Simulate the e-h pairs drift and calculate induced current     
@Date       : 2021/09/02 14:01:46
@Author     : tanyuhang
@version    : 1.0
'''

from array import array
import subprocess
import pyraser
import random
import numpy as np
import ROOT
import math
import sys
import os

  
#The drift of generated particles
class CalCurrent:
    def __init__(self, my_d, my_f, my_g4p, dset):
        #mobility related with the magnetic field (now silicon useless)
        self.muhh=1650.0   
        self.muhe=310.0
        self.BB=np.array([0,0,0])
        self.sstep=dset.steplength #drift step
        self.kboltz=8.617385e-5 #eV/K
        self.max_drift_len=1e9 #maximum diftlenght [um]
        if "3Dscan" in dset.det_model:
            self.scan_batch_drift(my_d, my_f, my_g4p, dset)
        else:
            self.drift(my_d, my_f, my_g4p)
            
    def drift(self, my_d, my_f, my_g4p, batch=0):
        self.d_dic_n = {}
        self.d_dic_p = {}
        if batch == 0:
            for j in range(len(my_g4p.p_steps)):
                if my_g4p.edep_devices[j]>0.2 and batch == 0:
                    self.beam_number = j
                    self.tracks_p = my_g4p.p_steps[j]
                    self.tracks_step_edep = my_g4p.energy_steps[j]
                    self.tracks_t_edep = my_g4p.edep_devices[j]
                    batch+=1
                    break
            if batch == 0:
                print("the sensor did have the hit particles")
                sys.exit()
        else:
            self.beam_number = batch
            self.tracks_p = my_g4p.p_steps[batch]
            self.tracks_step_edep = my_g4p.energy_steps[batch]
            self.tracks_t_edep = my_g4p.edep_devices[batch]
        for n in range(len(self.tracks_p)-1):
            self.d_dic_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            self.d_dic_p["tk_"+str(n+1)] = [ [] for n in range(5) ]
        self.ionized_drift(my_g4p,my_f,my_d)

    def scan_batch_drift(self, my_d, my_f, my_g4p, dset):
        """
        Description:
            Batch run some events to get time resolution
        Parameters:
        ---------
        start_n : int
            Start number of the event
        end_n : int
            end number of the event            
        @Returns:
        ---------
            None
        @Modify:
        ---------
            2021/09/02
        """      
        # from elereadout import Amplifier
        # from drawsave import savedata  

        start_n = dset.intance_number * dset.total_events
        end_n = (dset.intance_number + 1) * dset.total_events
        for event in range(start_n,end_n):
            print("run events number:%s"%event)
            self.drift(my_d, my_f, my_g4p,event)
            ele_current = pyraser.Amplifier(my_d, dset.amplifer)
            pyraser.savedata(my_d,dset.output,event,ele_current)
            del ele_current

    def ionized_drift(self,my_g4p,my_f,my_d):       
        for i in range(len(self.tracks_p)-1):
            self.n_track=i+1
            for j in range(2):
                if (j==0):
                    self.charg=1 #hole
                if (j==1):
                    self.charg=-1 #electron
                self.initial_parameter()
                #generated particles positions
                self.d_x=self.tracks_p[i+1][0]#initial position
                self.d_y=self.tracks_p[i+1][1]
                # the z position need to minus the position from geant4
                self.d_z=self.tracks_p[i+1][2] - my_g4p.init_tz_device 
                while (self.end_cond==0):
                    # electic field of the position
                    self.e_field = np.array(my_f.get_e_field(
                                                             self.d_x,
                                                             self.d_y,
                                                             self.d_z)) 
                    # modify this part										
                    if (self.d_y>=(my_d.l_y-1.0) or self.d_x>=(my_d.l_x-1.0) or self.d_z>=(my_d.l_z-1.0)):
                        self.end_cond=3  
                    elif (self.d_y<=(1.0) or self.d_x<=(1.0) or self.d_z<=(1.0)):
                        self.end_cond=8                    
                    elif (self.e_field[0]==0 and self.e_field[1]==0 and self.e_field[1]==0):
                        self.end_cond=9
                    else:                                                    
                        self.delta_p() #delta_poisiton                    
                        self.drift_v(my_d,my_f) #drift_position                    
                        self.drift_s_step(my_d) #drift_next_posiiton
                        #charge_collection
                        delta_Uw=my_f.get_w_p(self.d_cx,self.d_cy,self.d_cz)\
                                 - my_f.get_w_p(self.d_x,self.d_y,self.d_z)
                        self.charge=self.charg*delta_Uw
                        if(self.v_drift!=0):
                            self.d_time=self.d_time+self.sstep*1e-4/self.v_drift
                            self.path_len+=self.sstep
                        self.d_x=self.d_cx
                        self.d_y=self.d_cy
                        self.d_z=self.d_cz
                        self.wpot=my_f.get_w_p(self.d_x,self.d_y,self.d_z)
                        self.save_inf_track(my_d)  
                        self.drift_end_condition()                                                                         
                    self.n_step+=1
        self.get_current(my_d,my_g4p)

    def initial_parameter(self):
        self.end_cond=0
        self.d_time=1.0e-9
        self.path_len=0
        self.n_step=0
        self.charge=0

    def delta_p(self):
        # magnetic field effect
        if(self.charg)>0:
            FF=self.e_field+self.muhh*np.cross(self.e_field,self.BB)
        else:
            FF=self.e_field-self.muhe*np.cross(self.e_field,self.BB)
        #the delta x with electric field
        if(np.linalg.norm(FF)!=0):
            self.delta_x=-self.sstep*self.charg*FF[0]/np.linalg.norm(FF)
            self.delta_y=-self.sstep*self.charg*FF[1]/np.linalg.norm(FF)
            self.delta_z=-self.sstep*self.charg*FF[2]/np.linalg.norm(FF)
        else:
            self.delta_x=0
            self.delta_y=0
            self.delta_z=0

    def drift_v(self,my_d,my_f):
        e_delta_f = np.array(my_f.get_e_field(self.d_x+self.delta_x,
                                              self.d_y+self.delta_y,
                                            self.d_z+self.delta_z))
        aver_e=(np.linalg.norm(self.e_field)+np.linalg.norm(e_delta_f))/2.0*1e4
                                                                        #V/cm 
         # mobility cm2/(V s) v : cm/s
        self.v_drift=sic_mobility(self.charg,aver_e,my_d)*aver_e 
        #drift part
        if(self.v_drift==0):
            self.delta_x=0.0
            self.delta_y=0.0
            self.delta_z=0.0
            self.dif_x=0.0
            self.dif_y=0.0
            self.dif_z=0.0
            self.end_cond=9
        else:
            #off when the field gets large enough
            DiffOffField=100.0  # if the electric field  
                                # > 100V/um, the holes will multiplicat             
            if(np.linalg.norm(e_delta_f)<DiffOffField):
                self.s_time=self.sstep*1e-4/self.v_drift
                s_sigma=math.sqrt(2.0*self.kboltz*
                                  sic_mobility(self.charg,aver_e,my_d)
                                  *my_d.temperature*self.s_time)
                self.dif_x=random.gauss(0.0,s_sigma)*1e4
                self.dif_y=random.gauss(0.0,s_sigma)*1e4
                self.dif_z=random.gauss(0.0,s_sigma)*1e4       
            else:
                print("the eletric field is too big, \
                       the multiplication appear. The system shold end. ")
                sys.exit(0)

    def drift_s_step(self,my_d):
        # x axis   
        if((self.d_x+self.delta_x+self.dif_x)>=my_d.l_x): 
            self.d_cx = my_d.l_x
        elif((self.d_x+self.delta_x+self.dif_x)<0):
            self.d_cx = 0
        else:
            self.d_cx = self.d_x+self.delta_x+self.dif_x
        # y axis
        if((self.d_y+self.delta_y+self.dif_y)>=my_d.l_y): 
            self.d_cy = my_d.l_y
        elif((self.d_y+self.delta_y+self.dif_y)<0):
            self.d_cy = 0
        else:
            self.d_cy = self.d_y+self.delta_y+self.dif_y
        # z axis
        if((self.d_z+self.delta_z+self.dif_z)>=my_d.l_z): 
            self.d_cz = my_d.l_z
        elif((self.d_z+self.delta_z+self.dif_z)<0):
            self.d_cz = 0
        else:
            self.d_cz = self.d_z+self.delta_z+self.dif_z

    def drift_end_condition(self):    
        if(self.wpot>(1-1e-5)):
            self.end_cond=1
        if(self.d_x<=0):
            self.end_cond=2
        if(self.d_y<=0):
            self.end_cond=4
        if(self.d_z<=0):
            self.end_cond=5
        if(self.path_len>self.max_drift_len):
            self.end_cond=6
        if(self.n_step>10000):
            self.end_cond=7

    def save_inf_track(self,my_d):
        e0 = 1.60217733e-19
        if(((self.charge<0 and my_d.v_voltage<0) 
             or (self.charge>0 and my_d.v_voltage>0))): 
            if(self.charg>0):
                self.d_dic_p["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_p["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_p["tk_"+str(self.n_track)][2].append(self.d_z)
                self.d_dic_p["tk_"+str(self.n_track)][3].append(self.charge)
                self.d_dic_p["tk_"+str(self.n_track)][4].append(self.d_time)
            else:
                self.d_dic_n["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_n["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_n["tk_"+str(self.n_track)][2].append(self.d_z)
                self.d_dic_n["tk_"+str(self.n_track)][3].append(self.charge)
                self.d_dic_n["tk_"+str(self.n_track)][4].append(self.d_time)

    def get_current(self,my_d,my_g4p):
        my_d.positive_cu.Reset()
        my_d.negtive_cu.Reset()
        my_d.sum_cu.Reset()
        self.sum_p_current = []
        test_p = ROOT.TH1F("test+","test+",my_d.n_bin,my_d.t_start,my_d.t_end)
        test_n = ROOT.TH1F("test-","test-",my_d.n_bin,my_d.t_start,my_d.t_end)
        test_sum = ROOT.TH1F("test sum","test sum",
                             my_d.n_bin,my_d.t_start,my_d.t_end)
        total_pairs=0
        if (my_d.mater == 1): # silicon carbide
            sic_loss_e = 8.4 #ev
        elif (my_d.mater == 0):   # silicon
            sic_loss_e = 3.6 #ev
        for j in range(len(self.tracks_p)-1):
            for i in range(len(self.d_dic_p["tk_"+str(j+1)][2])):
                test_p.Fill(self.d_dic_p["tk_"+str(j+1)][4][i],
                            self.d_dic_p["tk_"+str(j+1)][3][i])
            test_p = self.get_current_his(test_p)           
            for i in range(len(self.d_dic_n["tk_"+str(j+1)][2])):
                test_n.Fill(self.d_dic_n["tk_"+str(j+1)][4][i],
                            self.d_dic_n["tk_"+str(j+1)][3][i])
            test_n = self.get_current_his(test_n)
            n_pairs=self.tracks_step_edep[j]*1e6/sic_loss_e
            total_pairs+=n_pairs
            test_p.Scale(n_pairs)
            test_n.Scale(n_pairs)            
            my_d.positive_cu.Add(test_p)
            my_d.negtive_cu.Add(test_n)
            test_p.Reset()
            test_n.Reset()
        self.laudau_t_pairs = self.tracks_t_edep*1e6/sic_loss_e
        if total_pairs != 0:
            n_scale = self.laudau_t_pairs/total_pairs
        else:
            n_scale=0
        my_d.sum_cu.Add(my_d.positive_cu)
        my_d.sum_cu.Add(my_d.negtive_cu)
        my_d.sum_cu.Scale(n_scale)

    def get_current_his(self,histo):
        e0 = 1.60217733e-19
        hist = ROOT.TH1F()
        histo.Copy(hist)
        for i in range(hist.GetNbinsX()):
            histo.SetBinContent(i, hist.GetBinContent(i) \
                /((hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) \
                / hist.GetNbinsX())*e0)
        return histo

    def draw_drift_path(self,my_d,my_f,sensor_model):
        ROOT.gStyle.SetOptStat(0)
        # # ROOT.gROOT.SetBatch(1)
        c1 = ROOT.TCanvas("c", "canvas1", 200,10,1000, 1000)
        c1.Divide(1,2)

        if "plugin3D" in sensor_model:
            n_bin=[int((my_f.sx_r-my_f.sx_l)/5),
                   int((my_f.sy_r-my_f.sy_l)/5),int((my_d.l_z)/10)]
            structrue = ROOT.TH3D("","",n_bin[0],my_f.sx_l,my_f.sx_r,
                                        n_bin[1],my_f.sy_l,my_f.sy_r,
                                        n_bin[2],0,my_d.l_z)
        elif "planar3D" in sensor_model:
            n_bin=[int(my_d.l_x/50),int(my_d.l_y/50),int(my_d.l_z)]
            structrue = ROOT.TH3D("","",n_bin[0],0,my_d.l_x,
                                        n_bin[1],0,my_d.l_y,
                                        n_bin[2],0,my_d.l_z)
        c1.cd(1)
        for k in range(n_bin[2]):
            for j in range (n_bin[1]):
                for i in range(n_bin[0]):
                    if "plugin3D" in sensor_model:
                        x_v = (i+1)*((my_f.sx_r-my_f.sx_l)/n_bin[0])+my_f.sx_l
                        y_v = (j+1)*((my_f.sx_r-my_f.sx_l)/n_bin[1])+my_f.sx_l
                        z_v = (k+1)*(my_d.l_z/n_bin[2])
                    elif "planar3D" in sensor_model:
                        x_v = (i+1)*(my_d.l_x/n_bin[0])
                        y_v = (j+1)*(my_d.l_y/n_bin[1])
                        z_v = (k+1)*(my_d.l_z/n_bin[2])
                    try:
                        x_value,y_value,z_value = my_f.get_e_field(x_v,y_v,z_v)
                        if x_value==0 and y_value==0 and z_value ==0:
                            structrue.SetBinContent(i+1,j+1,k+1,1)
                        else:                       
                            structrue.SetBinContent(i+1,j+1,k+1,0)
                    except RuntimeError:
                        structrue.SetBinContent(i+1,j+1,k+1,1)
        structrue.SetFillColor(1)
        structrue.Draw("ISO")

        mg = ROOT.TMultiGraph("mg","")
        x_array=array('f')
        y_array=array('f')
        z_array=array('f')
        for i in range(len(self.d_dic_p)):
            n=len(self.d_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.d_dic_p["tk_"+str(i+1)][0])
                y_array.extend(self.d_dic_p["tk_"+str(i+1)][1]) 
                z_array.extend(self.d_dic_p["tk_"+str(i+1)][2])              
                gr_p = ROOT.TPolyLine3D(n,x_array,y_array,z_array)
                gr_p.SetLineColor(2)
                gr_p.SetLineStyle(1)
                gr_p.Draw("SAME")
                gr_2D_p=ROOT.TGraph(n,x_array,z_array)
                gr_2D_p.SetMarkerColor(2)
                gr_2D_p.SetLineColor(2)
                gr_2D_p.SetLineStyle(1)
                mg.Add(gr_2D_p)
                del x_array[:]
                del y_array[:]
                del z_array[:]
        for j in range(len(self.d_dic_n)):
            m=len(self.d_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.d_dic_n["tk_"+str(j+1)][0])
                y_array.extend(self.d_dic_n["tk_"+str(j+1)][1])
                z_array.extend(self.d_dic_n["tk_"+str(j+1)][2])                
                gr_n = ROOT.TPolyLine3D(m,x_array,y_array,z_array)
                gr_n.SetLineColor(4)
                gr_n.SetLineStyle(1)
                gr_n.Draw("SAME")
                gr_2D_n=ROOT.TGraph(m,x_array,z_array)
                gr_2D_n.SetMarkerColor(4)
                gr_2D_n.SetLineColor(4)
                gr_2D_n.SetLineStyle(1)
                mg.Add(gr_2D_n)
                del x_array[:]
                del y_array[:]
                del z_array[:]
        c1.cd(2)
        mg.Draw("APL")
        mg.GetXaxis().SetTitle("x aixs")
        mg.GetYaxis().SetTitle("z aixs")
        c1.SaveAs("fig/drift_path.root")
        del c1
        
# # # mobility model
def sic_mobility(charge,aver_e,my_d):
    T=my_d.temperature
    E=aver_e
    Neff=abs(my_d.d_neff)
    #silicon
    if my_d.mater == 0:
        alpha = 0.72*math.pow(T/300.0,0.065)
        if(charge>0):
            ulp = 460.0 * math.pow(T / 300.0, -2.18)
            uminp = 45.0*math.pow(T / 300.0, -0.45)
            Crefp = 2.23e17*math.pow(T / 300.0, 3.2)
            betap = 1.0
            vsatp = 9.05e6 * math.sqrt(math.tanh(312.0/T))
            lfm = uminp + (ulp-uminp)/(1.0 + math.pow(Neff*1e12 / Crefp, alpha))
            hfm = 2*lfm / (1.0+math.pow(1.0 + math.pow(2*lfm * E / vsatp, betap), 1.0 / betap))                        
        else:
            uln = 1430.0 * math.pow(T / 300.0, -2.0)
            uminn = 80.0*math.pow(T / 300.0, -0.45)
            Crefn = 1.12e17*math.pow(T/300.0,3.2)
            betan = 2
            vsatn = 1.45e7 * math.sqrt(math.tanh(155.0/T))
            lfm = uminn + (uln-uminn)/ (1.0 + math.pow(Neff*1e12 / Crefn, alpha))
            hfm = 2*lfm / (1.0+math.pow(1.0 + math.pow(2*lfm * E / vsatn, betan), 1.0/betan))
    #silicon carbide
    elif my_d.mater == 1:
        if(charge>0):
            alpha = 0.34
            ulp = 124.0 * math.pow(T / 300.0, -2.0)
            uminp = 15.9
            Crefp = 1.76e19
            betap = 1.213 * math.pow(T / 300.0, 0.17)
            vsatp = 2e7 * math.pow(T / 300.0, 0.52)
            lfm = uminp + ulp/(1.0 + math.pow(Neff*1e12 / Crefp, alpha))
            hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0 / betap))                        
        else:
            alpha = 0.61
            ulp = 947.0 * math.pow(T / 300.0, -2)
            Crefp = 1.94e19
            betap = 1.0 * math.pow(T / 300.0, 0.66)
            vsatp = 2e7 * math.pow(T / 300.0, 0.87)
            lfm = ulp/ (1.0 + math.pow(Neff*1e12 / Crefp, alpha))
            hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0/betap))                      
    return hfm

def runcmd(command):
    """ subprocess run the shell command """
    ret = subprocess.run([command],shell=True)