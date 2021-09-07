#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Description:  Simulate the e-h pairs drift and calculate induced current     
@Date       : 2021/09/02 14:01:46
@Author     : tanyuhang
@version    : 1.0
'''

import random
import numpy as np
import ROOT
import math
import sys


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
                    self.e_field = my_f.get_e_field(self.d_x,
                                                    self.d_y,
                                                    self.d_z)
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
                        self.wpot = my_f.get_w_p(self.d_cx,self.d_cy,self.d_cz)
                        delta_Uw = (self.wpot 
                                   - my_f.get_w_p(self.d_x,self.d_y,self.d_z))
                        self.charge=self.charg*delta_Uw
                        if(self.v_drift!=0):
                            self.d_time=self.d_time+self.sstep*1e-4/self.v_drift
                            self.path_len+=self.sstep
                        self.d_x=self.d_cx
                        self.d_y=self.d_cy
                        self.d_z=self.d_cz
                        
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
            FF=self.list_add(self.e_field,
                             self.cross(self.e_field,self.BB,self.muhh))
        else:
            FF=self.list_sub(self.e_field,
                             self.cross(self.e_field,self.BB,self.muhe))
        
        total_ef = self.root_mean_square(FF)
        if(total_ef!=0):
            self.delta_x=-self.sstep*self.charg*FF[0]/total_ef
            self.delta_y=-self.sstep*self.charg*FF[1]/total_ef
            self.delta_z=-self.sstep*self.charg*FF[2]/total_ef
        else:
            self.delta_x=0.0
            self.delta_y=0.0
            self.delta_z=0.0

    def cross(self,p1,p2,scale):
        """ Get vector cross product of p1, p2 """
        o1 = p1[1]*p2[2]-p1[2]*p2[1]
        o2 = p1[2]*p2[0]-p1[0]*p2[2]
        o3 = p1[0]*p2[1]-p1[1]*p2[0]
        return [scale*o1,scale*o2,scale*o3]

    def root_mean_square(self,p1):
        " Return root_mean_square of p1"
        return math.sqrt(p1[0]*p1[0]+p1[1]*p1[1]+p1[2]*p1[2])

    def list_add(self,p1,p2):
        " Return the added two lists. eg:[1,2,3]+[1,2,3] = [2,4,6]"
        return [ a+b for a,b in zip(p1,p2) ]

    def list_sub(self,p1,p2):
        " Return the added two lists. eg:[1,2,3]-[1,2,3] = [0,0,0]"
        return [ a-b for a,b in zip(p1,p2) ]

    def drift_v(self,my_d,my_f):
        """ The drift of e-h pairs at electric field """
        e_delta_f = my_f.get_e_field(self.d_x+self.delta_x,
                                     self.d_y+self.delta_y,
                                     self.d_z+self.delta_z)
        te_delta_f = self.root_mean_square(e_delta_f)
        aver_e = (self.root_mean_square(self.e_field) 
                  + te_delta_f)/2.0*1e4            # V/cm                                                        
        mobility = sic_mobility(self.charg,aver_e,my_d)  # mobility cm2/(V s) v : cm/s
        self.v_drift = mobility*aver_e 
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
            if(te_delta_f < DiffOffField):
                self.s_time = self.sstep*1e-4/self.v_drift
                s_sigma = math.sqrt(2.0*self.kboltz*mobility
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


