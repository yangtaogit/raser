#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Description: 
    Simulate induced current through BB or CSA amplifier 
@Date       : 2021/09/02 14:11:57
@Author     : tanyuhang
@version    : 1.0
'''

import math
import ROOT

class Amplifier:
    def __init__(self,my_d,ampl_par,mintstep="50e-12"):
        CSA_par = ampl_par[0]
        BB_par = ampl_par[1]
        self.max_num = my_d.sum_cu.GetNbinsX()
        self.max_hist_num = my_d.n_bin
        self.undersampling = int(float(mintstep)/my_d.t_bin)
        self.time_unit = my_d.t_bin*self.undersampling
        self.CDet_j = 0

        self.t_rise    = CSA_par['t_rise']
        self.t_fall    = CSA_par['t_fall']
        self.trans_imp = CSA_par['trans_imp']
        self.CDet      = CSA_par['CDet']
        self.BBW       = BB_par['BBW']
        self.BBGain    = BB_par['BBGain']
        self.BB_imp    = BB_par['BB_imp']
        self.OscBW     = BB_par['OscBW'] 

        ##BB simualtion parameter
        tau_C50 = 1.0e-12*50.*self.CDet          #Oscil. RC
        tau_BW = 0.35/(1.0e9*self.OscBW)/2.2      #Oscil. RC
        tau_BB_RC = 1.0e-12*self.BB_imp*self.CDet     #BB RC
        tau_BB_BW = 0.35/(1.0e9*self.BBW)/2.2    #BB Tau
        self.tau_scope = math.sqrt(pow(tau_C50,2)+pow(tau_BW,2))
        self.tau_BBA =  math.sqrt(pow(tau_BB_RC,2)+pow(tau_BB_BW,2))    #BB_out

        self.itot = [0.0]*self.max_num
        self.qtot = 0.0
        # get total charge
        i=0
        for j in range(0,self.max_hist_num,self.undersampling):
            self.itot[i] = my_d.sum_cu.GetBinContent(j)
            self.qtot = self.qtot + self.itot[i]*self.time_unit
            i+=1
        max_hist_num = int(self.max_hist_num/self.undersampling)
        IintTime = max(2.0*(self.t_rise+self.t_fall)*1e-9/self.time_unit,
                       3.0*self.tau_BBA/self.time_unit)
        self.IMaxSh = int(max_hist_num + IintTime)

        self.charge_t=my_d.sum_cu.Integral() \
            * ((my_d.sum_cu.GetXaxis().GetXmax() \
            - my_d.sum_cu.GetXaxis().GetXmin()) \
            / my_d.sum_cu.GetNbinsX()) * 1e15
        self.CSA_amp()
        self.BB_amp()

    def CSA_amp(self):   
        IMaxSh = self.IMaxSh
        t_rise = self.t_rise
        t_fall = self.t_fall
        sh_max = 0.0

        tau_rise = t_rise / 2.2*1e-9
        tau_fall = t_fall / 2.2*1e-9   
        if (tau_rise == tau_fall):
            tau_rise *= 0.9

        preamp_Q     = [0.0]*IMaxSh
        shaper_out_Q = [0.0]*IMaxSh
        shaper_out_V = [0.0]*IMaxSh
        step=1
        for i in range(IMaxSh-step):
            if(i>0 and i <self.max_hist_num-step):
                preamp_Q[i] = 0.0
                for il in range(i,i+step):
                    preamp_Q[i] += self.itot[il]*self.time_unit
            elif (i != 0):
                preamp_Q[i]=0.0

        for i in range(IMaxSh-step):
            if i >= step:
                dif_shaper_Q = preamp_Q[i]
            else:
                dif_shaper_Q = 0
            if (dif_shaper_Q != 0):
                for j in range(IMaxSh-i):
                    shaper_out_Q[i+j] += tau_fall/(tau_fall+tau_rise) \
                                         *dif_shaper_Q*(math.exp(
                                         -j*self.time_unit/tau_fall)
                                         -math.exp(-j*self.time_unit/tau_rise))                
            if (abs(shaper_out_Q[i]) > abs(sh_max)):
                sh_max = shaper_out_Q[i]                
        Ci = 3.5e-11  #fF
        Qfrac = 1.0/(1.0+self.CDet*1e-12/Ci)
        self.CSA_ele = ROOT.TH1F("electronics","electronics",
                                 IMaxSh,0,IMaxSh*self.time_unit)
        for i in range(IMaxSh):
            if sh_max == 0.0:
                shaper_out_V[i] = 0.0
            elif self.CDet_j == 0:
                shaper_out_V[i] = shaper_out_Q[i]*self.trans_imp*\
                                  1e15*self.qtot*Qfrac/sh_max            
            elif self.CDet_j ==1:
                shaper_out_V[i] = shaper_out_Q[i]*self.trans_imp*\
                                  1e15*self.qtot/sh_max
            self.CSA_ele.SetBinContent(i,shaper_out_V[i])

        max_CSA_height = min(shaper_out_V)
        time_t = shaper_out_V.index(max_CSA_height)
        # print("CSA_time=%s"%(time_t*self.time_unit))
        return self.CSA_ele

    def BB_amp(self):
  
        IMaxSh = self.IMaxSh
        preamp_Q = [0.0]*IMaxSh
        Vout_scope = [0.0]*IMaxSh
        Iout_BB_RC = [0.0]*IMaxSh
        Iout_C50 = [0.0]*IMaxSh   
        BBGraph = [0.0]*IMaxSh

        step=1
        self.BB_ele = ROOT.TH1F("electronics BB","electronics BB",
                                IMaxSh,0,IMaxSh*self.time_unit)
        for i in range(IMaxSh-step):
            if(i>0 and i <self.max_hist_num-step):
                preamp_Q[i] = 0.0
                for il in range(i,i+step):
                    preamp_Q[i] += self.itot[il]*self.time_unit
            elif (i != 0):
                preamp_Q[i]=0.0

        for i in range(IMaxSh-step):
            if i >= step:
                dif_shaper_Q = preamp_Q[i]
            else:
                dif_shaper_Q = 0
            # if (dif_shaper_Q != 0):
            for j in range(IMaxSh-i):
                Iout_C50[i+j]   += (dif_shaper_Q)/self.tau_scope \
                                    *math.exp(-j*self.time_unit/self.tau_scope)
                Iout_BB_RC[i+j] += (dif_shaper_Q)/self.tau_BBA \
                                   *math.exp(-j*self.time_unit/self.tau_BBA)
                BBGraph[i+j] =  1e+3*self.BBGain*Iout_BB_RC[i+j]
                Vout_scope[i+j] = 50*Iout_C50[i+j]
                if (abs(BBGraph[i+j])>800):
                    BBGraph[i+j] = 800*BBGraph[i+j]/abs(BBGraph[i+j])
        for i in range(len(BBGraph)+1):
            # random_gauss = ROOT.gRandom.Gaus
            # noise_2=random_gauss(0.293,2.606)
            if i == 0:
                self.BB_ele.SetBinContent(i,0)
            else:
                self.BB_ele.SetBinContent(i,BBGraph[i-1])
        max_BB_height = min(BBGraph)
        time_t =  BBGraph.index(max_BB_height)
        
        # print("BB_time=%s"%(time_t*self.time_unit))

        return self.BB_ele

    def save_ele(self,number,output_path="none"):
        print('number=%s'%number)
        output_file = output_path + "/t_" +str(number)+"_events.csv"
        f1 = open(output_file,"w")
        f1.write("charge:%s fc\n"%(self.qtot*1e15))
        f1.write("time[ns],CSA Amplitude [mV], BB Amplitude [mV] \n")
        for i in range(self.BB_ele.GetNbinsX()):
            f1.write("%s,%s,%s \n"%(i*self.time_unit,
                                    self.CSA_ele[i],
                                    self.BB_ele[i]))
        f1.close()
        print("output_file:%s"%output_file)

        del self.BB_ele
        del self.CSA_ele
    
    def save_charge(self,output_path):
    
        with open(output_path+'.txt','a') as f:
            f.write(str(self.charge_t)+','+str(self.qtot)+'\n')
    def __del__(self):
        pass