#!/usr/bin/env python3

"""
author: tanyuhang
time: 2021.3.8
Use: 1.Read the data of Raser induced current
     2.Add the noise on induced current
     3.Get the time resolution with constance CFD   
"""

from array import array
import contextlib
import os
import sys
import re
import ROOT
import math

# ROOT file parameters difinition
Events=array('i',[0])
h_pulse_time=ROOT.std.vector(float)()
h_BB_pulse_height=ROOT.std.vector(float)()
h_BB_nps_height=ROOT.std.vector(float)()
h_BB_max_nps_height=array('f',[0])
h_BB_max_pulse_time=array('f',[0])
h_BB_time_resolution=array('f',[0])

h_CSA_pulse_height=ROOT.std.vector(float)()
h_CSA_nps_height=ROOT.std.vector(float)()
h_CSA_max_nps_height=array('f',[0])
h_CSA_max_pulse_time=array('f',[0])
h_CSA_time_resolution=array('f',[0])

h_CSA_dVdt=array('f',[0])
h_BB_dVdt=array('f',[0])

h_noise_CSA_height_jitter=array('f',[0])
h_noise_BB_height_jitter=array('f',[0])
h_noise_height_RMS=array('f',[0])
h_noise_height=ROOT.std.vector(float)()

h_pulse_time.clear()
h_CSA_pulse_height.clear()
h_CSA_nps_height.clear()
h_BB_pulse_height.clear()
h_BB_nps_height.clear()
h_CSA_max_nps_height[0]=0.0
h_CSA_max_pulse_time[0]=0.0
h_CSA_time_resolution[0]=0.0

h_BB_max_nps_height[0]=0.0
h_BB_max_pulse_time[0]=0.0
h_BB_time_resolution[0]=0.0
h_CSA_dVdt[0]=0.0
h_BB_dVdt[0]=0.0

h_noise_BB_height_jitter[0]=0.0
h_noise_CSA_height_jitter[0]=0.0
h_noise_height_RMS[0]=0.0
h_noise_height.clear()


# The judge parameter configuration and read data
class Settting:
    def __init__(self):
        pass
    
    def loop_out_p(self):
        """
        @description: The judge parameter configuration
         
        @param:
            thre_vth - Signal threshold
            CFD - Constant ratio timing ratio		 
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        self.thre_vth=18 # mv
        self.CFD=0.5
        self.i=1
        self.CFD_BB_time=[]
        self.CFD_CSA_time=[]

    def create_outpath(self,path):
        """
        @description: If the path is not exit, then create the path
         
        @param:
            None
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        if not os.access(path, os.F_OK):
            os.makedirs(path)

    @contextlib.contextmanager
    def open_func(self,file_name):
        """
        @description: Open file with context manager
         
        @param:
            None
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        print('open file:', file_name, 'in __enter__')
        file_handler = open(file_name, 'r')

        try:
            yield file_handler
        except Exception as exc:
            print('the exception was thrown')
        finally:
            print('close file:', file_name, 'in __exit__')
            file_handler.close()
        return

    def write_list(self,path,list_c):
        """
        @description: Save the file contents in list_c
         
        @param:
            None
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        with self.open_func(path) as file_in:
            for line in file_in:
                if not (is_number(line.split(",")[0])):
                    continue
                list_c.append(line)


# Add noise in the wavefroms and save the data in dictionaries
class AddNoise:
    def __init__(self):
        self.time=0.0
        self.ampl_CSA_nps = 0.0
        self.ampl_CSA_s = 0.0
        self.ampl_BB_nps = 0.0
        self.ampl_BB_s = 0.0
        self.time_list = []
        self.ampl_paras = {}
        self.ampl_CSA_nps_list = []
        self.ampl_BB_nps_list = []
        self.ampl_CSA_s_list = []
        self.ampl_BB_s_list = []
        self.noise_height_list = []

        self.list_c = []
        self.BB_Fv = 0
        self.CSA_Fv = 0        
        self.dVdt={"BB":0, "CSA":0}
        self.CFD_time_r = {"BB":0, "CSA":0}
        self.noist_height_jitter = {"BB":0, "CSA":0}

    def add_n(self,list_c):
        """
        @description: 
            Add gauss noise from data fitting at the waveforms
        @param:
            nps -- Noise plus signal
            s -- Signal
            CSA -- Charge sensitive amplifier
            BB -- Current sensitive amplifier
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        ROOT.gRandom.SetSeed(0)
        random_gauss = ROOT.gRandom.Gaus
        for j in range (0,len(list_c)):
            time= float(list(filter(None,list_c[j].split(",")))[0])
            noise_height=random_gauss(-0.133,2.671)
            ampl_CSA_nps=-float(list(filter(None,list_c[j].split(",")))[1])+noise_height
            ampl_BB_nps=-float(list(filter(None,list_c[j].split(",")))[2])+noise_height
            ampl_CSA_s=-float(list(filter(None,list_c[j].split(",")))[1])
            ampl_BB_s=-float(list(filter(None,list_c[j].split(",")))[2])
            self.time_list.append(time)
            self.noise_height_list.append(noise_height)
            self.ampl_CSA_nps_list.append(ampl_CSA_nps)
            self.ampl_BB_nps_list.append(ampl_BB_nps)
            self.ampl_CSA_s_list.append(ampl_CSA_s)
            self.ampl_BB_s_list.append(ampl_BB_s)
        self.noise_height_RMS= math.sqrt(sum([x**2 
                                              for x in self.noise_height_list])
                                              /len(self.noise_height_list))
        h_noise_height_RMS[0]=self.noise_height_RMS
        self.get_max()
        self.fill_dic()

    def get_max(self):
        """
        @description: 
            Get the max waveform height and height index
        @param:
            None
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        self.max_CSA_nps_height=max(self.ampl_CSA_nps_list)
        self.max_CSA_index=self.ampl_CSA_nps_list.index(max(self.ampl_CSA_nps_list))
        self.max_CSA_pulse_time=self.time_list[self.max_CSA_index]
        self.max_CSA_s_height=max(self.ampl_CSA_s_list)
        self.max_CSA_s_index=self.ampl_CSA_s_list.index(max(self.ampl_CSA_s_list))
        self.max_CSA_s_time=self.time_list[self.max_CSA_s_index]

        self.max_BB_nps_height=max(self.ampl_BB_nps_list)
        self.max_BB_index=self.ampl_BB_nps_list.index(max(self.ampl_BB_nps_list))
        self.max_BB_pulse_time=self.time_list[self.max_BB_index]
        self.max_BB_s_height=max(self.ampl_BB_s_list)
        self.max_BB_s_index=self.ampl_BB_s_list.index(max(self.ampl_BB_s_list))
        self.max_BB_s_time=self.time_list[self.max_BB_s_index]

    def fill_dic(self):
        """
        @description: 
            Fill the parameter of BB and CSA in the ampl_paras dictionaries
        @param:
            None
        @Returns:
            None
        @Modify:
            2021/08/31
        """
        self.ampl_paras["max_CSA_nps_height"] = self.max_CSA_nps_height
        self.ampl_paras["max_CSA_pulse_time"] = self.max_CSA_pulse_time
        self.ampl_paras["ampl_CSA_nps_list"] =  self.ampl_CSA_nps_list
        self.ampl_paras["ampl_CSA_s_list"] =    self.ampl_CSA_s_list
        self.ampl_paras["max_CSA_s_height"] =   self.max_CSA_s_height
        self.ampl_paras["max_CSA_s_time"] =     self.max_CSA_s_time

        self.ampl_paras["max_BB_nps_height"] = self.max_BB_nps_height
        self.ampl_paras["max_BB_pulse_time"] = self.max_BB_pulse_time
        self.ampl_paras["ampl_BB_nps_list"] =  self.ampl_BB_nps_list
        self.ampl_paras["ampl_BB_s_list"] =    self.ampl_BB_s_list
        self.ampl_paras["max_BB_s_height"] =   self.max_BB_s_height
        self.ampl_paras["max_BB_s_time"] =     self.max_BB_s_time
        self.ampl_paras["time_list"] =         self.time_list

        self.ampl_paras["noise_height_list"] = self.noise_height_list


# Root file init definition and fill
class RootFile:
    def init_parameter(self):
        h_BB_max_nps_height[0] = 0.0
        h_BB_max_pulse_time[0] = 0.0  
        h_BB_time_resolution[0] = 0.0
        h_CSA_max_nps_height[0] = 0.0
        h_CSA_max_pulse_time[0] = 0.0  
        h_CSA_time_resolution[0] = 0.0
        h_noise_BB_height_jitter[0] = 0.0
        h_noise_CSA_height_jitter[0] = 0.0
        h_BB_dVdt[0] = 0.0
        h_CSA_dVdt[0] = 0.0
        h_pulse_time.clear()
        h_BB_nps_height.clear()
        h_CSA_nps_height.clear()	
        h_BB_pulse_height.clear()
        h_CSA_pulse_height.clear()	
        h_noise_height.clear()

    def root_define(self):
        """ Root tree branceh definition """
        self.tree_out=ROOT.TTree('tree','tree')
        self.tree_out.Branch('Events',Events,'Events/I')
        self.tree_out.Branch('h_pulse_time',h_pulse_time)
        self.tree_out.Branch('h_BB_pulse_height',h_BB_pulse_height)  # BB
        self.tree_out.Branch('h_BB_nps_height',h_BB_nps_height)  
        self.tree_out.Branch('h_BB_dVdt',h_BB_dVdt,'h_BB_dVdt/F')
        self.tree_out.Branch('h_CSA_pulse_height',h_CSA_pulse_height)  # CSA
        self.tree_out.Branch('h_CSA_nps_height',h_CSA_nps_height) 
        self.tree_out.Branch('h_CSA_dVdt',h_CSA_dVdt,'h_CSA_dVdt/F')
        self.tree_out.Branch('h_noise_height',h_noise_height)  # noise

        self.tree_out.Branch('h_BB_max_nps_height',
                              h_BB_max_nps_height,'h_BB_max_nps_height/F')    
        self.tree_out.Branch('h_BB_max_pulse_time',
                              h_BB_max_pulse_time,'h_BB_max_pulse_time/F')
        self.tree_out.Branch('h_BB_time_resolution',
                              h_BB_time_resolution,'h_BB_time_resolution/F')
        self.tree_out.Branch('h_CSA_max_nps_height',
                              h_CSA_max_nps_height,'h_CSA_max_nps_height/F')    
        self.tree_out.Branch('h_CSA_max_pulse_time',
                              h_CSA_max_pulse_time,'h_CSA_max_pulse_time/F')
        self.tree_out.Branch('h_CSA_time_resolution',
                              h_CSA_time_resolution,'h_CSA_time_resolution/F')

        self.tree_out.Branch('h_noise_CSA_height_jitter',
                              h_noise_CSA_height_jitter,
                             'h_noise_CSA_height_jitter/F')
        self.tree_out.Branch('h_noise_BB_height_jitter',
                              h_noise_BB_height_jitter,
                             'h_noise_BB_height_jitter/F')
        self.tree_out.Branch('h_noise_height_RMS',
                              h_noise_height_RMS,
                             'h_noise_height_RMS/F')
    
    def fill_ampl_BB(self,addNoise,rset,max_height,max_time):
        """ Fill parameters from BB """
        h_BB_max_nps_height[0]=addNoise.ampl_paras[max_height]
        h_BB_max_pulse_time[0]=addNoise.ampl_paras[max_time]                      
        h_BB_time_resolution[0]=addNoise.CFD_time_r["BB"]
        h_noise_BB_height_jitter[0]=addNoise.noist_height_jitter["BB"]
        h_BB_dVdt[0]=addNoise.dVdt["BB"]
        rset.CFD_BB_time.append(addNoise.CFD_time_r["BB"]) 
        addNoise.BB_Fv=1

    def fill_ampl_CSA(self,addNoise,rset,max_height,max_time):
        """ Fill parameters from CSA """
        h_CSA_max_nps_height[0]=addNoise.ampl_paras[max_height]
        h_CSA_max_pulse_time[0]=addNoise.ampl_paras[max_time]                      
        h_CSA_time_resolution[0]=addNoise.CFD_time_r["CSA"]
        h_noise_CSA_height_jitter[0]=addNoise.noist_height_jitter["CSA"]
        h_CSA_dVdt[0]=addNoise.dVdt["CSA"]
        rset.CFD_CSA_time.append(addNoise.CFD_time_r["CSA"]) 
        addNoise.CSA_Fv=1 

    def fill_verctor(self,rset,addNoise):
        """ Fill parameters from CSA and CSA vector """ 
        if (addNoise.BB_Fv==1 or addNoise.CSA_Fv==1):
            for j in range(0,len(addNoise.time_list)):
                h_pulse_time.push_back(addNoise.time_list[j])
                h_BB_pulse_height.push_back(addNoise.ampl_paras["ampl_BB_s_list"][j])
                h_CSA_pulse_height.push_back(addNoise.ampl_paras["ampl_CSA_s_list"][j])
                h_BB_nps_height.push_back(addNoise.ampl_paras["ampl_BB_nps_list"][j])
                h_CSA_nps_height.push_back(addNoise.ampl_paras["ampl_CSA_nps_list"][j])
                h_noise_height.push_back(addNoise.ampl_paras["noise_height_list"][j])


def judge_threshold(addNoise,rset,tree_class,model):
    """
    @description: 
        Judge the pluse height of waveform is larger than threshold or not.
        If Yes, the parameters of the waveform will be saved.
    @param:
        None
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    max_height = "max_repampl_nps_height".replace("repampl",model)
    max_time ="max_repampl_pulse_time".replace("repampl",model)
    if (addNoise.ampl_paras[max_height]>rset.thre_vth and addNoise.ampl_paras[max_time]<80):
        get_CFD_time(addNoise,addNoise.ampl_paras,rset,model)
        if model == "BB":
            tree_class.fill_ampl_BB(addNoise,rset,max_height,max_time)
        elif model == "CSA":
            tree_class.fill_ampl_CSA(addNoise,rset,max_height,max_time)


def get_CFD_time(addNoise,Ampl_paras,rset,model):
    """
    @description: 
        Get the time resolution for constance CFD value
        Time resolution is jitter + Landau fluctuation
    @param:
        None
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    noise_from_sensor_capacitance()
    random_gauss = ROOT.gRandom.Gaus

    CFD_time=0.0
    jitter=0.0
    dVdt=0.0
    time_list=[]
    time_list=addNoise.time_list
    model_list="ampl_repampl_nps_list".replace("repampl",model)
    max_time="max_repampl_pulse_time".replace("repampl",model)
    model_height="max_repampl_nps_height".replace("repampl",model)
   
    for i in range (0,len(time_list)):
        if Ampl_paras[model_list][i]>=Ampl_paras[model_height]*rset.CFD \
           and time_list[i]<Ampl_paras[max_time] \
           and time_list[i+2]<Ampl_paras[max_time] and time_list[i-2]>1.0e-9:
            
            dVdt=(Ampl_paras[model_list][i+2]
                  -Ampl_paras[model_list][i-2]) \
                  /(time_list[i+2]-time_list[i-2])/1e9/1.38       
            if (dVdt!=0):
                jitter=addNoise.noise_height_RMS/dVdt
                # 3.85 is the initial time by personal customization 
                CFD_time = 3.85+time_list[i]*1e9+random_gauss(0,jitter)
            else:
                CFD_time=0
            break
    # Fill the parameter to the AddNoise class
    addNoise.CFD_time_r[model]=CFD_time
    addNoise.dVdt[model]=dVdt
    addNoise.noist_height_jitter[model]=jitter        


def noise_from_sensor_capacitance():
    """
    @description: 
        This is the noise from the sensor capacitance.
        The model does not complete.
    @param:
        None
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    perm_sic = 9.76  # Permittivity SiC
    DCap = 5000*5000/100*perm_sic*8.851e-3 # fF backplane
    DCap += 0.014*perm_sic*4*5000  # fF n layer
    DCap +=50 # fF fixed 
    noise_sen = 2.0*DCap/math.sqrt(1) 


def save_waveform_threshold(output_file,event_n,addNoise):
    """ Save waveform in the outputfile """ 
    print(output_file)
    output_path = output_file + "out_txt/"
    if not os.access(output_path, os.F_OK):
        os.makedirs(output_path) 
    f1 = open(output_path+"t_"+str(event_n)+".csv","w")
    f1.write("time[ns],CSA Amplitude [mV], BB Amplitude [mV] \n")
    for i in range(len(addNoise.ampl_paras["time_list"])):
        time = addNoise.ampl_paras["time_list"][i]
        BB =   addNoise.ampl_paras["ampl_BB_nps_list"][i]
        CSA =  addNoise.ampl_paras["ampl_CSA_nps_list"][i]
        f1.write("%s,%s,%s \n"%(time,BB,CSA))
    f1.close()
    return event_n+1        


def FormatLegend(leg):
    """ ROOT Lengend setting """ 
    leg.SetBorderSize(1)
    leg.SetTextFont(43)
    leg.SetTextSize(40)
    leg.SetFillStyle(1)
    leg.SetFillColor(1)
    leg.SetLineColor(2) 


def set_color_marker(color,marker,i,gr):
    """ ROOT color and marker setting """ 
    f=marker[i]
    gr.SetMarkerStyle(f)
    gr.SetMarkerSize(1)
    k=color[i]
    gr.SetLineColor(k)
    gr.SetLineWidth(2)
    gr.SetMarkerColor(k)
    return gr


def fill_legend(leg,gr,name):
    """ Fill graph name in lengend """ 
    leg.AddEntry(gr,name,"LP")
    return leg

def is_number(s):
    """ 
    Define the input s is number or not.
    if Yes, return True, else return False.
    """ 
    try:
        float(s)
        return True
    except ValueError:
        pass
    try:
        import unicodedata
        unicodedata.numeric(s)
        return True
    except (TypeError, ValueError):
        pass
    return False


def draw_2D_CFD_time(CFD_time,out_put,model):
    """
    @description: 
        Draw and fit time distribution.
        Get the time resolution with constant CFD
    @param:
        None
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    c1 =  ROOT.TCanvas("c1"+model,"c1"+model,200,10,800,600)
    c1.SetGrid()
    c1.SetLeftMargin(0.2)
    c1.SetTopMargin(0.12)
    c1.SetBottomMargin(0.2)
    # Define lengend th1f and root gstyle
    leg = ROOT.TLegend(0.25, 0.6, 0.40, 0.8)
    histo=ROOT.TH1F("","",20,4.6,5.6)
    root_set()
    # Get the data
    for i in range(0,len(CFD_time)):
        if CFD_time[i]>0:
            histo.Fill(CFD_time[i])
    # Fit data
    fit_func_1,sigma,error=fit_data(histo)
    histo=th1f_define(histo)
    # Lengend setting
    leg.AddEntry(fit_func_1,"Fit","L")
    leg.AddEntry(histo,"Sim","L")
    # Draw
    histo.Draw()
    fit_func_1.Draw("same")
    leg.Draw("same")
    # Text set
    root_tex(sigma,error,model)
    # Save
    c1.SaveAs(out_put+model+".pdf")
    del c1


def fit_data(histo):
    """ Fit data of time distribution and get the time resolution """
    fit_func_1 = ROOT.TF1('fit_func_1','gaus',4.6,5.6)
    histo.Fit("fit_func_1","ROQ+","",4.6,5.6)
    print("constant:%s"%fit_func_1.GetParameter(0))
    print("constant_error:%s"%fit_func_1.GetParError(0))
    print("mean:%s"%fit_func_1.GetParameter(1))
    print("mean_error:%s"%fit_func_1.GetParError(1))
    print("sigma:%s"%fit_func_1.GetParameter(2))
    print("sigma_error:%s"%fit_func_1.GetParError(2))
    sigma=fit_func_1.GetParameter(2)*1000
    error=fit_func_1.GetParError(2)*1000
    fit_func_1.SetLineWidth(2)
    return fit_func_1,sigma,error


def th1f_define(histo):
    """ TH1f definition """
    histo.GetXaxis().SetTitle("ToA [ns]")
    histo.GetYaxis().SetTitle("Events")
    histo.GetXaxis().CenterTitle()
    histo.GetYaxis().CenterTitle()
    histo.GetXaxis().SetTitleOffset(1.2)
    histo.GetXaxis().SetTitleSize(0.07)
    histo.GetXaxis().SetLabelSize(0.05)
    histo.GetXaxis().SetNdivisions(510)
    histo.GetYaxis().SetTitleOffset(1.1)
    histo.GetYaxis().SetTitleSize(0.07)
    histo.GetYaxis().SetLabelSize(0.05)
    histo.GetYaxis().SetNdivisions(505)
    histo.GetXaxis().CenterTitle()
    histo.GetYaxis().CenterTitle()
    histo.SetLineWidth(2)
    return histo


def root_tex(sigma,error,model):
    """　Latex definition """
    tex = ROOT.TLatex()
    tex.SetNDC(1)
    tex.SetTextFont(43)
    tex.SetTextSize(25)
    tex.DrawLatexNDC(0.6, 0.7, "CFD=0.5"+" "+model+"ampl")
    tex.DrawLatexNDC(0.6, 0.6, "#sigma = %.0f #pm %.0f ps"%(sigma,error))


def root_set():
    """ ROOT gstyle setting"""
    ROOT.gStyle.SetOptFit()
    ROOT.gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(1)


# Loop and add noise in the raser
def loop_addNoise(input_file,rset,tree_class,out_file):
    for root,dirs,files in os.walk(input_file):
        for file in files:    
            if rset.i<100000:    
                print(file) 
                print("................events:%s..............."%(Events[0])) 
                print("................Save events:%s..............."%rset.i)
                path = os.path.join(input_file, file)
                Events[0]+=1

                addNoise = AddNoise() 
                rset.write_list(path,addNoise.list_c)
                if len(addNoise.list_c)>5:
                    addNoise.add_n(addNoise.list_c) 
                    judge_threshold(addNoise,rset,tree_class,"BB") 
                    judge_threshold(addNoise,rset,tree_class,"CSA")
                    tree_class.fill_verctor(rset,addNoise) 
                    if (addNoise.CFD_time_r["BB"]>0 or addNoise.CFD_time_r["CSA"]>0):
                        rset.i=save_waveform_threshold(out_file,rset.i,addNoise)        
                    tree_class.tree_out.Fill()
                    tree_class.init_parameter()
            else:
                break


if __name__ == '__main__':
    args = sys.argv[1:]
    input_file=args[0]
    o_ls=input_file.split("/")[:]
    # Outfilename and init_parameter
    rset = Settting()
    out_file=o_ls[0]+"/"+o_ls[1]+"/"+"outfile:"+o_ls[2]+"/"
    rset.create_outpath(out_file)
    rset.loop_out_p()
    # Root defined
    out_root_f=ROOT.TFile(out_file+"out.root","RECREATE")
    tree_class=RootFile()
    tree_class.root_define()
    # Add noise
    loop_addNoise(input_file,rset,tree_class,out_file)
    # Draw get time resolution for constant CFD
    draw_2D_CFD_time(rset.CFD_BB_time,out_file,"BB")
    draw_2D_CFD_time(rset.CFD_CSA_time,out_file,"CSA")
    tree_class.tree_out.Write()
    out_root_f.Close()  