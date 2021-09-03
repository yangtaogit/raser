#!/usr/bin/env python3
'''
author: tanyuhang
time: 2021.3.8
Use: Read the data of KDetsim induced current
'''
import os
import sys
import re
from ROOT import TCanvas,TGraph,TMultiGraph,TLegend,gRandom,TF1,TLatex
# from tools import set_root_style
from array import array
import ROOT
from ROOT import gStyle
import math

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
h_CSA_time_resolution=array('f',[0])\

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


def main(): 
    args = sys.argv[1:]
    input_file=args[0]
    out_p=args[1]
    thre_vth=18 #mv
    CFD=0.5
    # root defined
    out_file=ROOT.TFile(out_p,"RECREATE")
    tree_out=ROOT.TTree('tree','tree')

    tree_out.Branch('Events',Events,'Events/I')
    tree_out.Branch('h_pulse_time',h_pulse_time)
    tree_out.Branch('h_BB_pulse_height',h_BB_pulse_height) 
    tree_out.Branch('h_BB_nps_height',h_BB_nps_height)  
    tree_out.Branch('h_BB_max_nps_height',h_BB_max_nps_height,'h_BB_max_nps_height/F')    
    tree_out.Branch('h_BB_max_pulse_time',h_BB_max_pulse_time,'h_BB_max_pulse_time/F')
    tree_out.Branch('h_BB_time_resolution',h_BB_time_resolution,'h_BB_time_resolution/F')

    tree_out.Branch('h_CSA_pulse_height',h_CSA_pulse_height) 
    tree_out.Branch('h_CSA_nps_height',h_CSA_nps_height) 
    tree_out.Branch('h_CSA_max_nps_height',h_CSA_max_nps_height,'h_CSA_max_nps_height/F')    
    tree_out.Branch('h_CSA_max_pulse_time',h_CSA_max_pulse_time,'h_CSA_max_pulse_time/F')
    tree_out.Branch('h_CSA_time_resolution',h_CSA_time_resolution,'h_CSA_time_resolution/F')

    tree_out.Branch('h_CSA_dVdt',h_CSA_dVdt,'h_CSA_dVdt/F')
    tree_out.Branch('h_BB_dVdt',h_BB_dVdt,'h_BB_dVdt/F')

    tree_out.Branch('h_noise_height',h_noise_height)
    tree_out.Branch('h_noise_CSA_height_jitter',h_noise_CSA_height_jitter,'h_noise_CSA_height_jitter/F')
    tree_out.Branch('h_noise_BB_height_jitter',h_noise_BB_height_jitter,'h_noise_BB_height_jitter/F')
    tree_out.Branch('h_noise_height_RMS',h_noise_height_RMS,'h_noise_height_RMS/F')
    #start_write_to_txt(out_path) 
    i=1
    CFD_BB_time=[]
    CFD_CSA_time=[]
    for root,dirs,files in os.walk(input_file):
        for file in files:    
            if i<100000:
                print(file)
                print("................events:%s..............."%(Events[0])) 
                print("................Save events:%s..............."%i)
                Events[0]+=1
                list_c=[]
                path = os.path.join(input_file, file)
                BB_Fv = 0
                CSA_Fv = 0
                BB_time_r = 0
                CSA_time_r = 0
                BB_noise_height_jitter = 0
                CSA_noise_height_jitter = 0

                for line in open(path):
                    if not (is_number(line.split(",")[0])):
                        continue
                    list_c.append(line)
                if len(list_c)>5:
                    time_list,noise_height_list,CSA_paras,BB_paras=add_noise(list_c)
                    noise_height_RMS= math.sqrt(sum([x**2 for x in noise_height_list])/len(noise_height_list))
                    h_noise_height_RMS[0]=noise_height_RMS
                    #BB amplifer
                    if (BB_paras["max_BB_nps_height"]>thre_vth and BB_paras["max_BB_nps_height"]<80):

                        h_BB_max_nps_height[0]=BB_paras["max_BB_nps_height"]
                        h_BB_max_pulse_time[0]=BB_paras["max_BB_pulse_time"]
                        BB_time_r,BB_noise_height_jitter,BB_dVdt=get_CFD_time(time_list,BB_paras,CFD,noise_height_RMS,"BB")
                        # if BB_noise_height_jitter > 0 and  BB_noise_height_jitter < 1:
                        CFD_BB_time.append(BB_time_r)                     
                        h_BB_time_resolution[0]=BB_time_r
                        h_noise_BB_height_jitter[0]=BB_noise_height_jitter
                        h_BB_dVdt[0]=BB_dVdt
                        BB_Fv=1
                    #CSA amplifer
                    if (CSA_paras["max_CSA_nps_height"]>thre_vth and CSA_paras["max_CSA_nps_height"]<80):
    
                        h_CSA_max_nps_height[0]=CSA_paras["max_CSA_nps_height"]
                        h_CSA_max_pulse_time[0]=CSA_paras["max_CSA_pulse_time"]
                        CSA_time_r,CSA_noise_height_jitter,CSA_dVdt=get_CFD_time(time_list,CSA_paras,CFD,noise_height_RMS,"CSA")
                        h_noise_CSA_height_jitter[0]=CSA_noise_height_jitter
                        h_CSA_dVdt[0]=CSA_dVdt
                        # if CSA_noise_height_jitter > 0 and  CSA_noise_height_jitter < 1:
                        CFD_CSA_time.append(CSA_time_r)                     
                        h_CSA_time_resolution[0]=CSA_time_r
                        CSA_Fv=1
                    
                    if (BB_Fv==1 or CSA_Fv==1):
                        for j in range(0,len(time_list)):
                            h_pulse_time.push_back(time_list[j])
                            h_BB_pulse_height.push_back(BB_paras["ampl_BB_s_list"][j])
                            h_CSA_pulse_height.push_back(CSA_paras["ampl_CSA_s_list"][j])
                            h_BB_nps_height.push_back(BB_paras["ampl_BB_nps_list"][j])
                            h_CSA_nps_height.push_back(CSA_paras["ampl_CSA_nps_list"][j])
                            h_noise_height.push_back(noise_height_list[j])
                        print(CSA_time_r)
                        print(BB_time_r)
                        if (CSA_time_r>0 or BB_time_r>0):
                            i=save_waveform_threshold(input_file,i,BB_paras,CSA_paras)
                        
                        tree_out.Fill()
                        h_BB_max_nps_height[0]=0.0
                        h_BB_max_pulse_time[0]=0.0  
                        h_BB_time_resolution[0]=0.0
                        h_CSA_max_nps_height[0]=0.0
                        h_CSA_max_pulse_time[0]=0.0  
                        h_CSA_time_resolution[0]=0.0
                        h_noise_BB_height_jitter[0]=0.0
                        h_noise_CSA_height_jitter[0]=0.
                        h_BB_dVdt[0]=0.0
                        h_CSA_dVdt[0]=0.0

                        h_pulse_time.clear()
                        h_BB_nps_height.clear()
                        h_CSA_nps_height.clear()	
                        h_BB_pulse_height.clear()
                        h_CSA_pulse_height.clear()	
                        h_noise_height.clear()
                else:
                    break
    tree_out.Write()
    out_file.Close()
    # draw_toa(bc_CFD_times,file_n,input_file,out_p)
    #draw_CFD_time(CFD_time,input_file,out_p,file_m)
    draw_2D_CFD_time(CFD_BB_time,input_file,"BB")
    draw_2D_CFD_time(CFD_CSA_time,input_file,"CSA")
            # i+=1
    # draw_mg(mg,leg,out_p,"BB")
    # draw_mg(mg1,leg,out_p,"CSA")
    # draw_mg(mg2,leg,out_p,"Itotal")
def draw_toa(CFD_time,input_n,input_file,out_p):
    output_file = input_file.split("/")[2]
    out_name = out_p.split("/")[1]
    c1 = TCanvas("c1", "canvas1", 800, 600)
    #gStyle.SetPalette(57)
    gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(1)
    #set_root_style(stat=0, grid=0)
    c1.SetGrid()
    c1.SetLeftMargin(0.12)
    # c.SetTopMargin(0.12)
    c1.SetBottomMargin(0.14)
    histo1=ROOT.TH2F("","",100,0,100,100,0,100)
    histo1.GetXaxis().SetTitleOffset(1.2)
    histo1.GetXaxis().SetTitleSize(0.05)
    histo1.GetXaxis().SetLabelSize(0.05)
    histo1.GetXaxis().SetNdivisions(510)
    histo1.GetYaxis().SetTitleOffset(1.1)
    histo1.GetYaxis().SetTitleSize(0.05)
    histo1.GetYaxis().SetLabelSize(0.05)
    histo1.GetYaxis().SetNdivisions(505)
    histo1.GetXaxis().SetTitle("x [#mum]");
    histo1.GetYaxis().SetTitle("y [#mum]")
    histo1.GetXaxis().CenterTitle()
    histo1.GetYaxis().CenterTitle()
    gap_v=float(output_file.split('_')[1])
    voltage_v=float(output_file.split('_')[3])
    # e_r=5.0
    # l_x = 100
    # l_y = 100
    # e_t_y = math.sqrt(gap_v*gap_v*0.75)
    # x_min = l_x*0.5 - gap_v*0.5
    # x_max = l_x*0.5 + gap_v*0.5
    # y_min = l_y*0.5 - e_t_y
    # y_max = l_y*0.5 + e_t_y
    for i in range(len(CFD_time)):
        n_x = float(input_n[i].split("_")[3])
        n_y = float(input_n[i].split("_")[5])
        r_v = math.sqrt((n_x-50)*(n_x-50) + (n_y-50)*(n_y-50))
        # if CFD_time[i]*1e9>=0 and CFD_time[i]*1e9<=2:
        #     time = CFD_time[i]*1e9
        # else:
        #     time = 0.0*1e9
        if CFD_time[i]!=0 and r_v <= gap_v: 
            histo1.Fill(n_x,n_y,CFD_time[i]*1e9)

    ROOT.gStyle.SetPalette(107)
    histo1.Draw("COLZ")
    # histo1.SetMaximum(1.5)
    # histo1.SetMinimum(0)
    c1.Update()
    c1.SaveAs("out/"+out_name+"/toa_"+output_file+".pdf")
    c1.SaveAs("out/"+out_name+"/toa_"+output_file+".C")
# def toa_bc(list_c,CFD):
#     bc_time_list = []
#     bc_pulse_height = []
#     for j in range (0,len(list_c)):
#         time= float(list(filter(None,list_c[j].split(",")))[0])
#         ampl_CSA=float(list(filter(None,list_c[j].split(",")))[1])
#         bc_time_list.append(time)
#         bc_pulse_height.append(ampl_CSA)
#     max_pulse_height=max(bc_pulse_height)
#     max_index=bc_pulse_height.index(max(bc_pulse_height))
#     max_pulse_time=bc_time_list[max_index]
#     bc_CFD_time = 0
#     for i in range (0,len(bc_time_list)):
#         if bc_pulse_height[i]>=max_pulse_height*CFD and bc_time_list[i]<max_pulse_time:
#             bc_CFD_time = bc_time_list[i]
#             break
#     return bc_CFD_time

def add_noise(list_c):
    time_list=[]

    CSA_paras = {}
    BB_paras = {}
    ampl_CSA_nps_list=[]
    ampl_BB_nps_list=[]
    ampl_CSA_s_list=[]
    ampl_BB_s_list=[]
    noise_height_list=[]
    time=0.0
    ampl_CSA_nps=0.0
    ampl_CSA_s=0.0
    ampl_BB_nps=0.0
    ampl_BB_s=0.0
    #noise
    gRandom.SetSeed(0)
    random_gauss = gRandom.Gaus
    for j in range (0,len(list_c)):
        time= float(list(filter(None,list_c[j].split(",")))[0])
        noise_height=random_gauss(-0.133,2.671)
        #nps noise plus signal
        ampl_CSA_nps=-float(list(filter(None,list_c[j].split(",")))[1])+noise_height
        ampl_BB_nps=-float(list(filter(None,list_c[j].split(",")))[2])+noise_height
        ampl_CSA_s=-float(list(filter(None,list_c[j].split(",")))[1])
        ampl_BB_s=-float(list(filter(None,list_c[j].split(",")))[2])
        time_list.append(time)
        noise_height_list.append(noise_height)
        ampl_CSA_nps_list.append(ampl_CSA_nps)
        ampl_BB_nps_list.append(ampl_BB_nps)
        ampl_CSA_s_list.append(ampl_CSA_s)
        ampl_BB_s_list.append(ampl_BB_s)

    max_CSA_nps_height=max(ampl_CSA_nps_list)
    max_CSA_index=ampl_CSA_nps_list.index(max(ampl_CSA_nps_list))
    max_CSA_pulse_time=time_list[max_CSA_index]

    max_CSA_s_height=max(ampl_CSA_s_list)
    max_CSA_s_index=ampl_CSA_s_list.index(max(ampl_CSA_s_list))
    max_CSA_s_time=time_list[max_CSA_s_index]

    max_BB_nps_height=max(ampl_BB_nps_list)
    max_BB_index=ampl_BB_nps_list.index(max(ampl_BB_nps_list))
    max_BB_pulse_time=time_list[max_BB_index]

    max_BB_s_height=max(ampl_BB_s_list)
    max_BB_s_index=ampl_BB_s_list.index(max(ampl_BB_s_list))
    max_BB_s_time=time_list[max_BB_s_index]

    CSA_paras["max_CSA_nps_height"] = max_CSA_nps_height
    CSA_paras["max_CSA_pulse_time"] = max_CSA_pulse_time
    CSA_paras["ampl_CSA_nps_list"] = ampl_CSA_nps_list
    CSA_paras["ampl_CSA_s_list"] = ampl_CSA_s_list
    CSA_paras["max_CSA_s_height"] = max_CSA_s_height
    CSA_paras["max_CSA_s_time"] = max_CSA_s_time

    BB_paras["max_BB_nps_height"] = max_BB_nps_height
    BB_paras["max_BB_pulse_time"] = max_BB_pulse_time
    BB_paras["ampl_BB_nps_list"] =  ampl_BB_nps_list
    BB_paras["ampl_BB_s_list"] =  ampl_BB_s_list
    BB_paras["max_BB_s_height"] = max_BB_s_height
    BB_paras["max_BB_s_time"] = max_BB_s_time
    BB_paras["time_list"] = time_list

    return time_list,noise_height_list,CSA_paras,BB_paras

def draw_CFD_time(CFD_time,out_put,out_p,file_m):
    print(out_put)
    out_put_1=out_put.split("/")[2]
    voltage_1 = out_put.split("/")[2]
    voltage_2 = float(voltage_1.split("_")[3])
    gap_v = float(voltage_1.split("_")[1])
    out_name = out_p.split("/")[1]
    c = ROOT.TCanvas("c", "canvas", 800, 600)
    #gStyle.SetPalette(57)
    # e_r=5.0
    # l_x = 100
    # l_y = 100
    # e_t_y = math.sqrt(gap_v*gap_v*0.75)
    # x_min = l_x*0.5 - gap_v*0.5
    # x_max = l_x*0.5 + gap_v*0.5
    # y_min = l_y*0.5 - e_t_y
    # y_max = l_y*0.5 + e_t_y
    gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(1)
    #set_root_style(stat=0, grid=0)
    c.SetGrid()
    c.SetLeftMargin(0.2)
    # c.SetTopMargin(0.12)
    c.SetBottomMargin(0.2)
    leg = ROOT.TLegend(0.25, 0.6, 0.40, 0.8)
    
    #FormatLegend(leg)

    histo=ROOT.TH1F("","",100,-0.2,0.6)
    gStyle.SetOptFit()
    #gStyle.SetStatY(0.6);
    for i in range(0,len(CFD_time)):
        n_x = float(file_m[i].split("_")[3])
        n_y = float(file_m[i].split("_")[5])
        r_v = math.sqrt((n_x-50)*(n_x-50) + (n_y-50)*(n_y-50))
        if CFD_time[i]!=0 and r_v <= gap_v:
            histo.Fill(CFD_time[i])
    histo.GetXaxis().SetTitle("ToA [ns]");
    histo.GetYaxis().SetTitle("Events")
    histo.GetXaxis().CenterTitle()
    histo.GetYaxis().CenterTitle()
    fit_func_1 = TF1('fit_func_1','gaus',4.8,5.5)
    histo.Fit("fit_func_1","ROQ+","",4.8,5.5)
    print("constant:%s"%fit_func_1.GetParameter(0))
    print("constant_error:%s"%fit_func_1.GetParError(0))
    print("mean:%s"%fit_func_1.GetParameter(1))
    print("mean_error:%s"%fit_func_1.GetParError(1))
    print("sigma:%s"%fit_func_1.GetParameter(2))
    print("sigma_error:%s"%fit_func_1.GetParError(2))
    sigma=fit_func_1.GetParameter(2)*1000
    error=fit_func_1.GetParError(2)*1000
    leg.AddEntry(fit_func_1,"Fit","L")
    leg.AddEntry(histo,"Sim","L")
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
    fit_func_1.SetLineWidth(2)
    histo.Draw()
    fit_func_1.Draw("same")
    leg.Draw("same")
    tex = ROOT.TLatex()
    tex.SetNDC(1)
    tex.SetTextFont(43)
    tex.SetTextSize(25)
    tex.DrawLatexNDC(0.6, 0.8, "Electron Spacing:%.0f#mum"%(gap_v))
    tex.DrawLatexNDC(0.6, 0.7, "Volatge=%.0fV,CFD=0.5"%(voltage_2))
    tex.DrawLatexNDC(0.6, 0.6, "#sigma = %.0f #pm %.0f ps"%(sigma,error))
    c.Update()
    c.SaveAs("out/"+out_name+"/"+out_name+".pdf")
    # c.SaveAs("out/"+out_name+"/"+out_put_1+".c")
    c.SaveAs("out/"+out_name+"/"+out_name+".root")

def draw_2D_CFD_time(CFD_time,out_put,model):

    # c1=ROOT.TCanvas()
    c1 = TCanvas("c1"+model,"c1"+model,200,10,800,600)

    gStyle.SetOptStat(0)
    ROOT.gROOT.SetBatch(1)
    #set_root_style(stat=0, grid=0)
    c1.SetGrid()
    c1.SetLeftMargin(0.2)
    # c.SetTopMargin(0.12)
    c1.SetBottomMargin(0.2)
    leg = ROOT.TLegend(0.25, 0.6, 0.40, 0.8)

    histo=ROOT.TH1F("","",20,4.6,5.6)
    gStyle.SetOptFit()
    #gStyle.SetStatY(0.6);
    for i in range(0,len(CFD_time)):
        if CFD_time[i]>0:
            histo.Fill(CFD_time[i])
    histo.GetXaxis().SetTitle("ToA [ns]")
    histo.GetYaxis().SetTitle("Events")
    histo.GetXaxis().CenterTitle()
    histo.GetYaxis().CenterTitle()
    fit_func_1 = TF1('fit_func_1','gaus',4.6,5.6)
    histo.Fit("fit_func_1","ROQ+","",4.6,5.6)
    print("constant:%s"%fit_func_1.GetParameter(0))
    print("constant_error:%s"%fit_func_1.GetParError(0))
    print("mean:%s"%fit_func_1.GetParameter(1))
    print("mean_error:%s"%fit_func_1.GetParError(1))
    print("sigma:%s"%fit_func_1.GetParameter(2))
    print("sigma_error:%s"%fit_func_1.GetParError(2))
    sigma=fit_func_1.GetParameter(2)*1000
    error=fit_func_1.GetParError(2)*1000
    leg.AddEntry(fit_func_1,"Fit","L")
    leg.AddEntry(histo,"Sim","L")
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
    fit_func_1.SetLineWidth(2)
    histo.Draw()
    fit_func_1.Draw("same")
    leg.Draw("same")
    tex = ROOT.TLatex()
    tex.SetNDC(1)
    tex.SetTextFont(43)
    tex.SetTextSize(25)
    tex.DrawLatexNDC(0.6, 0.7, "CFD=0.5"+" "+model+"ampl")
    tex.DrawLatexNDC(0.6, 0.6, "#sigma = %.0f #pm %.0f ps"%(sigma,error))
    c1.Update()
    c1.SaveAs(out_put+model+".pdf")

def get_CFD_time(time_list,Ampl_paras,CFD,noise_height_RMS,model):

    # perm_sic = 9.76  #Permittivity SiC
    # DCap = 5000*5000/100*perm_sic*8.851e-3 #fF backplane
    # DCap += 0.014*perm_sic*4*5000  #fF n layer
    # DCap +=50 #fF fixed 
    # noise_sen = 2.0*DCap/math.sqrt(1) #
    # print(noise_sen)
    random_gauss = gRandom.Gaus
    # noise_height_jitter=random_gauss(7.82,1.01)
    CFD_time=0.0
    jitter=0.0
    dVdt=0.0
    if model == "CSA":
        for i in range (0,len(time_list)):
            if Ampl_paras["ampl_CSA_nps_list"][i]>=Ampl_paras["max_CSA_nps_height"]*CFD and time_list[i]<Ampl_paras["max_CSA_pulse_time"] and time_list[i+2]<Ampl_paras["max_CSA_pulse_time"] and time_list[i-2]>1.0e-9:
                
                dVdt=(Ampl_paras["ampl_CSA_nps_list"][i+2]-Ampl_paras["ampl_CSA_nps_list"][i-2])/(time_list[i+2]-time_list[i-2])/1e9/1.38       
                if (dVdt!=0):

                    jitter=noise_height_RMS/dVdt # 1.54 is the correction of dV/dt, need to change from physics
                    CFD_time = 3.85+time_list[i]*1e9+random_gauss(0,jitter) #4.85 is the initial time by personal customization
                else:
                    # CFD_time=time_list[i]*1e9
                    CFD_time=0
        
                break
    elif model == "BB":
        for i in range (0,len(time_list)):
            if Ampl_paras["ampl_BB_nps_list"][i]>=Ampl_paras["max_BB_nps_height"]*CFD and time_list[i]<Ampl_paras["max_BB_pulse_time"] and time_list[i+2]<Ampl_paras["max_BB_pulse_time"] and time_list[i-2]>1.0e-9:
                
                dVdt=(Ampl_paras["ampl_BB_nps_list"][i+2]-Ampl_paras["ampl_BB_nps_list"][i-2])/(time_list[i+2]-time_list[i-2])/1e9/1.28

                if (dVdt!=0):              
                    jitter=noise_height_RMS/dVdt				
                    CFD_time = 3.85+time_list[i]*1e9+random_gauss(0,jitter)
                else:
                    # CFD_time=time_list[i]*1e9
                    CFD_time=0             
                break
    else:
        print("model is wrong")
        sys.exit()
    return CFD_time,jitter,dVdt

def get_CFD_s_time(time_list,Ampl_paras,CFD,noise_height_RMS,model):
    
    random_gauss = gRandom.Gaus
    # noise_height_jitter=random_gauss(7.82,1.01)
    CFD_time=0.0
    jitter=0.0
    if model == "CSA":
        for i in range (0,len(time_list)):
            if Ampl_paras["ampl_CSA_s_list"][i]>=Ampl_paras["max_CSA_s_height"]*CFD and time_list[i]<Ampl_paras["max_CSA_s_time"] and time_list[i+2]<Ampl_paras["max_CSA_s_time"] and time_list[i-2]>0:
                
                dVdt=(Ampl_paras["ampl_CSA_s_list"][i+2]-Ampl_paras["ampl_CSA_s_list"][i-2])/(time_list[i+2]-time_list[i-2])        
                if (dVdt!=0):

                    jitter=noise_height_RMS/dVdt*1e9
                    #CFD_time = time_list[i]*1e9+random_gauss(0,abs(jitter))
                    CFD_time = random_gauss(0,abs(jitter))
                else:
                    #CFD_time=time_list[i]*1e9
                    CFD_time=0            
                break
    elif model == "BB":
        for i in range (0,len(time_list)):
            if Ampl_paras["ampl_BB_s_list"][i]>=Ampl_paras["max_BB_s_height"]*CFD and time_list[i]<Ampl_paras["max_BB_s_time"] and time_list[i+2]<Ampl_paras["max_BB_s_time"] and time_list[i-2]>0:
                
                dVdt=(Ampl_paras["ampl_BB_s_list"][i+2]-Ampl_paras["ampl_BB_s_list"][i-2])/(time_list[i+2]-time_list[i-2]) 
                
                if (dVdt!=0):              
                    jitter=noise_height_RMS/dVdt*1e9				
                    #CFD_time = time_list[i]*1e9+random_gauss(0,abs(jitter))
                    CFD_time = random_gauss(0,abs(jitter))
                else:
                    #CFD_time=time_list[i]*1e9
                    CFD_time=0
                # CFD_time=time_list[i]*1e9                
                break    
    else:
        print("model is wrong")
        sys.exit()
    return CFD_time,jitter
def save_waveform_threshold(output_file,event_n,BB_paras,CSA_paras):
    output_path = output_file + "_txt/"
    if event_n ==1 :
        os.system("mkdir %s -p"%(output_path))
    f1 = open(output_path+"t_"+str(event_n)+".csv","w")
    f1.write("time[ns],CSA Amplitude [mV], BB Amplitude [mV] \n")
    for i in range(len(BB_paras["time_list"])):
        time = BB_paras["time_list"][i]
        BB = BB_paras["ampl_BB_nps_list"][i]
        CSA = CSA_paras["ampl_CSA_nps_list"][i]
        f1.write("%s,%s,%s \n"%(time,BB,CSA))
    f1.close()
    return event_n+1        

def FormatLegend(leg):
    
    leg.SetBorderSize(1)
    leg.SetTextFont(43)
    leg.SetTextSize(40)
    leg.SetFillStyle(1)
    leg.SetFillColor(1)
    leg.SetLineColor(2) 

def set_color_marker(color,marker,i,gr):
    f=marker[i]
    gr.SetMarkerStyle(f)
    gr.SetMarkerSize(1)
    k=color[i]
    gr.SetLineColor(k)
    gr.SetLineWidth(2)
    gr.SetMarkerColor(k)
    return gr

def fill_legend(leg,gr,name):
    # leg name define
    #leg_name="position:x"+str(x_number)+"_y"+str(y_number)
    # print c_number
    # if c_number == 0:
    #     leg_name="first_point"
    # if c_number == 1:
    #     leg_name="second_point"
    # if c_number == 2:
    #     leg_name="third_point"
    leg.AddEntry(gr,name,"LP")
    return leg

def defind_color_marker(marker,color):
    
    color.append(int(2))
    color.append(int(45))
    color.append(int(3))
    color.append(int(7))
    color.append(int(46))
    color.append(int(38))
    color.append(int(40))
    color.append(int(4))
    color.append(int(8))
    color.append(int(11))

    marker.append(int(20))
    marker.append(int(24))
    marker.append(int(23))
    marker.append(int(25))
    marker.append(int(33)) 
    marker.append(int(26))
    marker.append(int(30))
    marker.append(int(27))
    marker.append(int(28))
    marker.append(int(13))
    return marker,color
def is_number(s):
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
def FormatLegend(leg):
    
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(30)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetLineColor(0)    
if __name__ == '__main__':
    main()   