#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
@Description: Draw and plot drift path and induced current       
@Date       : 2021/08/31 11:09:40
@Author     : tanyuhang
@version    : 1.0
'''

import math
import ROOT
import sys
import os

def drawplot(my_d,ele_current,my_f,my_g4p,my_current):
    """
    @description:
        Draw electric field ,drift path and energy deposition
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    create_path("fig/")
    draw_ele_field(my_d,my_f,"xz",my_d.det_model,my_d.l_y*0.5) 
    draw_ele_field(my_d,my_f,"xy",my_d.det_model,my_d.l_z*0.5)
    draw_ele_field(my_d,my_f,"yz",my_d.det_model,my_d.l_x*0.5)
    draw_plot(my_d,ele_current.BB_ele) # Draw current
    #energy_deposition(my_g4p)   # Draw Geant4 depostion distribution
    #Draw Drift path
    #my_current.draw_drift_path(my_d,my_f,my_d.det_model)
     
def draw_unittest(my_d,ele_current,my_f,my_g4p,my_current):
    """
    @description:
        Draw electric field ,drift path and energy deposition
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    create_path("fig/")
    draw_plot(my_d,ele_current.BB_ele,unit_test=True) # Draw current

def savedata(my_d,output,batch_number,ele_current):
    if "scan" in my_d.det_model:
        output_path = (output + "_d:"+str(my_d.d_neff) 
        +"_v:"+str(my_d.v_voltage)+"_g:"+str(my_d.e_gap))
        create_path(output_path)
        ele_current.save_ele(batch_number,output_path)
        ele_current.save_charge(output_path)

def draw_ele_field(my_d,my_f,plane,sensor_model,depth):
    """
    @description:
        Draw eletric field
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    c1 = ROOT.TCanvas("c", "canvas",1000, 1000)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    # ROOT.gStyle.SetOptFit()
    c1.SetLeftMargin(0.12)
    c1.SetRightMargin(0.2)
    c1.SetBottomMargin(0.14)
    c1.SetRightMargin(0.12)
    c1.Divide(2,2)
    model = ["E","P","WP"]
    for i in range(1,4):
        c1.GetPad(i).SetRightMargin(0.2)
        c1.cd(i)
        e_field=fill_his(model[i-1],depth,my_d,my_f,plane,sensor_model)
        e_field.Draw("COLZ")
        c1.Update()
    c1.SaveAs("fig/ele_field"+plane+str(depth)+".root")
    del c1


def fill_his(model,depth,my_d,my_f,plane,sensor_model):
    """
    @description:
        Draw eletric field - Fill histrogram
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    nx_e=200
    ny_e=200
    d_r=confirm_range(my_d,my_f,plane,sensor_model,depth)
    e_v = ROOT.TH2F("","",nx_e,d_r[0],d_r[1],ny_e,d_r[2],d_r[3])
    for j in range (ny_e):
        for i in range(nx_e):
            x_v = (i+1)*((d_r[1]-d_r[0])/nx_e)+d_r[0]
            y_v = (j+1)*((d_r[3]-d_r[2])/ny_e)+d_r[2]
            f_v=0.0
            try:
                f_v,e_v = get_f_v(x_v,y_v,depth,model,my_f,plane,e_v,d_r)
                if model == "E":
                    f_v = math.sqrt(math.pow(f_v[0],2)
                                    +math.pow(f_v[1],2)
                                    +math.pow(f_v[2],2))                           
            except RuntimeError:
                f_v = 0.0
            e_v.SetBinContent(i+1,j+1,f_v)
    if plane == "xy":
        e_v.GetXaxis().SetTitle("x")
        e_v.GetYaxis().SetTitle("y")
    elif plane == "yz":
        e_v.GetXaxis().SetTitle("y")
        e_v.GetYaxis().SetTitle("z")
    elif plane == "xz":
        e_v.GetXaxis().SetTitle("x")
        e_v.GetYaxis().SetTitle("z") 
    return e_v


def get_f_v(i_x,i_y,i_z,model,my_f,plane,e_v,d_r):
    """
    @description:
        Draw eletric field - Get parameters
    @param:
        "E" -- electric
        "P" -- potential
        "WP" -- weigthing potential    
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    if plane == "xy":
        input_x=i_x
        input_y=i_y
        input_z=i_z
    elif plane == "yz":
        input_x=i_z
        input_y=i_x
        input_z=i_y
    elif plane == "xz":
        input_x=i_x
        input_y=i_z
        input_z=i_y
    if model == "E":
        e_v.SetTitle("electric field "+d_r[4])
        f_v=my_f.get_e_field(input_x,input_y,input_z)
    elif model == "P":
        e_v.SetTitle("potential "+d_r[4])
        f_v=my_f.get_potential(input_x,input_y,input_z)
    elif model =="WP":
        e_v.SetTitle("weigthing potential "+d_r[4]) 
        f_v=my_f.get_w_p(input_x,input_y,input_z)
    return f_v,e_v


def confirm_range(my_d,my_f,plane,sensor_model,depth):
    """
    @description:
        Draw eletric field - Confirm draw electric field detector range
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    if "plugin3D" in sensor_model:
        l_xl = my_f.sx_l
        l_xr = my_f.sx_r 
        if plane == "xy":
            l_yl = my_f.sy_l
            l_yr = my_f.sy_r
        elif plane == "yz" or plane == "xz":
            l_yl = 0
            l_yr = my_d.l_z
        else:
            print("the draw plane is not existing")
    elif "planar3D" in sensor_model:
        l_xl = 0
        l_yl = 0 
        if plane == "xy":
            l_xr = my_d.l_x 
            l_yr = my_d.l_y
        elif plane == "yz":
            l_xr = my_d.l_y
            l_yr = my_d.l_z
        elif plane == "xz":
            l_xr = my_d.l_x
            l_yr = my_d.l_z
        else:
            print("the draw plane is not existing")
    else:
        print("sensor model is wrrong")
        sys.exit()
    for x in "xyz":
        if x not in plane:
            t_name = plane + " at" + x + " = " + str(depth)
    return [l_xl,l_xr,l_yl,l_yr,t_name]


def draw_plot(my_detector, ele_current, unit_test=False):
    """
    @description:
        Save current in root file
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    c=ROOT.TCanvas("c","canvas1",1000,1000)
    c.cd()
    c.Update()
    c.SetLeftMargin(0.12)
    # c.SetTopMargin(0.12)
    c.SetBottomMargin(0.14)
    ROOT.gStyle.SetOptStat(0)
    current_canvas_setting(my_detector)
    rightmax=pad_setting(ele_current)
    axis_leng_setting(my_detector,ele_current,rightmax)
    c.Update()
    if unit_test: 
        c.SaveAs("pyraser/unittext/test.pdf")
    else:
        c.SaveAs("fig/basic_infor.root")
    del c


def current_canvas_setting(my_d):
    """
    @description:
        Current graph setting
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    my_d.sum_cu.GetXaxis().SetTitleOffset(1.2)
    my_d.sum_cu.GetXaxis().SetTitleSize(0.05)
    my_d.sum_cu.GetXaxis().SetLabelSize(0.04)
    my_d.sum_cu.GetXaxis().SetNdivisions(510)
    my_d.sum_cu.GetYaxis().SetTitleOffset(1.1)
    my_d.sum_cu.GetYaxis().SetTitleSize(0.05)
    my_d.sum_cu.GetYaxis().SetLabelSize(0.04)
    my_d.sum_cu.GetYaxis().SetNdivisions(505)
    my_d.sum_cu.GetXaxis().CenterTitle()
    my_d.sum_cu.GetYaxis().CenterTitle() 
    my_d.sum_cu.GetXaxis().SetTitle("Time [s]")
    my_d.sum_cu.GetYaxis().SetTitle("Current [A]")
    my_d.sum_cu.Draw("HIST")
    my_d.positive_cu.Draw("SAME HIST")
    my_d.negtive_cu.Draw("SAME HIST")
    my_d.sum_cu.SetLineColor(3)
    my_d.positive_cu.SetLineColor(2)
    my_d.negtive_cu.SetLineColor(4)
    my_d.sum_cu.SetLineWidth(2)
    my_d.positive_cu.SetLineWidth(2)
    my_d.negtive_cu.SetLineWidth(2)


def pad_setting(ele_current):
    """
    @description:
        Canvas pad setting
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    if ele_current.GetMinimum() < 0:
        rightmax = 1.1*ele_current.GetMinimum()
    else:
        rightmax = 1.1*ele_current.GetMaximum()
    if rightmax == 0:
        n_scale=0
    elif ele_current.GetMinimum() <0:
        n_scale = ROOT.gPad.GetUymin() / rightmax
    else:
        n_scale = ROOT.gPad.GetUymax() / rightmax
    ele_current.Scale(n_scale)
    ele_current.Draw("SAME HIST")
    ele_current.SetLineWidth(2)   
    ele_current.SetLineColor(6)
    return rightmax


def axis_leng_setting(my_detector,ele_current,rightmax):
    """
    @description:
        Canvas axis and legend setting
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), 
                       ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), 
                       rightmax, 0, 510, "+L")
    axis.SetLineColor(6)
    axis.SetTextColor(6)
    axis.SetTextSize(0.02)
    axis.SetLabelColor(6)
    axis.SetLabelSize(0.02)
    axis.SetTitle("Ampl [mV]")
    axis.CenterTitle()
    axis.Draw("same")

    legend = ROOT.TLegend(0.5, 0.3, 0.9, 0.6)
    legend.AddEntry(my_detector.negtive_cu, "electron", "l")
    legend.AddEntry(my_detector.positive_cu, "hole", "l")
    legend.AddEntry(my_detector.sum_cu, "e+h", "l")
    legend.AddEntry(ele_current, "electronics", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(45)
    legend.Draw("same")


def energy_deposition(my_g4v):
    """
    @description:
        Energy_deposition for multi events of Geant4 simulation
    @param:
        None     
    @Returns:
        None
    @Modify:
        2021/08/31
    """
    c1=ROOT.TCanvas("c1","canvas1",1000,1000)
    h1 = ROOT.TH1F("Edep_device", "Energy deposition in SiC", 100, 0., 0.1)
    for i in range (len(my_g4v.edep_devices)):
        h1.Fill(my_g4v.edep_devices[i])
    g1 = ROOT.TF1("m1","landau",0,0.1)
    h1.Fit(g1,"S")
    print("MPV:%s"%g1.GetParameter(1))
    h1.Draw()
    c1.SaveAs("fig/dep_SiC_energy.root")

def create_path(path):
    """ If the path does not exit, create the path"""
    if not os.access(path, os.F_OK):
        os.makedirs(path) 