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


def create_path(path):
    """ If the path does not exit, create the path"""
    if not os.access(path, os.F_OK):
        os.makedirs(path) 


'''
Draw 2D functions 
'''
def draw_2Dcurrent(det,fen,calcurrent,ele_current):

    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "canvas",1500,500)
    c.Divide(3,1)


    det.sum_cu.SetLineColor(1)
    det.positive_cu.SetLineColor(2)
    det.negtive_cu.SetLineColor(4)

    det.gain_positive_cu.SetLineColor(2)
    det.gain_negtive_cu.SetLineColor(4)

    det.gain_positive_cu.SetLineStyle(3)
    det.gain_negtive_cu.SetLineStyle(3)

    ele_current.BB_ele.SetLineColor(2)
    ele_current.CSA_ele.SetLineColor(4)

    det.negtive_cu.GetXaxis().SetTitleSize(0.05)
    det.negtive_cu.GetXaxis().SetLabelSize(0.04)
    det.negtive_cu.GetXaxis().SetNdivisions(510)
    # det.negtive_cu.GetYaxis().SetTitleOffset(1.1)
    det.negtive_cu.GetYaxis().SetTitleSize(0.05)
    det.negtive_cu.GetYaxis().SetLabelSize(0.04)
    det.negtive_cu.GetYaxis().SetNdivisions(505)
    det.negtive_cu.GetXaxis().CenterTitle()
    det.negtive_cu.GetYaxis().CenterTitle()
    det.negtive_cu.GetXaxis().SetTitle("Time [s]")
    det.negtive_cu.GetYaxis().SetTitle("Current [A]")
    det.negtive_cu.GetYaxis().SetRangeUser(0,12e-6)

    det.gain_positive_cu.GetXaxis().SetTitleSize(0.05)
    det.gain_positive_cu.GetXaxis().SetLabelSize(0.04)
    det.gain_positive_cu.GetXaxis().SetNdivisions(510)
    # det.gain_positive_cu.GetYaxis().SetTitleOffset(1.1)
    det.gain_positive_cu.GetYaxis().SetTitleSize(0.05)
    det.gain_positive_cu.GetYaxis().SetLabelSize(0.04)
    det.gain_positive_cu.GetYaxis().SetNdivisions(505)
    det.gain_positive_cu.GetXaxis().CenterTitle()
    det.gain_positive_cu.GetYaxis().CenterTitle()
    det.gain_positive_cu.GetXaxis().SetTitle("Time [s]")
    det.gain_positive_cu.GetYaxis().SetTitle("Current [A]")
    det.gain_positive_cu.GetYaxis().SetRangeUser(0,12e-6)

    det.sum_cu.GetXaxis().SetTitleSize(0.05)
    det.sum_cu.GetXaxis().SetLabelSize(0.04)
    det.sum_cu.GetXaxis().SetNdivisions(510)
    # det.sum_cu.GetYaxis().SetTitleOffset(1.1)
    det.sum_cu.GetYaxis().SetTitleSize(0.05)
    det.sum_cu.GetYaxis().SetLabelSize(0.04)
    det.sum_cu.GetYaxis().SetNdivisions(505)
    det.sum_cu.GetXaxis().CenterTitle()
    det.sum_cu.GetYaxis().CenterTitle()
    det.sum_cu.GetXaxis().SetTitle("Time [s]")
    det.sum_cu.GetYaxis().SetTitle("Current [A]")
    det.sum_cu.GetYaxis().SetRangeUser(0,12e-6)

    legend1 = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend1.AddEntry(det.negtive_cu, "initial electron", "l")
    legend1.AddEntry(det.positive_cu, "initial hole", "l")
    legend1.AddEntry(det.gain_negtive_cu, "gain electron", "l")
    legend1.AddEntry(det.gain_positive_cu, "gain hole", "l")

    legend2 = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend2.AddEntry(det.sum_cu, "e+h", "l")

    legend3 = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend3.AddEntry(ele_current.BB_ele, "BB Amp", "l")
    legend3.AddEntry(ele_current.CSA_ele, "CSA Amp", "l")

    legend1.SetBorderSize(0)
    legend1.SetTextFont(43)
    legend1.SetTextSize(20)

    legend2.SetBorderSize(0)
    legend2.SetTextFont(43)
    legend2.SetTextSize(20)

    legend3.SetBorderSize(0)
    legend3.SetTextFont(43)
    legend3.SetTextSize(20)

    c.cd(1)
    c.GetPad(1).SetLeftMargin(0.12)
    c.GetPad(1).SetRightMargin(0.12)
    # c.GetPad(1).SetTopMargin(0.12)
    c.GetPad(1).SetBottomMargin(0.14)
    # det.sum_cu.GetXaxis().SetTitleOffset(1.2)

    det.gain_positive_cu.Draw("HIST")
    det.gain_negtive_cu.Draw("SAME HIST")
    det.negtive_cu.Draw("SAME HIST")
    det.positive_cu.Draw("SAME HIST")
    legend1.Draw("SAME")
    c.Update()


    c.cd(2)
    c.GetPad(2).SetLeftMargin(0.12)
    c.GetPad(2).SetRightMargin(0.12)
    # c.GetPad(2).SetTopMargin(0.12)
    c.GetPad(2).SetBottomMargin(0.14)
    # det.sum_cu.GetXaxis().SetTitleOffset(1.2)
    det.sum_cu.Draw("HIST")
    legend2.Draw("SAME")
    c.Update()
    
    #ele_current.Draw("HIST")

    ele_current.BB_ele.GetXaxis().SetTitleSize(0.05)
    ele_current.BB_ele.GetXaxis().SetLabelSize(0.04)
    ele_current.BB_ele.GetXaxis().SetNdivisions(510)
    # ele_current.BB_ele.GetYaxis().SetTitleOffset(1.1)
    ele_current.BB_ele.GetYaxis().SetTitleSize(0.05)
    ele_current.BB_ele.GetYaxis().SetLabelSize(0.04)
    ele_current.BB_ele.GetYaxis().SetNdivisions(505)
    ele_current.BB_ele.GetXaxis().CenterTitle()
    ele_current.BB_ele.GetYaxis().CenterTitle()
    ele_current.BB_ele.GetXaxis().SetTitle("Time [s]")
    ele_current.BB_ele.GetYaxis().SetTitle("Current [A]")
    #ele_current.BB_ele.GetYaxis().SetRangeUser(0,120)
 
    c.cd(3)
    c.GetPad(3).SetLeftMargin(0.12)
    c.GetPad(3).SetRightMargin(0.12)
    # c.GetPad(4).SetTopMargin(0.12)
    c.GetPad(3).SetBottomMargin(0.14)
    ele_current.CSA_ele.Draw("HIST")
    ele_current.BB_ele.Draw("SAME HIST")

    legend3.Draw("SAME")
    c.Update()

    create_path("./fig")
    c.SaveAs("./fig/silicon_lgad_2D_drift_current_contribution_150V.pdf")

    charge_t=det.sum_cu.Integral() \
        * ((det.sum_cu.GetXaxis().GetXmax() \
        - det.sum_cu.GetXaxis().GetXmin()) \
        / det.sum_cu.GetNbinsX()) * 1e15
    
    print(charge_t)
    # print(qtot*1e15)
    calcurrent.draw_drift_path(det)
    fen.draw()



def draw_2D_gain_current_contribution(det):

    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c1", "canvas",500,500)
    c.Divide(3,1)


    det.gain_n_n_cu.SetLineColor(4)
    det.gain_n_p_cu.SetLineColor(2)
    det.gain_p_n_cu.SetLineColor(4)
    det.gain_p_p_cu.SetLineColor(2)


    det.gain_n_p_cu.SetLineStyle(3)
    det.gain_n_n_cu.SetLineStyle(3) 

    det.gain_n_p_cu.GetXaxis().SetTitleSize(0.05)
    det.gain_n_p_cu.GetXaxis().SetLabelSize(0.04)
    det.gain_n_p_cu.GetXaxis().SetNdivisions(510)
    # det.gain_n_p_cu.GetYaxis().SetTitleOffset(1.1)
    det.gain_n_p_cu.GetYaxis().SetTitleSize(0.05)
    det.gain_n_p_cu.GetYaxis().SetLabelSize(0.04)
    det.gain_n_p_cu.GetYaxis().SetNdivisions(505)
    det.gain_n_p_cu.GetXaxis().CenterTitle()
    det.gain_n_p_cu.GetYaxis().CenterTitle()
    det.gain_n_p_cu.GetXaxis().SetTitle("Time [s]")
    det.gain_n_p_cu.GetYaxis().SetTitle("Current [A]")
    det.gain_n_p_cu.GetYaxis().SetRangeUser(0,12e-6)


    c.cd()
    c.GetPad(0).SetLeftMargin(0.12)
    c.GetPad(0).SetRightMargin(0.12)
    # c.GetPad(0).SetTopMargin(0.12)
    c.GetPad(0).SetBottomMargin(0.14)


    det.gain_n_p_cu.Draw("HIST")
    det.gain_p_p_cu.Draw("SAME HIST")
    det.gain_n_n_cu.Draw("SAME HIST")
    det.gain_p_n_cu.Draw("SAME HIST")
    c.Update()

# 
    legend = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend.AddEntry(det.gain_n_p_cu, "e_h current", "l")
    legend.AddEntry(det.gain_p_p_cu, "h_h current", "l")
    legend.AddEntry(det.gain_n_n_cu, "e_e current", "l")
    legend.AddEntry(det.gain_p_n_cu, "h_e current", "l")
    # legend.AddEntry(det.sum_cu, "e+h", "l")
    # legend.AddEntry(ele_current, "electronics", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(20)

    c.cd()
    legend.Draw("same")

    c.SaveAs("silicon_lgad_2D_drift_current_contribution_150V.pdf")
