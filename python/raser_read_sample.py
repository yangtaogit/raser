'''
author: tanyuhang
time: 2021.3.8
Use: Read the data of KDetsim induced current
'''
import os
import sys
import re
from ROOT import TCanvas,TGraph,TMultiGraph,TLegend
# from tools import set_root_style
from array import array
import ROOT
from ROOT import gStyle

time_delay=0  #20ns
def main(): 
    args = sys.argv[1:]
    input_file=args[0]
    out_p=args[1]
    model=args[2]
    list = os.listdir(input_file) 
    mg = TMultiGraph("mg","")
    leg=TLegend(0.55,0.65,0.75,0.85)
    color = []
    marker = []
    #start_write_to_txt(out_path) 
    k=int(0) 
    c_number=0
    i=0
    marker,color=defind_color_marker(marker,color)
    for root,dirs,files in os.walk(input_file):
        for file in files:
            #if (i == 17 and j ==40) or (i == 34 and j ==40) or (i == 28 and j ==29):
            #data_name="x"+str(i)+"_y"+str(j)+".C"
            path = os.path.join(input_file, file) 
            if ".C" in path:
                if i <100000:
                    i+=1
                    with open(path) as file_obj:
                        out_name=out_p+file
                        out_path=out_name.split(".C")[0]+".txt"
                        print(out_path)
                        write_zero_to_txt(out_path,k)
                        k=k+1
                        pattern = re.findall(r'tent[(](.*?)[)]',file_obj.read())
                        if len(pattern)<2:
                            continue
                        time_list=[]
                        current_list=[]
                        time_list,current_list=save_txt(pattern,out_path,k,model)
                    write_end_to_txt(out_path,time_list)
                    # gr=save_graph(time_list,current_list)
                #gr=set_color_marker(color,marker,c_number,gr)              
                #leg=fill_legend(leg,gr,out_path)   
                    # c_number=c_number+1  
                    # mg.Add(gr)
    print("event number: %s"%i)
    #draw the waveform together
    
    # draw_mg(mg,leg,out_p,model)

def save_txt(pattern,name,k,model):
    #change .C file to txt file and multiply a factor
    time_factor=3e-12
    e_0 = 1.60217733e-19
    #current_factor=1e-6
    k=k-1
    f = open(name,"a")
    time_list=[]
    current_list=[]
    time=0.0
    current=0.0
    for j in range (0,len(pattern)):
        #time=float(pattern[j].split(',')[0])*time_factor+k*time_delay
        #current=float(pattern[j].split(',')[1])/time_factor*100*50*1000*1e-15 
        time=float(pattern[j].split(',')[0])*time_factor
        if model == "I":
            current=float(pattern[j].split(',')[1])*1e6
        if model == "V":
            current=float(pattern[j].split(',')[1])			
        #current=(charge/delta time)*10 10 is factor, because we only simulate 0.1 electron-hole pairs
        time_list.append(float(pattern[j].split(',')[0])*time_factor)
        current_list.append(abs(current))
        f.write(str(time)+","+str(abs(current))+"\n")
    f.close()
    return time_list,current_list

def start_write_to_txt(out_path):
    f = open(out_path,"w")
    f.close()
def write_zero_to_txt(out_path,k):
    f = open(out_path,"a")
    time=k*time_delay
    f.write(str(time)+",0.0\n")
    f.close()
def write_end_to_txt(out_path,time_list):
    f = open(out_path,"a")
    time=1.5*time_list[len(time_list)-1]
    f.write(str(time)+",0.0\n")
    f.close()
def save_graph(time_list,current_list):
    time_a=array( 'f' )
    current_a=array( 'f' )
    time_a.extend(time_list)
    current_a.extend(current_list)
    gr=TGraph(len(time_a),time_a,current_a)
    return gr

def draw_mg(mg,leg,out_path,model):
    
    c = TCanvas("c", "canvas", 800, 600)
    FormatLegend(leg)
    gStyle.SetPalette(55)
    c.SetLeftMargin(0.12)
    c.SetTopMargin(0.12)
    # set_root_style(stat=0, grid=0)
    out_name=out_path.split("/")[1]	
    if model == "I":
        mg.SetTitle("; time (s); current [#mu A]")
    if model == "V":
        mg.SetTitle("; time (s); Amplitude [mV]")
    mg.GetXaxis().CenterTitle()
    mg.GetYaxis().CenterTitle()
    mg.Draw("APL")
    leg.Draw("same")
    c.SaveAs("result/"+out_name+model+".pdf")

def FormatLegend(leg):
    
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(30)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetLineColor(0) 

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
def FormatLegend(leg):
    
    leg.SetBorderSize(0)
    leg.SetTextFont(43)
    leg.SetTextSize(30)
    leg.SetFillStyle(0)
    leg.SetFillColor(0)
    leg.SetLineColor(0)    
if __name__ == '__main__':
    main()   