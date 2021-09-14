"""

Author: Tao Yang <yangtao@ihep.ac.cn> and Yuhang Tan <tanyuhang@ihep.ac.cn>
Created [2021-07-07 Wed 08:47] 
Based on Raser C++ https://github.com/dt-np/raser


This program solves 2D Poisson's equation for Si & SiC PIN & LGAD doping profile

    - div grad u(x,y) = f(x,y)

on the unit interval with source f given by

    f(x,y) = e0*Q(x,y)/perm0/perm

    PIN:
    f(x,y) = e0*1e13*1e6/perm0/perm_sic

    LGAD:
    f(x,y) = (6.1942*pow(10,16)/sqrt(2*3.1415926)/0.13821*exp(-pow(x-0.67,2))/2/0.13821/0.13821 + pow(10,13))*(-(1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))

and boundary conditions given by
    u(x,y) = Bias_voltage for x = 0
    u(x,y) = 0 for x = thin

######################################################################
Warning:    N-type bulk is positive Spacecharge, negative Bias_voltage
            P-type bulk is negative Spacecharge, positive Bias_voltage
######################################################################

"""


import fenics
import numpy as np
import math
import random
import time
import sys
import ROOT
from array import array
import matplotlib.pyplot as plt
from raser import twoD_time 

""" Global Constant """
e0 = 1.60217733e-19 # C
perm0 = 8.854187817e-12 # F/m



""" Define Material """

class Material:

    def __init__(self,mat_name):
        self.mat_name = mat_name

    def mat_database(self):

        # material par
        self.si_par_dict = {'Permittivity' : 11.5,\
                             'Avalanche': 'vanOverstraeten',\
                             'Mobility' : 'unknown'\
                            }

        self.sic_par_dict = {'Permittivity' : 9.76,\
                             'Avalanche': 'unknown',\
                             'Mobility' : 'unknown'\
                            }

        # global data base
        self.mat_db_dict = {'SiC' : self.sic_par_dict,\
                            'Si' : self.si_par_dict\
                            }

        return self.mat_db_dict[self.mat_name]


""" Define Mobility """

class Mobility:
    def __init__(self,mat_name):
        self.mat_name = mat_name

    def cal_mobility(self, det, position, charge, electric_field):

        x = position[0]
        y = position[1]
        T = det.temperature # K
        E = electric_field  # V/cm

        doping_expr = det.doping_epr
        doping_expr = doping_expr.replace("x[1]","y")
        doping_expr = doping_expr.replace("sqrt","math.sqrt")
        doping_expr = doping_expr.replace("exp","math.exp")
        #print(doping_expr)
        Neff = abs(eval(doping_expr))

        # SiC mobility
        if(self.mat_name == 'SiC'):
            if(charge>0):
                alpha = 0.34
                ulp = 124 * math.pow(T / 300, -2)
                uminp = 15.9
                Crefp = 1.76e19
                betap = 1.213 * math.pow(T / 300.0, 0.17)
                vsatp = 2e7 * math.pow(T / 300.0, 0.52)
                lfm = uminp + ulp/(1.0 + math.pow(Neff*1e12 / Crefp, alpha))
                hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0 / betap))  

            if(charge<0):
                alpha = 0.61
                ulp = 947 * math.pow(T / 300, -2)
                Crefp = 1.94e19
                betap = 1 * math.pow(T / 300, 0.66)
                vsatp = 2e7 * math.pow(T / 300, 0.87)
                lfm = ulp/ (1 + math.pow(Neff*1e12 / Crefp, alpha))
                hfm = lfm / (math.pow(1.0 + math.pow(lfm * E / vsatp, betap), 1.0/betap))

        # Si mobility
        if(self.mat_name == 'Si'):
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

        return hfm



""" Define Avalanche """

class Avalanche:
        
    def __init__(self,model_name):
        self.model_name = model_name

    def cal_coefficient(self, electric_field, charges, temperature):

        coefficient = 0.

        E = electric_field # V/cm
        T = temperature # K

        # van Overstraeten – de Man Model
        if(self.model_name == 'vanOverstraeten'):

            hbarOmega = 0.063 # eV
            E0 = 4.0e5 # V/cm
            T0 = 293.0 # K
            k_T0 = 0.0257 # eV

            # electron
            if( charges < 0 ): 

                a_low = 7.03e5 # cm-1
                a_high = 7.03e5 # cm-1

                b_low = 1.232e6 # cm-1
                b_high = 1.232e6 # cm-1

                #
                # For BandgapDependence parameters
                #

                # Glambda = 62e-8 #cm
                # beta_low = 0.678925 # 1
                # beta_high = 0.678925 # 1

            # hole
            if( charges > 0 ): 

                a_low = 1.582e6 # cm-1
                a_high = 6.71e5 # cm-1

                b_low = 2.036e6 # cm-1
                b_high = 1.693e6 # cm-1

                Glambda = 45e-8 #cm

                beta_low = 0.815009 # 1
                beta_high =  0.677706 # 1

            Ggamma = math.tanh(hbarOmega/(2*k_T0))/math.tanh(hbarOmega/(2*k_T0*T/T0))
            
            if(E>1.75e05):
                if(E>E0):
                    coefficient = Ggamma*a_high*math.exp(-(Ggamma*b_high)/E)
                else:
                    coefficient = Ggamma*a_low*math.exp(-(Ggamma*b_low)/E)
            else:
                coefficient = 0.

        return coefficient




""" Define Detector Geometry """

class R2dDetector:

    def __init__(self,det_width,det_thin):

        self.det_width = det_width #um X-axis
        self.det_thin = det_thin #um Y-axis

        self.n_bin = 1000
        self.t_end = 3e-9
        self.positive_cu = ROOT.TH1F("charge+","Positive Current",self.n_bin,0,self.t_end)
        self.negtive_cu = ROOT.TH1F("charge-","Negative Current",self.n_bin,0,self.t_end)

        self.gain_positive_cu = ROOT.TH1F("gain_charge+","Current Contribution",self.n_bin,0,self.t_end)
        self.gain_negtive_cu = ROOT.TH1F("gain_charge-","Gain Negative Current",self.n_bin,0,self.t_end)

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



""" Define Fenics Possion Solver """

class FenicsPossion:

    def __init__(self,det):
        
        self.det = det

        # poential & field
        self.potential_value_2d = []
    
        self.electric_field_x_value = [ [] for n in range(self.det.ny+1) ]
        self.electric_field_y_value = [ [] for n in range(self.det.ny) ]

        self.electric_field_x_position = [ [] for n in range(self.det.ny+1) ]
        self.electric_field_y_position = [ [] for n in range(self.det.ny) ]

        # weighting poential & weighting field
        self.weighting_potential_value_2d = []
    
        self.weighting_electric_field_x_value = [ [] for n in range(self.det.ny+1) ]
        self.weighting_electric_field_y_value = [ [] for n in range(self.det.ny) ]

        self.weighting_electric_field_x_position = [ [] for n in range(self.det.ny+1) ]
        self.weighting_electric_field_y_position = [ [] for n in range(self.det.ny) ]

    def cal_possion(self):
        
        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        # Create mesh and function space
        mesh = fenics.RectangleMesh(fenics.Point(0, 0), fenics.Point(width, thin), nx, ny)
        V = fenics.FunctionSpace(mesh, "P", 1)

        # Define boundary condition
        u_D = fenics.Expression('x[1] < tol? det_voltage : 0', degree = 2,tol = 1E-14,det_voltage = self.det.bias_voltage)

        def boundary(x, on_boundary):
            return abs(x[1])<1E-14 or abs(x[1]-thin)<1E-14

        bc = fenics.DirichletBC(V, u_D, boundary)

        # Define variational problem
        u = fenics.TrialFunction(V)
        v = fenics.TestFunction(V)

        f = fenics.Expression(self.det.doping_epr,degree=2)
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx #+ g*v*ds

        # Compute solution
        u = fenics.Function(V)
        fenics.solve(a == L, u, bc)

        potential_value_1d = u.compute_vertex_values()
        potential_value_2d = np.array(potential_value_1d).reshape(ny+1,nx+1)

        self.potential_value_2d = potential_value_2d

        # print(self.potential_value_2d)


    def cal_weighting_possion(self):

        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        # Create mesh and function space
        mesh = fenics.RectangleMesh(fenics.Point(0, 0), fenics.Point(width, thin), nx, ny)
        V = fenics.FunctionSpace(mesh, "P", 1)

        # Define boundary condition
        u_D = fenics.Expression('x[1] < tol? det_voltage : 0', degree = 2,tol = 1E-14,det_voltage = self.det.bias_voltage/abs(self.det.bias_voltage))

        def boundary(x, on_boundary):
            return abs(x[1])<1E-14 or abs(x[1]-thin)<1E-14

        bc = fenics.DirichletBC(V, u_D, boundary)

        # Define variational problem
        u = fenics.TrialFunction(V)
        v = fenics.TestFunction(V)

        f = fenics.Constant(0)
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx #+ g*v*ds

        # Compute solution
        u = fenics.Function(V)
        fenics.solve(a == L, u, bc)

        weighting_potential_value_1d = u.compute_vertex_values()
        weighting_potential_value_2d = np.array(weighting_potential_value_1d).reshape(ny+1,nx+1)

        self.weighting_potential_value_2d = weighting_potential_value_2d

        return weighting_potential_value_2d

    def cal_electric_field(self):

        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        x_step = width/nx
        y_step = thin/ny

        # x direction
        for j in range(ny+1):
            for i in range(nx):

                # electric field
                tmp_xpos = 0.5*x_step*(2*i+1)
                tmp_xef = (self.potential_value_2d[j][i] - self.potential_value_2d[j][i+1])/x_step
                self.electric_field_x_position[j].append(tmp_xpos)
                self.electric_field_x_value[j].append(tmp_xef)

                # weighting field
                tmp_wxpos = 0.5*x_step*(2*i+1)
                tmp_xwef = (self.weighting_potential_value_2d[j][i] - self.weighting_potential_value_2d[j][i+1])/x_step
                self.weighting_electric_field_x_position[j].append(tmp_wxpos)
                self.weighting_electric_field_x_value[j].append(tmp_xwef)


        # y direction 
        for j in range(ny):
            for i in range(nx+1):

                # electric field
                tmp_ypos = 0.5*y_step*(2*j+1)
                tmp_yef = (self.potential_value_2d[j][i]- self.potential_value_2d[j+1][i])/y_step
                self.electric_field_y_position[j].append(tmp_ypos)
                self.electric_field_y_value[j].append(tmp_yef)

                # weighting field
                tmp_wypos = 0.5*y_step*(2*j+1)
                tmp_ywef = (self.weighting_potential_value_2d[j][i] - self.weighting_potential_value_2d[j+1][i])/y_step
                self.weighting_electric_field_y_position[j].append(tmp_wypos)
                self.weighting_electric_field_y_value[j].append(tmp_ywef)     



    def cal_point_field(self,px_point,py_point,input_value):

        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        x_step = width/nx
        y_step = thin/ny

        #Interpolation method 
        # if px_point == px_point:
        #     px_point = px_point
        # else:
        #     px_point = 0
        rex_value=px_point%x_step
        nx_value=int(px_point/x_step)
        rey_value=py_point%y_step
        ny_value=int(py_point/y_step)

        if(rex_value>x_step/2):
            e_v_x1=rex_value-x_step/2
            nx1_v=nx_value
            nx2_v=nx_value+1
        else:
            e_v_x1=rex_value+x_step/2
            e_v_x2=x_step-e_v_x1
            nx1_v=nx_value-1
            nx2_v=nx_value

        if(rey_value>y_step/2):
            e_v_y1=rey_value-y_step/2
            ny1_v=ny_value
            ny2_v=ny_value+1
        else:
            e_v_y1=rey_value+y_step/2
            e_v_y2=y_step-e_v_y1
            ny1_v=ny_value-1
            ny2_v=ny_value

        if (nx_value<=0):
            r_u=0
            nx1_v=nx2_v
        elif (nx_value>=nx-1):
            r_u=0
            nx2_v=nx1_v
        else:
            r_u=e_v_x1/x_step

        if (ny_value<=0):
            r_t=0
            ny1_v=ny2_v
        elif (ny_value>=ny-1):
            r_t=0
            ny2_v=ny1_v
        else:
            r_t=e_v_y1/y_step

        value_11=input_value[ny1_v][nx1_v]
        value_21=input_value[ny1_v][nx2_v]
        value_12=input_value[ny1_v][nx1_v]
        value_22=input_value[ny1_v][nx2_v]
        out_field=0.0
        out_field=(1-r_u)*(1-r_t)*value_11
        out_field+=r_u*(1-r_t)*value_21
        out_field+=r_u*r_t*value_22
        out_field+=(1-r_u)*r_t*value_12

        return out_field        


    def solve(self):

        self.cal_possion()
        self.cal_weighting_possion()
        self.cal_electric_field()

    def draw(self):

        cutline = int(self.det.nx/2.0)

        # plot electric field at x = middle
        ep = array( 'd' )
        ev = array( 'd' )

        for i in range(self.det.ny):

            ep.append(self.electric_field_y_position[i][cutline])
            efx = self.cal_point_field(self.det.det_width/2.0, self.electric_field_y_position[i][cutline], self.electric_field_x_value)
            efy = self.cal_point_field(self.det.det_width/2.0, self.electric_field_y_position[i][cutline], self.electric_field_y_value)
            ef = np.array([efx,efy])

            ef_value = np.linalg.norm(ef)#*1e4 #V/cm

            ev.append(ef_value)

        # print(ep)
        # print(ev)

        g_e = ROOT.TGraph(self.det.ny,ep,ev)

        g_e.SetLineColor(600)
        g_e.SetLineWidth(4)
        g_e.SetTitle( 'Electric Field at Cut Line' )
        g_e.GetXaxis().SetTitle( 'Dpeth [um]' )
        g_e.GetXaxis().SetRangeUser(0,self.det.det_thin)
        g_e.GetYaxis().SetTitle( 'E [V/um]' )

        g_e.GetXaxis().CenterTitle()
        g_e.GetXaxis().SetTitleOffset(1.8)
        g_e.GetXaxis().SetTitleSize(0.05)
        g_e.GetXaxis().SetLabelSize(0.05)
        #g_e.GetXaxis().SetNdivisions(505)

        g_e.GetYaxis().CenterTitle()
        g_e.GetYaxis().SetTitleOffset(1.8)
        g_e.GetYaxis().SetTitleSize(0.05)
        g_e.GetYaxis().SetLabelSize(0.05)
        #g_e.GetYaxis().SetNdivisions(505)

        # plot weighting electric field at x = middle
        # wep = [depth[250] for depth in self.weighting_electric_field_y_position]
        # wev = [wef[250] for wef in self.weighting_electric_field_y_value]
        # g_we = ROOT.TGraph(len(wep),wep,wev)
        
        c = ROOT.TCanvas( 'c', 'c',500, 500 )
        c.SetGrid()
        c.SetLeftMargin(0.18)
        c.SetBottomMargin(0.18)
        # c.Divide(1,2)

        c.cd()
        g_e.Draw()
        c.Modified()
        c.Update()
        c.SaveAs("electricfield.pdf")


        # plt.figure(figsize=(8,4), dpi=80)
        # plt.figure(1)
# # 
        # ax1 = plt.subplot(121)        
        # plt.plot( [depth[250] for depth in self.electric_field_y_position], [ef[250] for ef in self.electric_field_y_value])
        # plt.xlabel('depth [um]')
        # plt.ylabel('electric field [V/um]')
# # 
        # ax2 = plt.subplot(122)
        # plt.plot( [wdepth[250] for wdepth in self.weighting_electric_field_y_position], [wef[250] for wef in self.weighting_electric_field_y_value])
        # plt.xlabel('depth [um]')
        # plt.ylabel('weighting electric field [V/um]')
        # plt.tight_layout()
        # plt.show()
        



""" Define Track """

class Tracks:

    # MIPs particle
    def mips(self,track_entry,track_exit,n_div):

        # initial
        self.track_entry = track_entry
        self.track_exit = track_exit
        self.track_position = [ [] for n in range(n_div-1) ]
        self.delta_track = []

        # start position
        track_position_x=track_entry[0]
        track_position_y=track_entry[1]

        for i in range(n_div-1):

            x_div_point = (track_exit[0]-track_entry[0])/n_div*i+track_entry[0]+(track_exit[0]-track_entry[0])/(2*n_div)
            y_div_point = (track_exit[1]-track_entry[1])/n_div*i+track_entry[1]+(track_exit[1]-track_entry[1])/(2*n_div)

            self.track_position[i].append(x_div_point)
            self.track_position[i].append(y_div_point)

            self.delta_track.append(math.sqrt(math.pow(x_div_point-track_position_x,2)+math.pow(y_div_point-track_position_y,2)))

            track_position_x = x_div_point
            track_position_y = y_div_point
    
    def mips_ionized(self):

        self.ionized_pairs = []
        self.ionized_total_pairs = 0.

        #Drift distance: self.delta_track
        loss_energy = 8.4 #silicon carbide /um
        Myf = ROOT.TFile("SiC_5um.root")
        hist = Myf.Get("Edep_device")

        for i in range(len(self.track_position)):
            
            gRandom = ROOT.TRandom3(0)
            ran_energy = hist.GetRandom(gRandom)
            n_pairs = self.delta_track[i]*ran_energy*1e6/(5*loss_energy)
            self.ionized_pairs.append(n_pairs)
            self.ionized_total_pairs += n_pairs


    def mips_laudau(self):

        loss_energy = 8.4
        d_x = math.pow(self.track_exit[0]-self.track_entry[0],2)
        d_y = math.pow(self.track_exit[1]-self.track_entry[1],2)
        d_dis = math.sqrt(d_x+d_y)
        d_mpv = 0.04 * math.log(d_dis) + 0.27
        d_FWHM = 0.31 * math.pow(d_dis, -0.17)
        d_da = d_FWHM / 4.
        gRandom = ROOT.TRandom3(0)
        LanConst = gRandom.Landau(d_mpv, d_da)
        if (LanConst > 5.0 * d_mpv):
            LanConst = gRandom.Landau(d_mpv, d_da)
        Lan_pairs = LanConst*1000/loss_energy*d_dis

        return Lan_pairs


    # Laser-TCT...
    def tct_mips(self,i_z, i_r, zlen, rlen, nocarrier):
        # initial
        self.track_position = [ [ [] for n in range(i_z) ] for m in range (i_r)]

        for i in range(i_r):
            self.track_entry = [0, i]
            self.track_exit  = [zlen, i]
            track_position_x = self.track_entry[0]
            track_position_y = self.track_entry[1]
            for j in range(i_z):
                x_div_point = j*(zlen/i_z)+ (zlen/(i_z*2))
                y_div_point = i*(rlen/i_r)+ (rlen/(i_r*2))
                self.track_position[i][j].append(x_div_point)
                self.track_position[i][j].append(y_div_point)
    
    # TPA-TCT...




'''TCT Drifts'''

class TCTDrifts:

    def __init__(self,track):

        # initial tracks
        self.delta_track_info_dic_n = {}
        self.delta_track_info_dic_p = {}

        # gain tracks
        # self.delta_gain_track_info_dic_n_n = {}
        # self.delta_gain_track_info_dic_n_p = {}
        # self.delta_gain_track_info_dic_p_n = {}
        # self.delta_gain_track_info_dic_p_p = {}
        self.delta_gain_track_info_dic = {}

        self.track_time = 0.
        self.track_x = 0.
        self.track_y = 0.
        self.track_charges = 0.
        self.track_current = 0.

        self.end_condition = 0

        self.delta_x=0.
        self.delta_y=0.
        self.diff_x=0.
        self.diff_y=0.

        self.s_time = 0.

        self.s_gain = 1.

        for n in range(len(track.track_position)*len(track.track_position[0])):

            # initial tracks
            self.delta_track_info_dic_n["tk_"+str(n+1)] = [ [] for n in range(6) ]   # track_time, track_x, track_y, track_charges, track_current, track_gain
            self.delta_track_info_dic_p["tk_"+str(n+1)] = [ [] for n in range(6) ] 

            # gain tracks
            # self.delta_gain_track_info_dic_n_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            # self.delta_gain_track_info_dic_n_p["tk_"+str(n+1)] = [ [] for n in range(5) ] 
            # self.delta_gain_track_info_dic_p_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            # self.delta_gain_track_info_dic_p_p["tk_"+str(n+1)] = [ [] for n in range(5) ] 

    
    def drift_diffusion(self,det,fen):

        my_mobility = Mobility(det.mat_name)

        self.muhh=1650   # mobility related with the magnetic field (now silicon useless)
        self.muhe=310
        self.BB=np.array([0,0])
        self.sstep=1 # drift step
        self.kboltz=8.617385e-5 # eV/K
        self.max_drift_len=1e9 # maximum driftlength [um]

        self.delta_x=0.
        self.delta_y=0. 
        

        #
        # drift
        #

        # magnetic field effect

        # for hole
        if(self.charges)>0:
            FF=self.e_field+self.muhh*np.cross(self.e_field,self.BB)
        # for electron    
        else:
            FF=self.e_field-self.muhe*np.cross(self.e_field,self.BB)

        # the delta x with electric field and magnetic field for unit e or h

        if(np.linalg.norm(FF)!=0):
            self.delta_x = self.sstep*(self.charges/abs(self.charges))*FF[0]/np.linalg.norm(FF)
            self.delta_y = self.sstep*(self.charges/abs(self.charges))*FF[1]/np.linalg.norm(FF)
        else:
            self.delta_x=0
            self.delta_y=0 

        # cal current step ef & wef
        if self.delta_x < -1e-07:
            self.delta_x = 0
        efx = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.electric_field_x_value)
        efy = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.electric_field_y_value)
        ef = np.array([efx,efy])

        ef_value = np.linalg.norm(ef)*1e4 #V/cm

        wefx = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.weighting_electric_field_x_value)
        wefy = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.weighting_electric_field_y_value)
        wef = np.array([wefx,wefy])

        self.wef_value = np.linalg.norm(wef)


        pos = [self.track_x+self.delta_x, self.track_y+self.delta_y]

        self.drift_velocity = my_mobility.cal_mobility(det, pos, self.charges, ef_value)*ef_value
        # print(self.drift_velocity)

        self.e_field = ef
        self.we_field = wef

        #
        # diffusion
        #
        if(self.drift_velocity == 0):
            self.delta_x=0
            self.delta_y=0
            self.diff_x=0
            self.diff_y=0
            self.end_condition = 9
        else:
            DiffOffField=8*1e4 #V/cm
            if(ef_value<DiffOffField):
                self.s_time=self.sstep*1e-4/self.drift_velocity
                s_sigma= math.sqrt(2*self.kboltz*my_mobility.cal_mobility(det, pos, self.charges, ef_value)*det.temperature*self.s_time)
                # self.dif_x=random.gauss(0,s_sigma)*1e4
                # self.dif_y=random.gauss(0,s_sigma)*1e4
                self.dif_x=0.
                self.dif_y=0.

            else:
                self.dif_x=0.0
                self.dif_y=0.0 


        #
        # multiplication
        #
        my_avalanche = Avalanche('vanOverstraeten')
        tmp_coefficient = my_avalanche.cal_coefficient(ef_value,self.charges,det.temperature)

        self.s_gain = math.exp(self.sstep*1e-4*tmp_coefficient)

        


    def update_track_info(self):

        if(self.charges>0):
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][0].append(self.track_time)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][1].append(self.track_x)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][2].append(self.track_y)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][3].append(self.track_charges)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][4].append(self.track_current)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][5].append(self.track_gain)
 
        if(self.charges<0):
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][0].append(self.track_time)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][1].append(self.track_x)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][2].append(self.track_y)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][3].append(self.track_charges)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][4].append(self.track_current)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][5].append(self.track_gain)

    def update_gain_track_info(self):
        
        self.delta_gain_track_info_dic[self.track_name][0].append(self.track_time)
        self.delta_gain_track_info_dic[self.track_name][1].append(self.track_x)
        self.delta_gain_track_info_dic[self.track_name][2].append(self.track_y)
        self.delta_gain_track_info_dic[self.track_name][3].append(self.track_charges)
        self.delta_gain_track_info_dic[self.track_name][4].append(self.track_current)


    def update_step(self,det):

        self.track_time = self.track_time + self.s_time
        print(self.delta_x)
        self.track_x = self.track_x + self.delta_x + self.diff_x
        self.track_y = self.track_y + self.delta_y + self.diff_y
        self.track_charges = self.charges
        self.track_gain *= self.s_gain


        if(self.track_x>=det.det_width):
            self.track_x = det.det_width
        
        if(self.track_x<0):
            self.track_x = 0

        if(self.track_y>=det.det_thin):
            self.track_y = det.det_thin
        
        if(self.track_y<0):
            self.track_y = 0
        
    def update_end_condition(self):

        # if(self.we_field>(1-1e-5)):
        #     self.end_condition=1
        if(self.track_x<=0):
            self.end_condition=2
        if(self.track_y<=0):
            self.end_condition=3
        # else: 
        #     self.end_condition=0
            
        # if(self.path_len>self.max_drift_len):
        #     self.end_condition=6
        # if(self.n_step>10000):
        #     self.end_condition=7
        

    def cal_current(self, track, fen, det, nocarrier):


        #
        # initial carrier track
        #

        for i in range(len(track.track_position)):

            for m in range(len(track.track_position[i])):
                
                self.track_number = i+1

                for j in range(2):
                    
                    if(j==0):
                        
                        self.charges=1*nocarrier[m,i] # hole

                    if(j==1):
                        self.charges=-1*nocarrier[m,i] # electron
                
                    self.track_time = 0.
                    self.track_x = track.track_position[i][m][0]
                    self.track_y = track.track_position[i][m][1]
                    self.track_charges = 0.
                    self.track_current = 0.
                    self.track_gain = 1.

                    #print("current position "+str(i)+" : "+str(self.track_x)+","+str(self.track_y))
                
                
                    # print(self.end_condition)
                    self.end_condition = 0
                    while(self.end_condition==0):
                        if(self.track_y>=(det.det_thin-1) or self.track_x>=(det.det_width-1)):
                            self.end_condition=4

                        else:
                            efx = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_x_value)
                            efy = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_y_value)

                            ef = np.array([efx,efy])
                            ef_value = np.linalg.norm(ef)*1e4

                            self.e_field = np.array([efx,efy])

                            wefx = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_x_value)
                            wefy = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_y_value)

                            wef = np.array([wefx,wefy])
                            wef_value = np.linalg.norm(wef)*1e4

                            self.we_field = np.array([wefx,wefy])                                
                            print([i,m])
                            self.drift_diffusion(det,fen)

                            # SR current

                            self.track_current = abs(self.charges*e0*self.drift_velocity*wef_value)

                            self.update_track_info()
                            self.update_step(det)
                            self.update_end_condition()



        #
        # gian carrier track
        #

        # get gain tracks start info
        self.gain_track_info_list = [] #[[name,time,x,y,charges]]

        for i in range(len(track.track_position)):

            for m in range(len(track.track_position[i])):

                for j in range(len(self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][0])):
                    
                    if(self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][5][j]>1.0):
                        tmp_gain_time = self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][0][j]
                        tmp_gain_x = self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][1][j]
                        tmp_gain_y = self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][2][j]
                        tmp_gain_pairs = abs(self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][3][j]*(np.max(self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][5])-1))
                        tmp_gain_current = 0.

                        self.gain_track_info_list.append(["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_p_n",tmp_gain_time,tmp_gain_x,tmp_gain_y,-tmp_gain_pairs])

                        self.gain_track_info_list.append(["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_p_p",tmp_gain_time,tmp_gain_x,tmp_gain_y,tmp_gain_pairs])

                        self.delta_gain_track_info_dic["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_p_n"] = [ [] for n in range(5) ]
                        self.delta_gain_track_info_dic["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_p_p"] = [ [] for n in range(5) ]

                        break

                for k in range(len(self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][0])):

                    if(self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][5][k]>1.0):

                        tmp_gain_time = self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][0][k]
                        tmp_gain_x = self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][1][k]
                        tmp_gain_y = self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][2][k]
                        tmp_gain_pairs = abs(self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][3][k]*(np.max(self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][5])-1))
                        tmp_gain_current = 0.

                        self.gain_track_info_list.append(["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_n_n",tmp_gain_time,tmp_gain_x,tmp_gain_y,-tmp_gain_pairs])

                        self.gain_track_info_list.append(["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_n_p",tmp_gain_time,tmp_gain_x,tmp_gain_y,tmp_gain_pairs])

                        self.delta_gain_track_info_dic["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_n_n"] = [ [] for n in range(5) ]
                        self.delta_gain_track_info_dic["tk_"+str(i*(len(track.track_position[0]))+m+1)+"_n_p"] = [ [] for n in range(5) ]

                        break
        print("ttest3")
        #
        # cal gain current
        #

        for i in range(len(self.gain_track_info_list)):

            self.charges = self.gain_track_info_list[i][4]

            self.track_name = self.gain_track_info_list[i][0]
            self.track_time = self.gain_track_info_list[i][1]
            self.track_x = self.gain_track_info_list[i][2]
            self.track_y = self.gain_track_info_list[i][3]
            self.track_charges = self.gain_track_info_list[i][4]
            self.track_current = 0.

            self.end_condition = 0
            while(self.end_condition==0):

                if(self.track_y>=(det.det_thin-1) or self.track_x>=(det.det_width-1)):
                     
                    self.end_condition=4

                else:
                    efx = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_x_value)
                    efy = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_y_value)
                    ef = np.array([efx,efy])
                    ef_value = np.linalg.norm(ef)*1e4
                    self.e_field = np.array([efx,efy])
                    wefx = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_x_value)
                    wefy = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_y_value)
                    wef = np.array([wefx,wefy])
                    wef_value = np.linalg.norm(wef)*1e4
                    self.we_field = np.array([wefx,wefy])
                                                   
                    self.drift_diffusion(det,fen)
                    # SR current
                    self.track_current = abs(self.charges*e0*self.drift_velocity*wef_value)

                    self.update_gain_track_info()
                    self.update_step(det)
                    self.update_end_condition()      
        print("ttest4")      


        det.positive_cu.Reset()
        det.negtive_cu.Reset()

        det.gain_positive_cu.Reset()
        det.gain_negtive_cu.Reset()

        det.gain_n_n_cu.Reset()
        det.gain_n_p_cu.Reset()
        det.gain_p_n_cu.Reset()
        det.gain_p_p_cu.Reset()

        det.sum_cu.Reset()

        temp_positive_cu = ROOT.TH1F("temp+","temp+",det.n_bin,0,det.t_end)
        temp_negitive_cu = ROOT.TH1F("temp-","temp-",det.n_bin,0,det.t_end)
        temp_sum_cu = ROOT.TH1F("temp_sum","temp_sum",det.n_bin,0,det.t_end)
        
        #
        # initial current
        #
        for i in range(len(track.track_position)):

            for m in range(len(track.track_position[i])):
                
                for j in range(len(self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][0])):
                    temp_positive_cu.Fill(self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][0][j], self.delta_track_info_dic_p["tk_"+str(i*(len(track.track_position[0]))+m+1)][4][j])

                for k in range(len(self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][0])):
                    temp_negitive_cu.Fill(self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][0][k], self.delta_track_info_dic_n["tk_"+str(i*(len(track.track_position[0]))+m+1)][4][k])

                det.positive_cu.Add(temp_positive_cu)
                det.negtive_cu.Add(temp_negitive_cu)

                temp_positive_cu.Reset()
                temp_negitive_cu.Reset()
        print("ttest5")


        #
        # gain current
        #

        temp_gain_cu = ROOT.TH1F("temp_gain","temp_gain",det.n_bin,0,det.t_end)

        for i in range(len(self.gain_track_info_list)):

            tmp_track_name = self.gain_track_info_list[i][0]

            for j in range(len(self.delta_gain_track_info_dic[tmp_track_name][0])):
                
                temp_gain_cu.Fill(self.delta_gain_track_info_dic[tmp_track_name][0][j],self.delta_gain_track_info_dic[tmp_track_name][4][j])

            if(tmp_track_name[-1]=='p'):
                det.gain_positive_cu.Add(temp_gain_cu)

            if(tmp_track_name[-1]=='n'):
                det.gain_negtive_cu.Add(temp_gain_cu)


            if(tmp_track_name[-3:]=='n_n'):
                det.gain_n_n_cu.Add(temp_gain_cu)

            if(tmp_track_name[-3:]=='n_p'):
                det.gain_n_p_cu.Add(temp_gain_cu)

            if(tmp_track_name[-3:]=='p_n'):
                det.gain_p_n_cu.Add(temp_gain_cu)

            if(tmp_track_name[-3:]=='p_p'):
                det.gain_p_p_cu.Add(temp_gain_cu)


            temp_gain_cu.Reset()
        print("ttest6")


        det.sum_cu.Add(det.positive_cu)
        det.sum_cu.Add(det.negtive_cu)
        det.sum_cu.Add(det.gain_positive_cu)
        det.sum_cu.Add(det.gain_negtive_cu)
        
        
    def draw_drift_path(self,det):
        # ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "canvas1", 200,10,1000, 1000)
        mg = ROOT.TMultiGraph("mg","")
        x_array=array('f')
        y_array=array('f')
        for i in range(len(self.delta_track_info_dic_p)):
            n=len(self.delta_track_info_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.delta_track_info_dic_p["tk_"+str(i+1)][1])
                y_array.extend(self.delta_track_info_dic_p["tk_"+str(i+1)][2])             
                gr_p = ROOT.TGraph(n,x_array,y_array)
                gr_p.SetMarkerColor(4)
                gr_p.SetLineColor(4)
                gr_p.SetLineStyle(1)
                mg.Add(gr_p)
                del x_array[:]
                del y_array[:]
        for j in range(len(self.delta_track_info_dic_n)):
            m=len(self.delta_track_info_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.delta_track_info_dic_n["tk_"+str(j+1)][1])
                y_array.extend(self.delta_track_info_dic_n["tk_"+str(j+1)][2])             
                gr_n = ROOT.TGraph(m,x_array,y_array)
                gr_n.SetMarkerColor(2)
                gr_n.SetLineColor(2)
                gr_n.SetLineStyle(1)
                mg.Add(gr_n)
                del x_array[:]
                del y_array[:]
        mg.GetXaxis().SetRangeUser(0,det.det_width)
        mg.GetYaxis().SetRangeUser(0,det.det_thin)
        mg.Draw("APL")
        c1.SaveAs("drift_path.pdf")

""" Define Carriers Drifts """
class Drifts:

    def __init__(self,track):

        # initial tracks
        self.delta_track_info_dic_n = {}
        self.delta_track_info_dic_p = {}

        # gain tracks
        # self.delta_gain_track_info_dic_n_n = {}
        # self.delta_gain_track_info_dic_n_p = {}
        # self.delta_gain_track_info_dic_p_n = {}
        # self.delta_gain_track_info_dic_p_p = {}
        self.delta_gain_track_info_dic = {}

        self.track_time = 0.
        self.track_x = 0.
        self.track_y = 0.
        self.track_charges = 0.
        self.track_current = 0.

        self.end_condition = 0

        self.delta_x=0.
        self.delta_y=0.
        self.diff_x=0.
        self.diff_y=0.

        self.s_time = 0.

        self.s_gain = 1.

        for n in range(len(track.track_position)):
            # initial tracks
            self.delta_track_info_dic_n["tk_"+str(n+1)] = [ [] for n in range(6) ]   # track_time, track_x, track_y, track_charges, track_current, track_gain
            self.delta_track_info_dic_p["tk_"+str(n+1)] = [ [] for n in range(6) ] 

            # gain tracks
            # self.delta_gain_track_info_dic_n_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            # self.delta_gain_track_info_dic_n_p["tk_"+str(n+1)] = [ [] for n in range(5) ] 
            # self.delta_gain_track_info_dic_p_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            # self.delta_gain_track_info_dic_p_p["tk_"+str(n+1)] = [ [] for n in range(5) ] 

    
    def drift_diffusion(self,det,fen):

        my_mobility = Mobility(det.mat_name)

        self.muhh=1650   # mobility related with the magnetic field (now silicon useless)
        self.muhe=310
        self.BB=np.array([0,0])
        self.sstep=0.1 # drift step
        self.kboltz=8.617385e-5 # eV/K
        self.max_drift_len=1e9 # maximum driftlength [um]

        self.delta_x=0.
        self.delta_y=0. 
        

        #
        # drift
        #

        # magnetic field effect

        # for hole
        if(self.charges)>0:
            FF=self.e_field+self.muhh*np.cross(self.e_field,self.BB)
        # for electron    
        else:
            FF=self.e_field-self.muhe*np.cross(self.e_field,self.BB)

        # the delta x with electric field and magnetic field for unit e or h

        if(np.linalg.norm(FF)!=0):
            self.delta_x = self.sstep*(self.charges/abs(self.charges))*FF[0]/np.linalg.norm(FF)
            self.delta_y = self.sstep*(self.charges/abs(self.charges))*FF[1]/np.linalg.norm(FF)
        else:
            self.delta_x=0
            self.delta_y=0 

        # cal current step ef & wef
        efx = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.electric_field_x_value)
        efy = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.electric_field_y_value)
        ef = np.array([efx,efy])

        ef_value = np.linalg.norm(ef)*1e4 #V/cm

        wefx = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.weighting_electric_field_x_value)
        wefy = fen.cal_point_field(self.track_x+self.delta_x, self.track_y+self.delta_y,fen.weighting_electric_field_y_value)
        wef = np.array([wefx,wefy])

        self.wef_value = np.linalg.norm(wef)


        pos = [self.track_x+self.delta_x, self.track_y+self.delta_y]

        self.drift_velocity = my_mobility.cal_mobility(det, pos, self.charges, ef_value)*ef_value
        # print(self.drift_velocity)

        self.e_field = ef
        self.we_field = wef

        #
        # diffusion
        #
        if(self.drift_velocity == 0):
            self.delta_x=0
            self.delta_y=0
            self.diff_x=0
            self.diff_y=0
            self.end_condition = 9
        else:
            DiffOffField=8*1e4 #V/cm
            if(ef_value<DiffOffField):
                self.s_time=self.sstep*1e-4/self.drift_velocity
                s_sigma= math.sqrt(2*self.kboltz*my_mobility.cal_mobility(det, pos, self.charges, ef_value)*det.temperature*self.s_time)
                # self.dif_x=random.gauss(0,s_sigma)*1e4
                # self.dif_y=random.gauss(0,s_sigma)*1e4
                self.dif_x=0.
                self.dif_y=0.

            else:
                self.dif_x=0.0
                self.dif_y=0.0 


        #
        # multiplication
        #
        my_avalanche = Avalanche('vanOverstraeten')
        tmp_coefficient = my_avalanche.cal_coefficient(ef_value,self.charges,det.temperature)

        self.s_gain = math.exp(self.sstep*1e-4*tmp_coefficient)

        


    def update_track_info(self):

        if(self.charges>0):
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][0].append(self.track_time)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][1].append(self.track_x)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][2].append(self.track_y)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][3].append(self.track_charges)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][4].append(self.track_current)
            self.delta_track_info_dic_p["tk_"+str(self.track_number)][5].append(self.track_gain)
 
        if(self.charges<0):
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][0].append(self.track_time)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][1].append(self.track_x)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][2].append(self.track_y)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][3].append(self.track_charges)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][4].append(self.track_current)
            self.delta_track_info_dic_n["tk_"+str(self.track_number)][5].append(self.track_gain)

    def update_gain_track_info(self):
        
        self.delta_gain_track_info_dic[self.track_name][0].append(self.track_time)
        self.delta_gain_track_info_dic[self.track_name][1].append(self.track_x)
        self.delta_gain_track_info_dic[self.track_name][2].append(self.track_y)
        self.delta_gain_track_info_dic[self.track_name][3].append(self.track_charges)
        self.delta_gain_track_info_dic[self.track_name][4].append(self.track_current)


    def update_step(self,det):

        self.track_time = self.track_time + self.s_time
        self.track_x = self.track_x + self.delta_x + self.diff_x
        self.track_y = self.track_y + self.delta_y + self.diff_y
        self.track_charges = self.charges
        self.track_gain *= self.s_gain


        if(self.track_x>=det.det_width):
            self.track_x = det.det_width
        
        if(self.track_x<0):
            self.track_x = 0

        if(self.track_y>=det.det_thin):
            self.track_y = det.det_thin
        
        if(self.track_y<0):
            self.track_y = 0
        
    def update_end_condition(self):

        # if(self.we_field>(1-1e-5)):
        #     self.end_condition=1
        if(self.track_x<=0):
            self.end_condition=2
        if(self.track_y<=0):
            self.end_condition=3
        # else: 
        #     self.end_condition=0
            
        # if(self.path_len>self.max_drift_len):
        #     self.end_condition=6
        # if(self.n_step>10000):
        #     self.end_condition=7
        

    def cal_current(self, track, fen, det):


        #
        # initial carrier track
        #
        track.mips_ionized()
        
        for i in range(len(track.track_position)):

            self.track_number = i+1

            for j in range(2):

                if(j==0):
                    self.charges=1*track.ionized_pairs[i] # hole

                if(j==1):
                    self.charges=-1*track.ionized_pairs[i] # electron
                
                self.track_time = 0.
                self.track_x = track.track_position[i][0]
                self.track_y = track.track_position[i][1]
                self.track_charges = 0.
                self.track_current = 0.
                self.track_gain = 1.

                #print("current position "+str(i)+" : "+str(self.track_x)+","+str(self.track_y))
                
                
                # print(self.end_condition)
                self.end_condition = 0
                while(self.end_condition==0):

                    if(self.track_y>=(det.det_thin-1) or self.track_x>=(det.det_width-1)):
                         
                        self.end_condition=4

                    else:
                        efx = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_x_value)
                        efy = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_y_value)

                        ef = np.array([efx,efy])
                        ef_value = np.linalg.norm(ef)*1e4

                        self.e_field = np.array([efx,efy])

                        wefx = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_x_value)
                        wefy = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_y_value)

                        wef = np.array([wefx,wefy])
                        wef_value = np.linalg.norm(wef)*1e4

                        self.we_field = np.array([wefx,wefy])                                

                        self.drift_diffusion(det,fen)

                        # SR current

                        self.track_current = abs(self.charges*e0*self.drift_velocity*wef_value)

                        self.update_track_info()
                        self.update_step(det)
                        self.update_end_condition()


        #
        # gian carrier track
        #

        # get gain tracks start info
        self.gain_track_info_list = [] #[[name,time,x,y,charges]]

        for i in range(len(track.track_position)):

            for j in range(len(self.delta_track_info_dic_p["tk_"+str(i+1)][0])):
                
                if(self.delta_track_info_dic_p["tk_"+str(i+1)][5][j]>1.0):

                    tmp_gain_time = self.delta_track_info_dic_p["tk_"+str(i+1)][0][j]
                    tmp_gain_x = self.delta_track_info_dic_p["tk_"+str(i+1)][1][j]
                    tmp_gain_y = self.delta_track_info_dic_p["tk_"+str(i+1)][2][j]
                    tmp_gain_pairs = abs(self.delta_track_info_dic_p["tk_"+str(i+1)][3][j]*(np.max(self.delta_track_info_dic_p["tk_"+str(i+1)][5])-1))
                    tmp_gain_current = 0.

                    self.gain_track_info_list.append(["tk_"+str(i+1)+"_p_n",tmp_gain_time,tmp_gain_x,tmp_gain_y,-tmp_gain_pairs])

                    self.gain_track_info_list.append(["tk_"+str(i+1)+"_p_p",tmp_gain_time,tmp_gain_x,tmp_gain_y,tmp_gain_pairs])

                    self.delta_gain_track_info_dic["tk_"+str(i+1)+"_p_n"] = [ [] for n in range(5) ]
                    self.delta_gain_track_info_dic["tk_"+str(i+1)+"_p_p"] = [ [] for n in range(5) ]

                    break

            for k in range(len(self.delta_track_info_dic_n["tk_"+str(i+1)][0])):

                if(self.delta_track_info_dic_n["tk_"+str(i+1)][5][k]>1.0):

                    tmp_gain_time = self.delta_track_info_dic_n["tk_"+str(i+1)][0][k]
                    tmp_gain_x = self.delta_track_info_dic_n["tk_"+str(i+1)][1][k]
                    tmp_gain_y = self.delta_track_info_dic_n["tk_"+str(i+1)][2][k]
                    tmp_gain_pairs = abs(self.delta_track_info_dic_n["tk_"+str(i+1)][3][k]*(np.max(self.delta_track_info_dic_n["tk_"+str(i+1)][5])-1))
                    tmp_gain_current = 0.

                    self.gain_track_info_list.append(["tk_"+str(i+1)+"_n_n",tmp_gain_time,tmp_gain_x,tmp_gain_y,-tmp_gain_pairs])

                    self.gain_track_info_list.append(["tk_"+str(i+1)+"_n_p",tmp_gain_time,tmp_gain_x,tmp_gain_y,tmp_gain_pairs])

                    self.delta_gain_track_info_dic["tk_"+str(i+1)+"_n_n"] = [ [] for n in range(5) ]
                    self.delta_gain_track_info_dic["tk_"+str(i+1)+"_n_p"] = [ [] for n in range(5) ]

                    break

        #
        # cal gain current
        #

        for i in range(len(self.gain_track_info_list)):

            self.charges = self.gain_track_info_list[i][4]

            self.track_name = self.gain_track_info_list[i][0]
            self.track_time = self.gain_track_info_list[i][1]
            self.track_x = self.gain_track_info_list[i][2]
            self.track_y = self.gain_track_info_list[i][3]
            self.track_charges = self.gain_track_info_list[i][4]
            self.track_current = 0.

            self.end_condition = 0
            while(self.end_condition==0):

                if(self.track_y>=(det.det_thin-1) or self.track_x>=(det.det_width-1)):
                     
                    self.end_condition=4

                else:
                    efx = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_x_value)
                    efy = fen.cal_point_field(self.track_x, self.track_y,fen.electric_field_y_value)
                    ef = np.array([efx,efy])
                    ef_value = np.linalg.norm(ef)*1e4
                    self.e_field = np.array([efx,efy])
                    wefx = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_x_value)
                    wefy = fen.cal_point_field(self.track_x, self.track_y,fen.weighting_electric_field_y_value)
                    wef = np.array([wefx,wefy])
                    wef_value = np.linalg.norm(wef)*1e4
                    self.we_field = np.array([wefx,wefy])
                                                   
                    self.drift_diffusion(det,fen)
                    # SR current
                    self.track_current = abs(self.charges*e0*self.drift_velocity*wef_value)

                    self.update_gain_track_info()
                    self.update_step(det)
                    self.update_end_condition()            


        det.positive_cu.Reset()
        det.negtive_cu.Reset()

        det.gain_positive_cu.Reset()
        det.gain_negtive_cu.Reset()

        det.gain_n_n_cu.Reset()
        det.gain_n_p_cu.Reset()
        det.gain_p_n_cu.Reset()
        det.gain_p_p_cu.Reset()

        det.sum_cu.Reset()

        temp_positive_cu = ROOT.TH1F("temp+","temp+",det.n_bin,0,det.t_end)
        temp_negitive_cu = ROOT.TH1F("temp-","temp-",det.n_bin,0,det.t_end)
        temp_sum_cu = ROOT.TH1F("temp_sum","temp_sum",det.n_bin,0,det.t_end)
        
        #
        # initial current
        #
        for i in range(len(track.track_position)):

            for j in range(len(self.delta_track_info_dic_p["tk_"+str(i+1)][0])):
                temp_positive_cu.Fill(self.delta_track_info_dic_p["tk_"+str(i+1)][0][j], self.delta_track_info_dic_p["tk_"+str(i+1)][4][j])

            for k in range(len(self.delta_track_info_dic_n["tk_"+str(i+1)][0])):
                temp_negitive_cu.Fill(self.delta_track_info_dic_n["tk_"+str(i+1)][0][k], self.delta_track_info_dic_n["tk_"+str(i+1)][4][k])

            det.positive_cu.Add(temp_positive_cu)
            det.negtive_cu.Add(temp_negitive_cu)

            temp_positive_cu.Reset()
            temp_negitive_cu.Reset()


        #
        # gain current
        #

        temp_gain_cu = ROOT.TH1F("temp_gain","temp_gain",det.n_bin,0,det.t_end)

        for i in range(len(self.gain_track_info_list)):

            tmp_track_name = self.gain_track_info_list[i][0]

            for j in range(len(self.delta_gain_track_info_dic[tmp_track_name][0])):
                
                temp_gain_cu.Fill(self.delta_gain_track_info_dic[tmp_track_name][0][j],self.delta_gain_track_info_dic[tmp_track_name][4][j])

            if(tmp_track_name[-1]=='p'):
                det.gain_positive_cu.Add(temp_gain_cu)

            if(tmp_track_name[-1]=='n'):
                det.gain_negtive_cu.Add(temp_gain_cu)


            if(tmp_track_name[-3:]=='n_n'):
                det.gain_n_n_cu.Add(temp_gain_cu)

            if(tmp_track_name[-3:]=='n_p'):
                det.gain_n_p_cu.Add(temp_gain_cu)

            if(tmp_track_name[-3:]=='p_n'):
                det.gain_p_n_cu.Add(temp_gain_cu)

            if(tmp_track_name[-3:]=='p_p'):
                det.gain_p_p_cu.Add(temp_gain_cu)


            temp_gain_cu.Reset()


        det.sum_cu.Add(det.positive_cu)
        det.sum_cu.Add(det.negtive_cu)
        det.sum_cu.Add(det.gain_positive_cu)
        det.sum_cu.Add(det.gain_negtive_cu)

        #laudau_pairs = track.mips_laudau()
        laudau_pairs = 4000

        current_scale = laudau_pairs/track.ionized_total_pairs

        det.positive_cu.Scale(current_scale)
        det.negtive_cu.Scale(current_scale)
        det.gain_positive_cu.Scale(current_scale)
        det.gain_negtive_cu.Scale(current_scale)

        det.gain_n_n_cu.Scale(current_scale)
        det.gain_n_p_cu.Scale(current_scale)
        det.gain_p_n_cu.Scale(current_scale)
        det.gain_p_p_cu.Scale(current_scale)

        det.sum_cu.Scale(current_scale)
        
                            
        
    def draw_drift_path(self,det):
        # ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "canvas1", 200,10,1000, 1000)
        mg = ROOT.TMultiGraph("mg","")
        x_array=array('f')
        y_array=array('f')
        for i in range(len(self.delta_track_info_dic_p)):
            n=len(self.delta_track_info_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.delta_track_info_dic_p["tk_"+str(i+1)][1])
                y_array.extend(self.delta_track_info_dic_p["tk_"+str(i+1)][2])             
                gr_p = ROOT.TGraph(n,x_array,y_array)
                gr_p.SetMarkerColor(4)
                gr_p.SetLineColor(4)
                gr_p.SetLineStyle(1)
                mg.Add(gr_p)
                del x_array[:]
                del y_array[:]
        for j in range(len(self.delta_track_info_dic_n)):
            m=len(self.delta_track_info_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.delta_track_info_dic_n["tk_"+str(j+1)][1])
                y_array.extend(self.delta_track_info_dic_n["tk_"+str(j+1)][2])             
                gr_n = ROOT.TGraph(m,x_array,y_array)
                gr_n.SetMarkerColor(2)
                gr_n.SetLineColor(2)
                gr_n.SetLineStyle(1)
                mg.Add(gr_n)
                del x_array[:]
                del y_array[:]
        mg.GetXaxis().SetRangeUser(0,det.det_width)
        mg.GetYaxis().SetRangeUser(0,det.det_thin)
        mg.Draw("APL")
        c1.SaveAs("drift_path.pdf")



class Amplifier:
    def CSA_amp(self,det,t_rise,t_fall,trans_imp):
        hist = ROOT.TH1F()
        hist.SetNameTitle("After Amplifier","After Amplifier")
        det.sum_cu.Copy(hist)
        max_num = hist.GetNbinsX()
        preamp_Q = [0.0]*max_num
        itot = [0.0]*max_num
        shaper_out_Q = [0.0]*max_num
        shaper_out_V = [0.0]*max_num
        qtot = 0.0
        sh_max = 0.0

        tau_rise = t_rise / 2.2*1e-9
        tau_fall = t_fall / 2.2*1e-9   
        if (tau_rise == tau_fall):
            tau_rise *= 0.9
        time_unit = (hist.GetXaxis().GetXmax()- hist.GetXaxis().GetXmin()) / hist.GetNbinsX()
        for i in range(max_num-1):
            if (i>0):
                preamp_Q[i] = 0.0
                itot[i] = hist.GetBinContent(i)
                preamp_Q[i] = itot[i] * time_unit
                qtot += itot[i] *time_unit
            elif (i==0):
                preamp_Q[i]==0.0
        for i in range(max_num-1):
            dif_shaper_Q = preamp_Q[i]
            if(dif_shaper_Q != 0):
                for j in range(max_num-i):
                    shaper_out_Q[i+j] += tau_fall/(tau_fall+tau_rise)*dif_shaper_Q \
                        *(math.exp(-j*time_unit/tau_fall)-math.exp(-j*time_unit/tau_rise))
            if (abs(shaper_out_Q[i]) > abs(sh_max)):
                sh_max = shaper_out_Q[i]
        for i in range(max_num):
            if(sh_max==0.): shaper_out_V[i] = 0
            else:
                shaper_out_V[i] = shaper_out_Q[i] * trans_imp * 1e15 *qtot / sh_max
            hist.SetBinContent(i,shaper_out_V[i])
        return qtot,hist

def TCT_draw_current(det,ele_current,qtot,drift):

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
    det.negtive_cu.GetYaxis().SetRangeUser(0,2000e-6)           ###12e-6

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
    det.gain_positive_cu.GetYaxis().SetRangeUser(0,2000e-6)      ###12e-6

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
    det.sum_cu.GetYaxis().SetRangeUser(0,2000e-6)               ###12e-6


    c.cd(1)
    c.GetPad(1).SetLeftMargin(0.12)
    c.GetPad(1).SetRightMargin(0.12)
    # c.GetPad(1).SetTopMargin(0.12)
    c.GetPad(1).SetBottomMargin(0.14)
    # det.sum_cu.GetXaxis().SetTitleOffset(1.2)


    # det.gain_positive_cu.Draw("HIST")
    # det.gain_negtive_cu.Draw("SAME HIST")
    det.negtive_cu.Draw("HIST")
    det.positive_cu.Draw("SAME HIST")
    c.Update()

    # rightmax = 1.1*ele_current.GetMinimum()
    # n_scale = ROOT.gPad.GetUymin() / rightmax
    # ele_current.Scale(n_scale)
    c.cd(2)
    c.GetPad(2).SetLeftMargin(0.12)
    c.GetPad(2).SetRightMargin(0.12)
    # c.GetPad(2).SetTopMargin(0.12)
    c.GetPad(2).SetBottomMargin(0.14)
    # det.sum_cu.GetXaxis().SetTitleOffset(1.2)
    det.sum_cu.Draw("HIST")
    #c.Update()
    
    #ele_current.Draw("HIST")


    ele_current.SetLineColor(6)
    det.sum_cu.SetLineWidth(2)
    det.positive_cu.SetLineWidth(2)
    det.negtive_cu.SetLineWidth(2)
    ele_current.SetLineWidth(2)

    ele_current.GetXaxis().SetTitleSize(0.05)
    ele_current.GetXaxis().SetLabelSize(0.04)
    ele_current.GetXaxis().SetNdivisions(510)
    # ele_current.GetYaxis().SetTitleOffset(1.1)
    ele_current.GetYaxis().SetTitleSize(0.05)
    ele_current.GetYaxis().SetLabelSize(0.04)
    ele_current.GetYaxis().SetNdivisions(505)
    ele_current.GetXaxis().CenterTitle()
    ele_current.GetYaxis().CenterTitle()
    ele_current.GetXaxis().SetTitle("Time [s]")
    ele_current.GetYaxis().SetTitle("Current [A]")
    ele_current.GetYaxis().SetRangeUser(0,4000)     ###120
 
    c.cd(3)
    c.GetPad(3).SetLeftMargin(0.12)
    c.GetPad(3).SetRightMargin(0.12)
    # c.GetPad(4).SetTopMargin(0.12)
    c.GetPad(3).SetBottomMargin(0.14)
    ele_current.Draw("HIST")
    #c.Update()


    # axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(),# rightmax, 0, 510, "+L")
    # axis.SetLineColor(6)
    # axis.SetTextColor(6)
    # axis.SetTextSize(0.02)
    # axis.SetLabelColor(6)
    # axis.SetLabelSize(0.02)
    # axis.SetTitle("Ampl [mV]")
    # axis.CenterTitle()
    # axis.Draw("same")
# 
    legend = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend.AddEntry(det.negtive_cu, "initial electron", "l")
    legend.AddEntry(det.positive_cu, "initial hole", "l")
    #legend.AddEntry(det.gain_negtive_cu, "gain electron", "l")
    #legend.AddEntry(det.gain_positive_cu, "gain hole", "l")
    # legend.AddEntry(det.sum_cu, "e+h", "l")
    # legend.AddEntry(ele_current, "electronics", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(20)

    c.cd(1)
    legend.Draw("same")

    #c.Update()
    #c.SaveAs("drift_current.pdf")
    c.SaveAs("drift_current_lgad_200V.png")

    charge_t=det.sum_cu.Integral() \
        * ((det.sum_cu.GetXaxis().GetXmax() \
        - det.sum_cu.GetXaxis().GetXmin()) \
        / det.sum_cu.GetNbinsX()) * 1e15
    
    print(charge_t)
    print(qtot*1e15)
    drift.draw_drift_path(det)

    elecu = [ []for n in range(ele_current.GetNbinsX())]
    for i in range(ele_current.GetNbinsX()):
        elecu[i].append(ele_current.GetBinContent(i))
    mincu = max(elecu)
    # velocity = 0
    # for i in range(8):
    #     vel_small = elecu[i][0]*0.05
    #     velocity = velocity + vel_small
    return(mincu)



def draw_current_contribution(det):

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
    det.gain_n_p_cu.GetYaxis().SetRangeUser(0,500e-6)    ###12e-6


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

    #c.Update()
    #c.SaveAs("drift_current.pdf")
    c.SaveAs("drift_current_contribution_lgad_200V.png")

class SPAGeneration():
    def __init__(self):
        ### for Si 1064 IR
        self.tau = 7.5e-12   #ps
        self.alfa = 987      #SI
        self.power = 1e-11   #J/s
        self.wavelength = 1.064e-6  #m
        self.widthBeamWaist = 2e-6  #m
        self.refractionIndex = 3.51

    def getWidthSquared(self,z):
        return (self.widthBeamWaist**2)*(1+((self.wavelength*z)/(np.pi* (self.widthBeamWaist**2)*self.refractionIndex))**2)
    def getWidth(self,z):
        return np.sqrt(self.getWidthSquared(z))
    def getIntensity(self,r,z,z_o = 0 ):
        widthSquared= self.getWidthSquared(z-z_o)
        intensity = ((2*self.power)/(np.pi*widthSquared))*np.exp((-2*(r**2)/(widthSquared)))*np.exp(-self.alfa*z)
        return intensity
    def getCarrierDensity(self,r,z,z_o = 0 ):
        I = self.getIntensity(r,z,z_o)
        return (self.alfa*I)/(3.6*1.60217657e-19)

class ProjGrid():
    def __init__(self,r_min,r_max,r_nBins,z_max,z_min,z_nBins,z_o):
        self.r_min = r_min
        self.r_max = r_max
        self.r_nBins = r_nBins
        self.z_min = z_min
        self.z_max = z_max
        self.z_nBins = z_nBins
        self.z_o = z_o
        self.rArray = np.linspace(r_min*1e-6,r_max*1e-6,r_nBins)
        self.zArray = np.linspace(z_min*1e-6,z_max*1e-6,z_nBins)
        self.r_step = abs(r_max*1e-6-r_min*1e-6)/r_nBins
        self.z_step = abs(z_max*1e-6-z_min*1e-6)/z_nBins
        self.xArray = self.rArray
        self.yArray = self.rArray
        self.x_step = self.r_step
        self.y_step = self.r_step

class NNocarrier():
    def __init__(self, r, rlen, zlen, rbin, zbin):
        self.z_o = 0   # z focus position in um
        self.rlen = rlen
        self.zlen = zlen
        self.r_min = -r # um
        self.r_max = -r+ rlen # um
        self.r_nBins = rbin
        self.z_min = 0. # um
        self.z_max = zlen # um
        self.z_nBins = zbin

    def Nocarrier(self):
        my_pro = ProjGrid(self.r_min, self.r_max, self.r_nBins, self.z_max, self.z_min, self.z_nBins, self.z_o)
        rGrid, zGrid = np.meshgrid(my_pro.rArray,my_pro.zArray)
        carriergeneration = SPAGeneration()
        CGrid = carriergeneration.getCarrierDensity(rGrid, zGrid, self.z_o*1e-6)
        xGrid = rGrid.copy()
        projGrid = CGrid.copy()

        for i_z in range(self.z_nBins):
            for i_r in range(self.r_nBins):
                x_value = xGrid[i_z, i_r]
                r_valueArray = np.sqrt(x_value*x_value+ my_pro.yArray*my_pro.yArray)
                # Get z value
                z_value = zGrid[i_z, i_r]
                z_valueArray = np.ones_like(r_valueArray)*z_value
                carr_den_projArray = carriergeneration.getCarrierDensity(r_valueArray, z_valueArray, self.z_o*1e-6)
                # Project sum and take into account integral step
                projGrid[i_z, i_r] = (carr_den_projArray.sum()*my_pro.y_step)/2
        projGrid = projGrid*my_pro.x_step*my_pro.z_step
        sumprojGrid = np.sum(projGrid)
        return(projGrid, sumprojGrid)


def main():
    args = sys.argv[1:]
    model = "PIN"

    rlen = 50
    zlen = 1300
    rbin = 50
    zbin = 130
    # define geometry
    my_det = R2dDetector(zlen,rlen) # [um]

    # set material
    my_sic_mat = Material('Si')
    my_det.set_mat(my_sic_mat)

    # set doping
    if model == "LGAD":
        # my_det.set_doping("(6.1942*(1e14)/sqrt(2*3.1415926)/0.13821*exp(-pow(x[1]-0.67,2))/2/0.13821/0.13821 + (1e13))*(-(1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))*1e-12") # [cm-3]
        my_det.set_doping("(2.2*(1e14)/sqrt(2*3.1415926)/0.13821*exp(-pow(x[1]-0.67,2))/2/0.13821/0.13821 + (1e13))*(-(1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))*1e-12") # [cm-3]
    if model == "PIN":
        my_det.set_doping("1e13*(-(1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))*1e-12")

    # set operating condition
    my_det.set_temperature(300) # [K]
    my_det.set_bias_voltage(-200) # [V]

    # mesh
    my_det.mesh(1, 1) # [um]

    # calculate electric field
    my_possion_solver = FenicsPossion(my_det)
    my_possion_solver.solve()
    my_possion_solver.draw()

    num = 50
    cvmin = ROOT.TCanvas("cvmin", "vmin", 100, 100, 800, 600)
    mg = ROOT.TMultiGraph("mg","")
    x_array=array('f')
    y_array=array('f')
    x = [[]for n in range(num+1)]
    y = [ ]

    for i in range(1):
        r = i*(50/num)
        carriergen = NNocarrier(25, rlen, zlen, rbin, zbin)    #Input the position of the laser
        nocarrier, sumprojGrid = carriergen.Nocarrier()
        my_track = Tracks()
        print("test1")
        my_track.tct_mips(zbin, rbin, zlen, rlen, nocarrier)
        # print(my_track.track_position)
        print("test2")
        my_drift = TCTDrifts(my_track)
        print("test3")
        my_drift.cal_current(my_track,my_possion_solver,my_det,nocarrier)
        print("test4")
        # print(my_drift.delta_track_info_dic_n['tk_80'][5])
        # print(my_drift.delta_track_info_dic_p)
        # print(my_drift.gain_track_info_list)

        # after the electronics
        my_electronics = Amplifier()
        qtot,ele_current=my_electronics.CSA_amp(my_det,t_rise=0.4,t_fall=0.2,trans_imp=10)
    
        my_drift.draw_drift_path(my_det)

        mincu = TCT_draw_current(my_det,ele_current,qtot,my_drift)
        vmin = [i * 50 for i in mincu]
        print(r)
        print(vmin)
        draw_current_contribution(my_det)

        y.append(vmin)
        x[i].append(r)
        x_array.extend(x[i])
        y_array.extend(vmin)
        n = 1
        gr = ROOT.TGraph(n,x_array,y_array)
        gr.SetMarkerColor(47)
        gr.SetMarkerStyle(42)
        gr.SetMarkerSize(3)
        mg.Add(gr)
        del x_array[:]
        del y_array[:]
    mg.Draw("APL")
    cvmin.SaveAs("vmax200V_lgad.pdf")

    outputfile = ROOT.TFile("dataout200V_lgad.root","RECREATE")
    tree = ROOT.TTree("tree","tree")
    v = array('f', [0.])
    z = array('f', [0.])
    tree.Branch("v", v, "v/F")
    tree.Branch("z", z, "z/F")
    for i in range(num+1):
        v[0] = y[i][0]
        z[0] = x[i][0]
        tree.Fill()
    tree.Write()
    outputfile.Close()



if __name__ == '__main__':

    #record run time
    starttime = time.time()
    main() 
    endtime=time.time()
    dtime = endtime-starttime
    print ("the process run time %.8s s" %dtime) 

