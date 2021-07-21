"""

Author: Tao Yang <yangtao@ihep.ac.cn> 
Created [2021-07-07 Wed 08:47] 
Based on Raser C++ https://github.com/dt-np/raser


This program solves 2D Poisson's equation for SiC PIN & LGAD doping profile

    - div grad u(x,y) = f(x,y)

on the unit interval with source f given by

    f(x,y) = e0*Q(x,y)/perm0/perm

    PIN:
    f(x,y) = e0*1e13*1e6/perm0/perm_sic

    LGAD:
    f(x,y) = (6.1942*pow(10,16)/sqrt(2*3.1415926)/0.13821*exp(-pow(x-0.67,2))/2/0.13821/0.13821 + pow(10,13))*(-(1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))

and boundary conditions given by
    u(x,y) = 0 for x = 0
    u(x,y) = Bias_voltage for x = thin

"""


import fenics
import numpy as np
import math
import random
import time
import sys
import ROOT
from array import array
from raser import twoD_time 

""" Global Constant """
E0 = 1.60217733e-19 # C
PERM0 = 8.854187817e-12 #F/m



""" Define Material """

class Material:

    def __init__(self,mat_name):
        self.mat_name = mat_name

    def mat_database(self):

        # material par
        self.si_par_dict = {'Permittivity' : 11.5,\
                             'Avalanche': 'unknown',\
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
    def __init__(self,model_name):
        self.model_name = model_name

    def cal_mobility(self, det, position, charge, electric_field):

        x = position[0]
        y = position[1]
        T = det.temperature
        E = electric_field

        doping_expr = det.doping_epr
        doping_expr = doping_expr.replace("x[1]","y")
        doping_expr = doping_expr.replace("sqrt","math.sqrt")
        doping_expr = doping_expr.replace("exp","math.exp")
        #print(doping_expr)
        Neff = eval(doping_expr)

        # SiC mobility
        if(self.model_name == 'SiC'):
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
        if(self.model_name == 'Si'):
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




""" Define Detector Geometry """

class R2dDetector:

    def __init__(self,det_width,det_thin):
        self.det_width = det_width #um X-axis
        self.det_thin = det_thin #um Y-axis

        self.n_bin = 1000
        self.t_end = 3e-9
        self.positive_cu = ROOT.TH1F("charge+","Positive Current",self.n_bin,0,self.t_end)
        self.negtive_cu = ROOT.TH1F("charge-","Negative Current",self.n_bin,0,self.t_end)
        self.sum_cu = ROOT.TH1F("charge","Current",self.n_bin,0,self.t_end)
    
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

        # print(potential_value_2d)


    def cal_weighting_possion(self):

        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        # Create mesh and function space
        mesh = fenics.RectangleMesh(fenics.Point(0, 0), fenics.Point(width, thin), nx, ny)
        V = fenics.FunctionSpace(mesh, "P", 1)

        # Define boundary condition
        u_D = fenics.Expression('x[1] < tol? det_voltage : 0', degree = 2,tol = 1E-14,det_voltage = -100)

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
                tmp_ef = (self.potential_value_2d[j][i] - self.potential_value_2d[j][i+1])/x_step
                self.electric_field_x_position[j].append(tmp_xpos)
                self.electric_field_x_value[j].append(tmp_ef)

                # weighting field
                tmp_wxpos = 0.5*x_step*(2*i+1)
                tmp_wef = (self.weighting_potential_value_2d[j][i] - self.weighting_potential_value_2d[j][i+1])/x_step
                self.weighting_electric_field_x_position[j].append(tmp_wxpos)
                self.weighting_electric_field_x_value[j].append(tmp_wef)


        # y direction 
        for j in range(ny):
            for i in range(nx+1):

                # electric field
                tmp_ypos = 0.5*y_step*(2*j+1)
                tmp_ef = (self.potential_value_2d[j][i]- self.potential_value_2d[j+1][i])/y_step
                self.electric_field_y_position[j].append(tmp_ypos)
                self.electric_field_y_value[j].append(tmp_ef)

                # weighting field
                tmp_wypos = 0.5*y_step*(2*j+1)
                tmp_wef = (self.weighting_potential_value_2d[j][i] - self.weighting_potential_value_2d[j+1][i])/y_step
                self.weighting_electric_field_y_position[j].append(tmp_wypos)
                self.weighting_electric_field_y_value[j].append(tmp_wef)     



    def cal_point_field(self,px_point,py_point,input_value):

        width = self.det.det_width
        thin = self.det.det_thin
        
        nx = self.det.nx
        ny = self.det.ny

        x_step = width/nx
        y_step = thin/ny

        #Interpolation method 
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

        value_11=input_value[nx1_v][ny1_v]
        value_21=input_value[nx2_v][ny1_v]
        value_12=input_value[nx1_v][ny2_v]
        value_22=input_value[nx2_v][ny2_v]
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

        cutline = int(self.det.ny/2)

        # plot electric field at x = middle
        ep = array( 'd' )
        ev = array( 'd' )

        for i in range(self.det.ny):
            ep.append(self.electric_field_y_position[i][cutline])
            ev.append(self.electric_field_y_value[i][cutline])

        #print(ep)
        #print(ev)

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
# 
        # ax1 = plt.subplot(121)        
        # plt.plot( [depth[250] for depth in self.electric_field_y_position], [ef[250] for ef in self.electric_field_y_value])
        # plt.xlabel('depth [um]')
        # plt.ylabel('electric field [V/um]')
# 
        # ax2 = plt.subplot(122)
        # plt.plot( [wdepth[250] for wdepth in self.weighting_electric_field_y_position], [wef[250] for wef in self.weighting_electric_field_y_value])
        # plt.xlabel('depth [um]')
        # plt.ylabel('weighting electric field [V/um]')
        # plt.tight_layout()
        # plt.show()
        



""" Define Track """

class Tracks:
    def mips(self,p_entry,p_exit,n_div):
        self.p_entry = p_entry
        self.p_exit = p_exit
        self.p_tracks = [ [] for n in range(n_div-1) ]
        self.d_tracks = []
        p_x=p_entry[0]
        p_y=p_entry[1]
        # self.py_entry=[]
        for i in range(n_div-1):
            x_div_point = (p_exit[0]-p_entry[0])/n_div*i+p_entry[0]+(p_exit[0]-p_entry[0])/(2*n_div)
            y_div_point = (p_exit[1]-p_entry[1])/n_div*i+p_entry[1]+(p_exit[1]-p_entry[1])/(2*n_div)
            self.p_tracks[i].append(x_div_point)
            self.p_tracks[i].append(y_div_point)
            self.d_tracks.append(math.sqrt(math.pow(x_div_point-p_x,2)+math.pow(y_div_point-p_y,2)))
            p_x = x_div_point
            p_y = y_div_point



""" Define Detector Geometry """

class Drifts:

    def __init__(self,track):
        self.muhh=1650   #mobility related with the magnetic field (now silicon useless)
        self.muhe=310
        self.BB=np.array([0,0])
        self.sstep=0.1 #drift step
        self.kboltz=8.617385e-5 #eV/K
        self.max_drift_len=1e9 #maximum diftlenght [um]
        self.d_dic_n = {}
        self.d_dic_p = {}
        for n in range(len(track.p_tracks)):
            self.d_dic_n["tk_"+str(n+1)] = [ [] for n in range(4) ]
            self.d_dic_p["tk_"+str(n+1)] = [ [] for n in range(4) ] 

    def initial_parameter(self):
        
        self.end_cond=0
        self.d_time=0
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
        else:
            self.delta_x=0
            self.delta_y=0 

    def drift_v(self,det,fen):

        my_mobility = Mobility(det.mat_name)

        ex_delta_f = fen.cal_point_field(self.d_x+self.delta_x,self.d_y+self.delta_y,fen.electric_field_x_value)
        ey_delta_f = fen.cal_point_field(self.d_x+self.delta_x,self.d_y+self.delta_y,fen.electric_field_y_value)

        e_delta_f = np.array([ex_delta_f,ey_delta_f])

        pos = [self.d_x+self.delta_x,self.d_y+self.delta_y]

        aver_e=(np.linalg.norm(self.e_field)+np.linalg.norm(e_delta_f))/2*1e4

        self.v_drift=my_mobility.cal_mobility(det, pos, self.charg, aver_e)*aver_e

        #drift part
        if(self.v_drift==0):
            self.delta_x=0.0
            self.delta_y=0.0
            self.dif_x=0.0
            self.dif_y=0.0
            self.end_cond=9
        else:
            #off when the field gets large enough
            DiffOffField=8  # the silicon value?              
            if(np.linalg.norm(e_delta_f)<DiffOffField):
                self.s_time=self.sstep*1e-4/self.v_drift
                s_sigma=math.sqrt(2*self.kboltz*my_mobility.cal_mobility(det, pos, self.charg, aver_e)*det.temperature*self.s_time)
                self.dif_x=random.gauss(0,s_sigma)*1e4
                self.dif_y=random.gauss(0,s_sigma)*1e4          
            else:
                self.dif_x=0.0
                self.dif_y=0.0

    def drift_s_step(self,det):
        # x axis   
        if((self.d_x+self.delta_x+self.dif_x)>=det.det_width): 
            self.d_cx = det.det_width
        elif((self.d_x+self.delta_x+self.dif_x)<0):
            self.d_cx = 0
        else:
            self.d_cx = self.d_x+self.delta_x+self.dif_x
        # y axis
        if((self.d_y+self.delta_y+self.dif_y)>=det.det_thin): 
            self.d_cy = det.det_thin
        elif((self.d_y+self.delta_y+self.dif_y)<0):
            self.d_cy = 0
        else:
            self.d_cy = self.d_y+self.delta_y+self.dif_y

    def drift_end_condition(self):    
        if(self.wpot>(1-1e-5)):
            self.end_cond=1
        if(self.d_x<=0):
            self.end_cond=2
        if(self.d_y<=0):
            self.end_cond=4
        if(self.path_len>self.max_drift_len):
            self.end_cond=6
        if(self.n_step>10000):
            self.end_cond=7
            
    def save_inf_track(self):
        e0 = 1.60217733e-19
        if(self.charge<0):
            if(self.charg>0):
                self.d_dic_p["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_p["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_p["tk_"+str(self.n_track)][2].append(self.charge)
                self.d_dic_p["tk_"+str(self.n_track)][3].append(self.d_time)
            else:
                self.d_dic_n["tk_"+str(self.n_track)][0].append(self.d_x)
                self.d_dic_n["tk_"+str(self.n_track)][1].append(self.d_y)
                self.d_dic_n["tk_"+str(self.n_track)][2].append(self.charge)
                self.d_dic_n["tk_"+str(self.n_track)][3].append(self.d_time)

    def cal_current(self,det,track):

        det.positive_cu.Reset()
        det.negtive_cu.Reset()
        det.sum_cu.Reset()
        self.sum_p_current = []
        test_p = ROOT.TH1F("test+","test+",det.n_bin,0,det.t_end)
        test_n = ROOT.TH1F("test-","test-",det.n_bin,0,det.t_end)
        test_sum = ROOT.TH1F("test sum","test sum",det.n_bin,0,det.t_end)
        total_pairs = 0

        for j in range(len(track.p_tracks)):
            for i in range(len(self.d_dic_p["tk_"+str(j+1)][2])):
                test_p.Fill(self.d_dic_p["tk_"+str(j+1)][3][i],self.d_dic_p["tk_"+str(j+1)][2][i])
            test_p = Drifts.get_current_his(self,test_p)           
            for i in range(len(self.d_dic_n["tk_"+str(j+1)][2])):
                test_n.Fill(self.d_dic_n["tk_"+str(j+1)][3][i],self.d_dic_n["tk_"+str(j+1)][2][i])
            test_n = Drifts.get_current_his(self,test_n)
            n_pairs=Drifts.get_energy_loss(self,track.d_tracks[j])
            total_pairs+=n_pairs
            test_p.Scale(n_pairs)
            test_n.Scale(n_pairs)            
            det.positive_cu.Add(test_p)
            det.negtive_cu.Add(test_n)
            test_p.Reset()
            test_n.Reset()

        laudau_t_pairs = Drifts.get_laudau_dis(self,track)
        n_scale = laudau_t_pairs/total_pairs
        det.positive_cu.Scale(n_scale)
        det.negtive_cu.Scale(n_scale)
        det.sum_cu.Add(det.positive_cu)
        det.sum_cu.Add(det.negtive_cu)

    def get_current_his(self,histo):
        e0 = 1.60217733e-19
        hist = ROOT.TH1F()
        histo.Copy(hist)
        for i in range(hist.GetNbinsX()):
            histo.SetBinContent(i, hist.GetBinContent(i) \
                /((hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) \
                / hist.GetNbinsX())*e0)
        return histo

    def get_energy_loss(self,d_drift):
        #Drift distance: d_drift
        loss_energy = 8.4 #silicon carbide /um
        Myf = ROOT.TFile("SiC_5um.root")
        hist = Myf.Get("Edep_device")
        gRandom = ROOT.TRandom3(0)
        ran_energy = hist.GetRandom(gRandom)
        n_pairs = d_drift*ran_energy*1e6/(5*loss_energy)
        return n_pairs

    def get_laudau_dis(self,track):
        loss_energy = 8.4
        d_x = math.pow(track.p_exit[0]-track.p_entry[0],2)
        d_y = math.pow(track.p_exit[1]-track.p_entry[1],2)
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

    def ionized_drift(self,track,fen,det):       
        for i in range(len(track.p_tracks)):
            self.n_track=i+1
            for j in range(2):
                if (j==0):
                    self.charg=1 #hole
                if (j==1):
                    self.charg=-1 #electron
                Drifts.initial_parameter(self)
                #generated particles positions
                self.d_x=track.p_tracks[i][0]#initial position
                self.d_y=track.p_tracks[i][1]
                while (self.end_cond==0):
                    if (self.d_y>=(det.det_thin-1) or self.d_x>=(det.det_width-1)):
                        self.end_cond=3  
                    else:                                       
                        # field of the position
                        ex_field = fen.cal_point_field(self.d_x,self.d_y,fen.electric_field_x_value)
                        ey_field = fen.cal_point_field(self.d_x,self.d_y,fen.electric_field_y_value)
                        self.e_field = np.array([ex_field,ey_field])
                        #delta_poisiton
                        Drifts.delta_p(self)
                        #drift_position
                        Drifts.drift_v(self,det,fen)
                        #drift_next_posiiton
                        Drifts.drift_s_step(self,det)
                        #charge_collection
                        delta_Ew=fen.cal_point_field(self.d_cx,self.d_cy,fen.weighting_potential_value_2d)-fen.cal_point_field(self.d_x,self.d_y,fen.weighting_potential_value_2d)
                        self.charge=self.charg*delta_Ew
                        if(self.v_drift!=0):
                            self.d_time=self.d_time+self.s_time
                            self.path_len+=self.sstep
                        self.d_x=self.d_cx
                        self.d_y=self.d_cy
                        self.wpot=fen.cal_point_field(self.d_x,self.d_y,fen.weighting_potential_value_2d)
                        Drifts.save_inf_track(self)  
                        Drifts.drift_end_condition(self)                                                                         
                    self.n_step+=1

        Drifts.cal_current(self,det,track)

    def draw_drift_path(self,det):
        # ROOT.gStyle.SetOptStat(0)
        c1 = ROOT.TCanvas("c1", "canvas1", 200,10,1000, 1000)
        mg = ROOT.TMultiGraph("mg","")
        x_array=array('f')
        y_array=array('f')
        for i in range(len(self.d_dic_p)):
            n=len(self.d_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.d_dic_p["tk_"+str(i+1)][0])
                y_array.extend(self.d_dic_p["tk_"+str(i+1)][1])             
                gr_p = ROOT.TGraph(n,x_array,y_array)
                gr_p.SetMarkerColor(4)
                gr_p.SetLineColor(4)
                gr_p.SetLineStyle(1)
                mg.Add(gr_p)
                del x_array[:]
                del y_array[:]
        for j in range(len(self.d_dic_n)):
            m=len(self.d_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.d_dic_n["tk_"+str(j+1)][0])
                y_array.extend(self.d_dic_n["tk_"+str(j+1)][1])             
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

def draw_current(det,ele_current,qtot,drift):

    ROOT.gStyle.SetOptStat(0)
    c = ROOT.TCanvas("c", "canvas", 200,10,1000, 1000)
    c.SetLeftMargin(0.12)
    c.SetRightMargin(0.12)
    # c.SetTopMargin(0.12)
    c.SetBottomMargin(0.14)
    det.sum_cu.GetXaxis().SetTitleOffset(1.2)
    det.sum_cu.GetXaxis().SetTitleSize(0.05)
    det.sum_cu.GetXaxis().SetLabelSize(0.04)
    det.sum_cu.GetXaxis().SetNdivisions(510)
    det.sum_cu.GetYaxis().SetTitleOffset(1.1)
    det.sum_cu.GetYaxis().SetTitleSize(0.05)
    det.sum_cu.GetYaxis().SetLabelSize(0.04)
    det.sum_cu.GetYaxis().SetNdivisions(505)
    det.sum_cu.GetXaxis().CenterTitle()
    det.sum_cu.GetYaxis().CenterTitle()
    det.sum_cu.GetXaxis().SetTitle("Time [s]")
    det.sum_cu.GetYaxis().SetTitle("Current [A]")
    det.sum_cu.Draw("HIST")
    #det.positive_cu.Draw("SAME HIST")
    #det.negtive_cu.Draw("SAME HIST")
    c.Update()
    rightmax = 1.1*ele_current.GetMinimum()
    n_scale = ROOT.gPad.GetUymin() / rightmax
    ele_current.Scale(n_scale)
    ele_current.Draw("SAME HIST")
    det.sum_cu.SetLineColor(3)
    det.positive_cu.SetLineColor(2)
    det.negtive_cu.SetLineColor(4)
    ele_current.SetLineColor(6)
    det.sum_cu.SetLineWidth(2)
    det.positive_cu.SetLineWidth(2)
    det.negtive_cu.SetLineWidth(2)
    ele_current.SetLineWidth(2)    
    axis = ROOT.TGaxis(ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymin(), ROOT.gPad.GetUxmax(), ROOT.gPad.GetUymax(), rightmax, 0, 510, "+L")
    axis.SetLineColor(6)
    axis.SetTextColor(6)
    axis.SetTextSize(0.02)
    axis.SetLabelColor(6)
    axis.SetLabelSize(0.02)
    axis.SetTitle("Ampl [mV]")
    axis.CenterTitle()
    axis.Draw("same")

    legend = ROOT.TLegend(0.5, 0.3, 0.8, 0.6)
    legend.AddEntry(det.negtive_cu, "electron", "l")
    legend.AddEntry(det.positive_cu, "hole", "l")
    legend.AddEntry(det.sum_cu, "e+h", "l")
    legend.AddEntry(ele_current, "electronics", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(40)
    legend.Draw("same")

    c.Update()
    c.SaveAs("drift_current.pdf")
    
    charge_t=det.sum_cu.Integral() \
        * ((det.sum_cu.GetXaxis().GetXmax() \
        - det.sum_cu.GetXaxis().GetXmin()) \
        / det.sum_cu.GetNbinsX()) * 1e15
    print(charge_t)
    print(qtot*1e15)
    drift.draw_drift_path(det)


def main():
    args = sys.argv[1:]
    model = args[0]

    # define geometry
    my_det = R2dDetector(100,100) # [um]

    # set material
    my_sic_mat = Material('SiC')
    my_det.set_mat(my_sic_mat)

    # set doping
    if model == "LGAD":
        my_det.set_doping("(6.1942*(1e14)/sqrt(2*3.1415926)/0.13821*exp(-pow(x[1]-0.67,2))/2/0.13821/0.13821 + (1e13))*((1.60217733e-19)*(1e6)/(8.854187817e-12)/(9.76))*1e-12") # [cm-3]
    if model == "PIN":
        twoD_time()

    # set operating condition
    my_det.set_temperature(300) # [K]
    my_det.set_bias_voltage(-200) # [V]

    # mesh
    my_det.mesh(0.1, 0.1) # [um]

    # calculate electric field
    my_possion_solver = FenicsPossion(my_det)
    my_possion_solver.solve()
    my_possion_solver.draw()

    my_track = Tracks()
    my_track.mips([50,0],[50,my_det.det_thin],100)

    drift = Drifts(my_track)
    drift.ionized_drift(my_track,my_possion_solver,my_det)
    ### after the electronics
    my_electronics = Amplifier()
    qtot,ele_current=my_electronics.CSA_amp(my_det,t_rise=0.4,t_fall=0.2,trans_imp=10)
    draw_current(my_det,ele_current,qtot,drift)



if __name__ == '__main__':

    #record run time
    starttime = time.time()
    main() 
    endtime=time.time()
    dtime = endtime-starttime
    print ("the process run time %.8s s" %dtime) 

