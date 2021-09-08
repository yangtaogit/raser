# -*- encoding: utf-8 -*-
'''
@Description: Calculate the weighting potential and electric field wiht fenics      
@Date       : 2021/08/31 15:04:25
@Author     : tanyuhang
@version    : 1.0
'''

import fenics
import mshr
import sys
import numpy as np
import ROOT
from array import array


#Calculate the weighting potential and electric field
class FenicsCal:

    def __init__(self,my_d,fen_dic):
        self.p_electric = []
        self.w_p_electric = []
        self.det_model = fen_dic['name']
        self.fl_x=my_d.l_x/fen_dic['xyscale']  
        self.fl_y=my_d.l_y/fen_dic['xyscale']  
        self.fl_z=my_d.l_z
        self.tol = 1e-14
        m_sensor_box=self.fenics_space(my_d)  
        self.mesh3D = mshr.generate_mesh(m_sensor_box,fen_dic['mesh'])  
        self.V = fenics.FunctionSpace(self.mesh3D, 'P', 1)
        self.fenics_p_electric(my_d)
        self.fenics_p_w_electric(my_d)

    def fenics_space(self,my_d):
        """
        @description: 
            Define the fenics solver space 
        @param:
            None
        @Returns:
            Fenics Box structure
        @Modify:
            2021/08/31
        """
        if "plugin3D" in self.det_model:
            self.sensor_range_confirm(my_d)
            m_sensor =  mshr.Box(fenics.Point(self.sx_l,self.sy_l, 0), 
                                 fenics.Point(self.sx_r,self.sy_r, self.fl_z))
            for i in range(len(my_d.e_tr)):
                e_t_i = my_d.e_tr[i]
                elec_n=mshr.Cylinder(fenics.Point(e_t_i[0],e_t_i[1],e_t_i[3]), 
                                     fenics.Point(e_t_i[0],e_t_i[1],e_t_i[4]),
                                     e_t_i[2],e_t_i[2])
                m_sensor =m_sensor - elec_n 
        elif "planar3D" in self.det_model:
            m_sensor =  mshr.Box(fenics.Point(0, 0, 0), 
                                 fenics.Point(self.fl_x, self.fl_y, self.fl_z))
        else:
            print("sensor model is wrrong")
            sys.exit()
        return m_sensor            
        
    def sensor_range_confirm(self,my_d):
        """
        @description:
            confirm the sensor range 
            at x,y axias to avoide the the big sensor size
            which will lead to the mesh complicate         
        @param:
            xv_min - fenics sensor x left value
            xv_max - fenics sensor x right value
            yv_min - fenics sensor y left value
            yv_max - fenics sensor y right value          
        @Modify:
            2021/08/31
        """    
        xv_list=[]
        yv_list=[]
        rest_length=50 #um
        length=0
        for i in range(len(my_d.e_tr)):
            e_t_i = my_d.e_tr[i]
            xv_list.append(e_t_i[0])
            yv_list.append(e_t_i[1])
            ele_radius= e_t_i[2]
        while length == 0:
            xv_max = max(xv_list)+ele_radius+rest_length 
            xv_min = min(xv_list)-ele_radius-rest_length 
            yv_max = max(yv_list)+ele_radius+rest_length
            yv_min = min(yv_list)-ele_radius-rest_length 
            if xv_max >= yv_max:
                yv_max = xv_max
            else:
                xv_max = yv_max
            if xv_min <= yv_min:
                yv_min = xv_min
            else:
                xv_min = yv_min
            if (xv_max > my_d.l_x or xv_min <0 
               or yv_max > my_d.l_y or yv_min < 0):
                rest_length -= 1
            else:
                length=1
        self.sx_l=xv_min 
        self.sx_r=xv_max 
        self.sy_l=yv_min 
        self.sy_r=yv_max 

    def fenics_p_electric(self,my_d):    
        """
        @description:
            Solve poisson equation to get potential and electric field
        @Modify:
            2021/08/31
        """
        if  "plugin3D" in self.det_model:
            bc_l=[]
            bc_l = self.boundary_definition_3D(my_d,"Possion")          
        elif "planar3D" in self.det_model:
            bc_l = self.boundary_definition_2D(my_d,"Possion")

        u = fenics.TrialFunction(self.V)
        v = fenics.TestFunction(self.V)
        f = fenics.Constant(self.f_value(my_d))
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx
        # Compute solution
        self.u = fenics.Function(self.V)
        fenics.solve(a == L, self.u, bc_l,
                     solver_parameters=dict(linear_solver='gmres',
                                            preconditioner='ilu'))
        #calculate electric field
        W = fenics.VectorFunctionSpace(self.mesh3D, 'P', 1)
        self.E_field = fenics.project(fenics.as_vector((self.u.dx(0),
                                                        self.u.dx(1),
                                                        self.u.dx(2))),W)

    def fenics_p_w_electric(self,my_d):  
        """
        @description:
            Solve Laplace equation to 
            get weighting potential and weighting electric field
        @Modify:
            2021/08/31
        """
        if  "plugin3D" in self.det_model:
            bc_l = []
            bc_l = self.boundary_definition_3D(my_d,"Laplace")
        elif "planar3D" in self.det_model:
            bc_l = self.boundary_definition_2D(my_d,"Laplace")
        # Define variational problem
        u_w = fenics.TrialFunction(self.V)
        v_w = fenics.TestFunction(self.V)
        f_w = fenics.Constant(0)
        a_w = fenics.dot(fenics.grad(u_w), fenics.grad(v_w))*fenics.dx
        L_w = f_w*v_w*fenics.dx
        # Compute solution
        self.u_w = fenics.Function(self.V)
        fenics.solve(a_w == L_w, self.u_w, bc_l)

    def boundary_definition_3D(self,my_d,model):
        """
        @description:
            Get boundary definition of 3D detector with Possion and Laplace
        @Modify:
            2021/08/31
        """
        bc_l = []
        p_ele,n_ele=self.model_para(my_d,model)
        for i in range (len(my_d.e_tr)):
            e_i = my_d.e_tr[i]
            str_e = "x[0]>={e_0}-{e_2} && x[0]<={e_0}+"\
                    +"{e_2} && x[1]>={e_1}-{e_2} && "\
                    +"x[1]<={e_1}+{e_2} && x[2]>={e_3} \
                    && x[2]<={e_4} && on_boundary"
            elec_p = str_e.format(e_0=e_i[0], e_1=e_i[1],
                                  e_2=e_i[2], e_3=e_i[3],
                                  e_4=e_i[4])
            if e_i[5] == "p":
                bc = fenics.DirichletBC(self.V, p_ele, elec_p)
            else:
                bc = fenics.DirichletBC(self.V, n_ele, elec_p)
            bc_l.append(bc)
        return bc_l

    def boundary_definition_2D(self,my_d,model):
        """
        @description:
            Get boundary definition of 2D detector with Possion and Laplace
        @Modify:
            2021/08/31
        """
        p_ele,n_ele=self.model_para(my_d,model)
        u_D = fenics.Expression('x[2]<tol ? p_1:p_2',
                                degree=2, tol=1E-14,
                                p_1=p_ele, p_2=n_ele)
        def boundary(x, on_boundary):
            return abs(x[2])<self.tol or abs(x[2]-self.fl_z)<self.tol
        bc_l = fenics.DirichletBC(self.V, u_D, boundary)
        return bc_l
    
    def model_para(self,my_d,model):
        """
        @description:
            Choose Possion and Laplace parameters
        @Modify:
            2021/08/31
        """       
        bc_l=[]
        if model == "Possion":
            p_ele = my_d.v_voltage
            n_ele = 0.0
        elif model == "Laplace":
            p_ele = 0.0
            n_ele = 1.0
        else:
            print("The input fenics solver model is wrong")
        return p_ele,n_ele

    def f_value(self,my_d):
        """
        @description: 
            Cal f_value of Poisson equation
        @param:
            perm_mat -- Dielectric constant of using material
                     -- 11.7 Silicon
                     -- 9.76 Silicon carbide
        @Modify:
            2021/08/31
        """
        if my_d.mater == 0:
            perm_mat = 11.7  
        elif my_d.mater == 1:
            perm_mat = 9.76  
        else:
            print("material is wrong")            
        e0 = 1.60217733e-19
        perm0 = 8.854187817e-12   #F/m
        f_value = e0*my_d.d_neff*1e6/perm0/perm_mat
        return f_value
        
    def get_e_field(self,px,py,pz):
        """
        @description: 
            Get eletric field at the px,py,pz position
        @param:
            out_range -- out_range = False
                      -- Position (x,y,z) don't exit in sensor fenics range
        @reture:
            Eletric field along x,y,z direction
        @Modify:
            2021/08/31
        """
        out_range=self.judge_fenics_range(px,py,pz)
        if out_range:   
            x_value,y_value,z_value = 0,0,0
        else:
            scale_px=px%self.fl_x
            scale_py=py%self.fl_y          
            try:
                x_value,y_value,z_value = self.E_field(scale_px,scale_py,pz)
            except RuntimeError:
                x_value,y_value,z_value = 0,0,0
        return x_value,y_value,z_value

    def get_w_p(self,px,py,pz):
        """
        @description: 
            Get weighting potential at the (x,y,z) position
        @param:
            out_range -- out_range = False
                      -- Position (x,y,z) don't exit in sensor fenics range
        @reture:
            Get weighting potential at (x,y,z) position
        @Modify:
            2021/08/31
        """
        out_range=self.judge_fenics_range(px,py,pz)
        if out_range:   
            f_w_p = 1.0
        else:
            scale_px=px%self.fl_x
            scale_py=py%self.fl_y   
            try:
                f_w_p = self.u_w(scale_px,scale_py,pz)
            except RuntimeError:
                f_w_p = 0.0
        return f_w_p

    def get_potential(self,px,py,pz):
        """
        @description: 
            Get potential at the (x,y,z) position
        @param:
            out_range -- out_range = False
                      -- Position (x,y,z) don't exit in sensor fenics range
        @reture:
            Get potential at (x,y,z) position
        @Modify:
            2021/08/31
        """
        out_range=self.judge_fenics_range(px,py,pz)
        if out_range:   #
            f_w_p = 0
        else:
            scale_px=px%self.fl_x
            scale_py=py%self.fl_y   
            try:
                f_w_p = self.u(scale_px,scale_py,pz)
            except RuntimeError:
                f_w_p = 0.0
        return f_w_p

    def judge_fenics_range(self,px,py,pz):
        """
        @description: 
           Judge whether (x,y,z) position is in sensor fenics range
        @reture:
            False: not
            Ture:  in
        @Modify:
            2021/08/31
        """
        if "plugin3D" in self.det_model:
            if (px < self.sx_l or px > self.sx_r 
                   or py < self.sy_l or py > self.sy_r):
                out_range=True
            else:
                out_range=False
        elif "planar3D" in self.det_model:
            out_range=False
        return out_range
    
    def __del__(self):
        pass


class FenicsCal2D:

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

        self.solve()
        #self.draw()

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
       
        c = ROOT.TCanvas( 'c', 'c',500, 500 )
        c.SetGrid()
        c.SetLeftMargin(0.18)
        c.SetBottomMargin(0.18)
        # c.Divide(1,2)

        c.cd()
        g_e.Draw()
        c.Modified()
        c.Update()
        c.SaveAs("./fig/silicon_lgad_2D_electricfield_150V.pdf")
