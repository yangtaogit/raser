# import matplotlib.pyplot as plt
from array import array
import fenics
import mshr
import numpy as np
import math
import random
import ROOT
import time
import sys
import os
import geant4_pybind as g4b

FACTOR_SIZE = 1.0
FACTOR_UNIT = 1.0  # only do unit transformation   cann't work now due to geant4 ?
#define the detector parameter 
class R3dDetector:
    def __init__(self,l_x,l_y,l_z):
        self.l_x = l_x / FACTOR_UNIT*FACTOR_SIZE
        self.l_y = l_y / FACTOR_UNIT*FACTOR_SIZE
        self.l_z = l_z / FACTOR_UNIT*FACTOR_SIZE
    def set_para(self,doping,voltage,temperature):
        self.d_neff = doping #dopingX1e12 cm^-3
        self.v_voltage = voltage #Voltage
        self.temperature = temperature #Voltage
        self.n_bin = 1000
        self.t_end = 3.0e-9
        self.mater = 1    # 0 is silicon, 1 is silicon carbide
        self.positive_cu = ROOT.TH1F("charge+","Positive Current",self.n_bin,0,self.t_end)
        self.negtive_cu = ROOT.TH1F("charge-","Negative Current",self.n_bin,0,self.t_end)
        self.sum_cu = ROOT.TH1F("charge","Total Current",self.n_bin,0,self.t_end)
    def set_3D_electrode(self,e_ir,e_gap):
        e_r = e_ir / FACTOR_UNIT * FACTOR_SIZE
        e_int = e_gap / FACTOR_UNIT * FACTOR_SIZE
        e_t_y = R3dDetector.infor_ele(self,e_r,e_int)
        self.e_tr=[]
        self.e_t_1 = [self.l_x*0.5          ,self.l_y*0.5      ,e_r,0,self.l_z,"p"]
        self.e_t_2 = [self.l_x*0.5-e_int    ,self.l_y*0.5      ,e_r,0,self.l_z,"n"]
        self.e_t_3 = [self.l_x*0.5+e_int    ,self.l_y*0.5      ,e_r,0,self.l_z,"n"]
        self.e_t_4 = [self.l_x*0.5-e_int*0.5,self.l_y*0.5+e_t_y,e_r,0,self.l_z,"n"]
        self.e_t_5 = [self.l_x*0.5+e_int*0.5,self.l_y*0.5+e_t_y,e_r,0,self.l_z,"n"]
        self.e_t_6 = [self.l_x*0.5-e_int*0.5,self.l_y*0.5-e_t_y,e_r,0,self.l_z,"n"]
        self.e_t_7 = [self.l_x*0.5+e_int*0.5,self.l_y*0.5-e_t_y,e_r,0,self.l_z,"n"]
        for i in range(7):
           n_e = eval('self.e_t_' + str(i+1))
           self.e_tr.append(n_e)

    def infor_ele(self,e_r,e_int):
        e_x_gap = self.l_x - 2*e_r - 2*e_int
        if e_x_gap < 0:
            print("the electrode at x position is large than sensor length")
            sys.exit(0)
        e_t_y = math.sqrt(e_int*e_int*0.75)
        if 2*e_t_y > self.l_y:
            print("the electrode at y position is large than sensor length")
            sys.exit(0)            
        return e_t_y

#Calculate the weighting potential and electric field
class Fenics_cal:
    #parameter of SiC
    def __init__(self,my_d,sensor_model,mesh_v):
        self.p_electric = []
        self.w_p_electric = []
        self.model = sensor_model
        if my_d.mater == 0:
            perm_sic = 11.7  #Permittivity Si
        elif my_d.mater == 1:
            perm_sic = 9.76  #Permittivity SiC
        else:
            print("material is wrong")
        e0 = 1.60217733e-19
        perm0 = 8.854187817e-12   #F/m
        self.f_value = e0*my_d.d_neff*1e6/perm0/perm_sic*FACTOR_UNIT*FACTOR_UNIT
        self.tol = 1e-14
        #fenics space        
        m_sensor =  mshr.Box(fenics.Point(0, 0, 0), fenics.Point(my_d.l_x, my_d.l_y, my_d.l_z))
        if self.model == "3D":
            for i in range(len(my_d.e_tr)):
                e_t_i = my_d.e_tr[i]
                elec_n=mshr.Cylinder(fenics.Point(e_t_i[0], e_t_i[1], e_t_i[3]), fenics.Point(e_t_i[0], e_t_i[1], e_t_i[4]),e_t_i[2],e_t_i[2])
                m_sensor =m_sensor - elec_n
        
        self.mesh3D = mshr.generate_mesh(m_sensor,mesh_v)      
        self.V = fenics.FunctionSpace(self.mesh3D, 'P', 1)

    def fenics_p_electric(self,my_d):    #get the electric potential
        # Define boundary condition
        if self.model == "3D":
            bc_u=[]
            for i in range (len(my_d.e_tr)):
                e_t_i = my_d.e_tr[i]
                str_e =  "x[0]>={elec_1_0}-{elec_1_2} && x[0]<={elec_1_0}+"+\
                    "{elec_1_2} && x[1]>={elec_1_1}-{elec_1_2} && "+\
                    "x[1]<={elec_1_1}+{elec_1_2} && x[2]>={elec_1_3} && x[2]<={elec_1_4} && on_boundary"
                elec_p = str_e.format(elec_1_0=e_t_i[0],elec_1_1=e_t_i[1],elec_1_2=e_t_i[2],elec_1_3=e_t_i[3],elec_1_4=e_t_i[4])
                if e_t_i[5] == "p":
                    bc = fenics.DirichletBC(self.V, my_d.v_voltage, elec_p)
                else:
                    bc = fenics.DirichletBC(self.V, 0.0, elec_p)
                bc_u.append(bc)
        else:
            u_D = fenics.Expression('x[2] < tol? det_voltage : 0', degree=2,tol=1E-14,det_voltage=my_d.v_voltage)
            def boundary(x, on_boundary):
                return abs(x[2])<self.tol or abs(x[2]-my_d.l_z)<self.tol
            bc_u = fenics.DirichletBC(self.V, u_D, boundary)
        # # Define variational problem
        u = fenics.TrialFunction(self.V)
        v = fenics.TestFunction(self.V)
        f = fenics.Constant(self.f_value)
        a = fenics.dot(fenics.grad(u), fenics.grad(v))*fenics.dx
        L = f*v*fenics.dx
        # # Compute solution
        self.u = fenics.Function(self.V)
        fenics.solve(a == L, self.u, bc_u,solver_parameters=dict(linear_solver='gmres', preconditioner='ilu'))
        W = fenics.VectorFunctionSpace(self.mesh3D, 'P', 1)
        self.E_field = fenics.project(fenics.as_vector((self.u.dx(0),self.u.dx(1),self.u.dx(2))),W)

    def fenics_p_w_electric(self,my_d):  #get the electric weighting potential
        #####Laplace's equation
        if self.model == "3D":
            bc_w=[]
            for i in range (len(my_d.e_tr)):
                e_t_i = my_d.e_tr[i]
                str_e =  "x[0]>={elec_1_0}-{elec_1_2} && x[0]<={elec_1_0}+"+\
                    "{elec_1_2} && x[1]>={elec_1_1}-{elec_1_2} && "+\
                    "x[1]<={elec_1_1}+{elec_1_2} && x[2]>={elec_1_3} && x[2]<={elec_1_4} && on_boundary"
                elec_p = str_e.format(elec_1_0=e_t_i[0],elec_1_1=e_t_i[1],elec_1_2=e_t_i[2],elec_1_3=e_t_i[3],elec_1_4=e_t_i[4])
                if e_t_i[5] == "p":
                    bc = fenics.DirichletBC(self.V, 0.0, elec_p)
                else:
                    bc = fenics.DirichletBC(self.V, 1.0, elec_p)
                bc_w.append(bc)
        else:
            u_w_D = fenics.Expression('x[2] < tol? 0 : 1', degree=2,tol=1E-14)
            def boundary_w(x, on_boundary):
                return abs(x[2])<self.tol or abs(x[2]-my_d.l_z)<self.tol
            bc_w = fenics.DirichletBC(self.V, u_w_D, boundary_w)            
        # # Define variational problem
        u_w = fenics.TrialFunction(self.V)
        v_w = fenics.TestFunction(self.V)
        f_w = fenics.Constant(0)
        a_w = fenics.dot(fenics.grad(u_w), fenics.grad(v_w))*fenics.dx
        L_w = f_w*v_w*fenics.dx
        # # Compute solution
        self.u_w = fenics.Function(self.V)
        fenics.solve(a_w == L_w, self.u_w, bc_w)

    def get_e_field(self,px,py,pz):

        try:
            x_value,y_value,z_value = self.E_field(px,py,pz)
        except RuntimeError:
            x_value,y_value,z_value = 0,0,0
        return x_value,y_value,z_value

    def get_w_p(self,px,py,pz):
        try:
            f_w_p = self.u_w(px,py,pz)
        except RuntimeError:
            f_w_p = 0.0
        return f_w_p
#Geant4 for particle drift path
class MyDetectorConstruction(g4b.G4VUserDetectorConstruction):
    "My Detector Construction"
    def __init__(self,my_detector,maxStep):
        g4b.G4VUserDetectorConstruction.__init__(self)
        self.solid = {}
        self.logical = {}
        self.physical = {}
        self.checkOverlaps = True
        self.create_world(my_detector)
        self.create_sic_box(
                        name = "Device",
                        sidex = my_detector.l_x*g4b.um,
                        sidey = my_detector.l_y*g4b.um,
                        sidez = my_detector.l_z*g4b.um,
                        translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,my_detector.l_z/2.0*g4b.um],
                        material_si = "Si",
                        material_c = "C",
                        colour = [0.,0.1,0.8],
                        mother = 'world')
        self.maxStep = maxStep
        self.fStepLimit = g4b.G4UserLimits(self.maxStep)
        self.logical["Device"].SetUserLimits(self.fStepLimit)

    def create_world(self,my_d):

        self.nist = g4b.G4NistManager.Instance()
        material = self.nist.FindOrBuildMaterial("G4_AIR")  
        self.solid['world'] = g4b.G4Box("world", my_d.l_x*4.0*g4b.um,  my_d.l_y*4.0*g4b.um,  my_d.l_z*4.0*g4b.um)
        self.logical['world'] = g4b.G4LogicalVolume(self.solid['world'], 
                                                material, 
                                                "world")
        self.physical['world'] = g4b.G4PVPlacement(None, g4b.G4ThreeVector(0,0,0), 
                                               self.logical['world'], 
                                               "world", None, False, 0,self.checkOverlaps)
        visual = g4b.G4VisAttributes()
        visual.SetVisibility(False)
        self.logical['world'].SetVisAttributes(visual)

    def create_sic_box(self, **kwargs):
        name = kwargs['name']
        material_si = self.nist.FindOrBuildElement(kwargs['material_si'],False)
        material_c = self.nist.FindOrBuildElement(kwargs['material_c'],False)
        sic_density = 3.2*g4b.g/g4b.cm3
        self.SiC = g4b.G4Material("SiC",sic_density,2) 
        self.SiC.AddElement(material_si,50*g4b.perCent)
        self.SiC.AddElement(material_c,50*g4b.perCent)
        translation = g4b.G4ThreeVector(*kwargs['translation'])
        visual = g4b.G4VisAttributes(g4b.G4Color(0.,0.1,0.8))
        mother = self.physical[kwargs['mother']]
        self.sidex = kwargs['sidex']
        self.sidey = kwargs['sidey']
        self.sidez = kwargs['sidez']

        self.solid["Device"] = g4b.G4Box("Device", self.sidex/2., self.sidey/2., self.sidez/2.)
        
        self.logical["Device"] = g4b.G4LogicalVolume(self.solid[name], 
                                                self.SiC, 
                                                "Device")
        self.physical["Device"] = g4b.G4PVPlacement(None,translation,                                                
                                            "Device",self.logical["Device"],
                                            mother, False, 0,self.checkOverlaps)
        self.logical["Device"].SetVisAttributes(visual)
        
    def Construct(self): # return the world volume
        self.fStepLimit.SetMaxAllowedStep(self.maxStep)
        return self.physical['world']
        
class MyPrimaryGeneratorAction(g4b.G4VUserPrimaryGeneratorAction):
    "My Primary Generator Action"
    def __init__(self,par_position,par_direction):
        g4b.G4VUserPrimaryGeneratorAction.__init__(self)     
        particle_table = g4b.G4ParticleTable.GetParticleTable()
        electron = particle_table.FindParticle("e-")      # define the beta electron
        beam = g4b.G4ParticleGun(1)
        beam.SetParticleEnergy(2.28*g4b.MeV)
        beam.SetParticleMomentumDirection(g4b.G4ThreeVector(par_direction[0],par_direction[1],par_direction[2]))
        beam.SetParticleDefinition(electron)
        beam.SetParticlePosition(g4b.G4ThreeVector(par_position[0]*g4b.um,par_position[1]*g4b.um,par_position[2]*g4b.um))  
        self.particleGun = beam

    def GeneratePrimaries(self, event):
        self.particleGun.GeneratePrimaryVertex(event)
# ------------------------------------------------------------------
class MyRunAction(g4b.G4UserRunAction):
    def __init__(self):
        g4b.G4UserRunAction.__init__(self)
        milligray = 1.e-3*g4b.gray
        microgray = 1.e-6*g4b.gray
        nanogray = 1.e-9*g4b.gray
        picogray = 1.e-12*g4b.gray

        g4b.G4UnitDefinition("milligray", "milliGy", "Dose", milligray)
        g4b.G4UnitDefinition("microgray", "microGy", "Dose", microgray)
        g4b.G4UnitDefinition("nanogray", "nanoGy", "Dose", nanogray)
        g4b.G4UnitDefinition("picogray", "picoGy", "Dose", picogray)
        self.eventIDs=[]
        self.edep_devices=[]
        self.p_steps=[]
        self.energy_steps=[]

    def BeginOfRunAction(self, run):
        g4b.G4RunManager.GetRunManager().SetRandomNumberStore(False)

    def EndOfRunAction(self, run):
        nofEvents = run.GetNumberOfEvent()
        if nofEvents == 0:
            print("nofEvents=0")
            return
# ------------------------------------------------------------------
class MyEventAction(g4b.G4UserEventAction):
    "My Event Action"
    def __init__(self, runAction):
        g4b.G4UserEventAction.__init__(self)
        self.fRunAction = runAction

    def BeginOfEventAction(self, event):
        self.edep_device=0.
        self.x_EdepC_device=0.
        self.y_EdepC_device=0.
        self.z_EdepC_device=0.
        self.p_step = []
        self.energy_step = []

    def EndOfEventAction(self, event):
        eventID = event.GetEventID()
        self.fRunAction.eventIDs.append(eventID)
        self.fRunAction.edep_devices.append(self.edep_device)
        self.fRunAction.p_steps.append(self.p_step)
        self.fRunAction.energy_steps.append(self.energy_step)

    def RecordDevice(self, edep,point_in,point_out):
        self.edep_device += edep
        self.x_EdepC_device += edep*(point_out.getX()+point_in.getX())*0.5
        self.y_EdepC_device += edep*(point_out.getY()+point_in.getY())*0.5
        self.z_EdepC_device += edep*(point_out.getZ()+point_in.getZ())*0.5
        self.p_step.append([point_in.getX()*1000,point_in.getY()*1000,point_in.getZ()*1000])
        self.energy_step.append(edep) 

# ------------------------------------------------------------------
class MySteppingAction(g4b.G4UserSteppingAction):
    "My Stepping Action"
    def __init__(self, eventAction):
        g4b.G4UserSteppingAction.__init__(self)
        self.fEventAction = eventAction

    def UserSteppingAction(self, step):
        edep = step.GetTotalEnergyDeposit()
        point_pre  = step.GetPreStepPoint()
        point_post = step.GetPostStepPoint() 
        point_in   = point_pre.GetPosition()
        point_out  = point_post.GetPosition()
        volume = step.GetPreStepPoint().GetTouchable().GetVolume().GetLogicalVolume()
        volume_name = volume.GetName()
        if(volume_name == "Device"):
            self.fEventAction.RecordDevice(edep,point_in,point_out)

class Geant4:
    def __init__(self,my_g4d,par_position,par_direction,seed,step_n):
        gRunManager = g4b.G4RunManagerFactory.CreateRunManager(g4b.G4RunManagerType.Default)
        rand_engine= g4b.RanecuEngine()
        g4b.HepRandom.setTheEngine(rand_engine)
        g4b.HepRandom.setTheSeed(seed)
        gRunManager.SetUserInitialization(my_g4d)	
        # set physics list
        physics_list =  g4b.FTFP_BERT()
        physics_list.SetVerboseLevel(1)
        physics_list.RegisterPhysics(g4b.G4StepLimiterPhysics())
        gRunManager.SetUserInitialization(physics_list)
        #define action
        myPGA = MyPrimaryGeneratorAction(par_position,par_direction)
        gRunManager.SetUserAction(myPGA)
        self.myRA= MyRunAction()                              # define run
        gRunManager.SetUserAction(self.myRA)
        self.myEA= MyEventAction(self.myRA)                   # define event
        gRunManager.SetUserAction(self.myEA)
        mySA= MySteppingAction(self.myEA)                     # define step
        gRunManager.SetUserAction(mySA)
        gRunManager.Initialize()
        gRunManager.BeamOn(step_n)
        del gRunManager          
# mobility model
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
    
#The drift of generated particles
class Drifts:
    def __init__(self,my_g4v,i):
        self.muhh=1650.0   #mobility related with the magnetic field (now silicon useless)
        self.muhe=310.0
        self.BB=np.array([0,0,0])
        self.sstep=1/FACTOR_UNIT*FACTOR_SIZE #drift step
        self.kboltz=8.617385e-5 #eV/K
        self.max_drift_len=1e9/FACTOR_UNIT #maximum diftlenght [um]
        self.d_dic_n = {}
        self.d_dic_p = {}
        self.beam_number = i
        self.tracks_p = my_g4v.myRA.p_steps[i]
        self.tracks_step_edep = my_g4v.myRA.energy_steps[i]
        self.tracks_t_edep = my_g4v.myRA.edep_devices[i]
        for n in range(len(self.tracks_p)-1):
            self.d_dic_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            self.d_dic_p["tk_"+str(n+1)] = [ [] for n in range(5) ]
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
            self.delta_z=-self.sstep*self.charg*FF[2]/np.linalg.norm(FF)
        else:
            self.delta_x=0
            self.delta_y=0
            self.delta_z=0

    def drift_v(self,my_d,my_f):
        e_delta_f = np.array(my_f.get_e_field(self.d_x+self.delta_x,self.d_y+self.delta_y,self.d_z+self.delta_z))
        aver_e=(np.linalg.norm(self.e_field)+np.linalg.norm(e_delta_f))/2.0*1e4/FACTOR_UNIT  #V/cm
        self.v_drift=sic_mobility(self.charg,aver_e,my_d)*aver_e  # mobility cm2/(V s) v : cm/s
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
            DiffOffField=100.0*FACTOR_UNIT  # the silicon value ???           
            if(np.linalg.norm(e_delta_f)<DiffOffField):
                self.s_time=self.sstep*1e-4*FACTOR_UNIT/self.v_drift
                s_sigma=math.sqrt(2.0*self.kboltz*sic_mobility(self.charg,aver_e,my_d)*my_d.temperature*self.s_time)
                self.dif_x=random.gauss(0.0,s_sigma)*1e4/FACTOR_UNIT
                self.dif_y=random.gauss(0.0,s_sigma)*1e4/FACTOR_UNIT
                self.dif_z=random.gauss(0.0,s_sigma)*1e4/FACTOR_UNIT          
            else:
                self.dif_x=0.0
                self.dif_y=0.0
                self.dif_z=0.0

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
    def save_inf_track(self):
        e0 = 1.60217733e-19
        if(self.charge<0):
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

    def cal_current(self,my_d,my_g4v):
        my_d.positive_cu.Reset()
        my_d.negtive_cu.Reset()
        my_d.sum_cu.Reset()
        self.sum_p_current = []
        test_p = ROOT.TH1F("test+","test+",my_d.n_bin,0,my_d.t_end)
        test_n = ROOT.TH1F("test-","test-",my_d.n_bin,0,my_d.t_end)
        test_sum = ROOT.TH1F("test sum","test sum",my_d.n_bin,0,my_d.t_end)
        total_pairs=0
        for j in range(len(self.tracks_p)-1):
            for i in range(len(self.d_dic_p["tk_"+str(j+1)][2])):
                test_p.Fill(self.d_dic_p["tk_"+str(j+1)][4][i],self.d_dic_p["tk_"+str(j+1)][3][i])
            test_p = Drifts.get_current_his(self,test_p)           
            for i in range(len(self.d_dic_n["tk_"+str(j+1)][2])):
                test_n.Fill(self.d_dic_n["tk_"+str(j+1)][4][i],self.d_dic_n["tk_"+str(j+1)][3][i])
            test_n = Drifts.get_current_his(self,test_n)
            if (my_d.mater == 1): # silicon carbide
                sic_loss_e = 8.4 #ev
            elif (my_d.mater == 0):   # silicon
                sic_loss_e = 3.6 #ev
            n_pairs=self.tracks_step_edep[j]*1e6/sic_loss_e
            total_pairs+=n_pairs
            test_p.Scale(n_pairs)
            test_n.Scale(n_pairs)            
            my_d.positive_cu.Add(test_p)
            my_d.negtive_cu.Add(test_n)
            test_p.Reset()
            test_n.Reset()
        laudau_t_pairs = self.tracks_t_edep*1e6/sic_loss_e
        # print("laudau:%s"%laudau_t_pairs)
        n_scale = laudau_t_pairs/total_pairs
        my_d.positive_cu.Scale(n_scale)
        my_d.negtive_cu.Scale(n_scale)
        my_d.sum_cu.Add(my_d.positive_cu)
        my_d.sum_cu.Add(my_d.negtive_cu)
    def get_current_his(self,histo):
        e0 = 1.60217733e-19
        hist = ROOT.TH1F()
        histo.Copy(hist)
        for i in range(hist.GetNbinsX()):
            histo.SetBinContent(i, hist.GetBinContent(i) \
                /((hist.GetXaxis().GetXmax() - hist.GetXaxis().GetXmin()) \
                / hist.GetNbinsX())*e0)
        return histo

    def ionized_drift(self,my_g4v,my_f,my_d):       
        for i in range(len(self.tracks_p)-1):
            self.n_track=i+1
            for j in range(2):
                if (j==0):
                    self.charg=1 #hole
                if (j==1):
                    self.charg=-1 #electron
                Drifts.initial_parameter(self)
                #generated particles positions
                self.d_x=self.tracks_p[i+1][0]#initial position
                self.d_y=self.tracks_p[i+1][1]
                self.d_z=self.tracks_p[i+1][2] 
                while (self.end_cond==0):
                    if (self.d_y>=(my_d.l_y-1.0/FACTOR_UNIT) or self.d_x>=(my_d.l_x-1.0/FACTOR_UNIT) or self.d_z>=(my_d.l_z-1.0/FACTOR_UNIT)):
                        self.end_cond=3  
                    else:                                     
                        # field of the position
                        self.e_field = np.array(my_f.get_e_field(self.d_x,self.d_y,self.d_z))
                        #delta_poisiton
                        Drifts.delta_p(self)
                        #drift_position
                        Drifts.drift_v(self,my_d,my_f)
                        #drift_next_posiiton
                        Drifts.drift_s_step(self,my_d)
                        #charge_collection
                        delta_Ew=my_f.get_w_p(self.d_cx,self.d_cy,self.d_cz)-my_f.get_w_p(self.d_x,self.d_y,self.d_z)
                        self.charge=self.charg*delta_Ew
                        if(self.v_drift!=0):
                            self.d_time=self.d_time+self.sstep*1e-4*FACTOR_UNIT/self.v_drift
                            self.path_len+=self.sstep
                        self.d_x=self.d_cx
                        self.d_y=self.d_cy
                        self.d_z=self.d_cz
                        self.wpot=my_f.get_w_p(self.d_x,self.d_y,self.d_z)
                        Drifts.save_inf_track(self)  
                        Drifts.drift_end_condition(self)                                                                         
                    self.n_step+=1
        Drifts.cal_current(self,my_d,my_g4v)

    def draw_drift_path(self,my_d,my_f,n_e_v):
        ROOT.gStyle.SetOptStat(0)
        # # ROOT.gROOT.SetBatch(1)
        c1 = ROOT.TCanvas("c", "canvas1", 200,10,1000, 1000)
        nx_e = n_e_v[0]
        ny_e = n_e_v[1]
        nz_e = n_e_v[2]
        structrue = ROOT.TH3D("","",nx_e,0,my_d.l_x,ny_e,0,my_d.l_y,nz_e,0,my_d.l_z)
        for k in range(nz_e):
            for j in range (ny_e):
                for i in range(nx_e):
                    x_v = (i+1)*(my_d.l_x/nx_e)
                    y_v = (j+1)*(my_d.l_y/ny_e)
                    z_v = (k+1)*(my_d.l_z/nz_e)
                    try:
                        x_value,y_value,z_value = my_f.E_field(x_v,y_v,z_v)
                        structrue.SetBinContent(i+1,j+1,k+1,0)
                    except RuntimeError:
                        structrue.SetBinContent(i+1,j+1,k+1,2)
        structrue.Draw("ISO")

        x_array=array('f')
        y_array=array('f')
        z_array=array('f')
        for i in range(len(self.d_dic_p)):
            n=len(self.d_dic_p["tk_"+str(i+1)][0])
            if(n>0):
                x_array.extend(self.d_dic_p["tk_"+str(i+1)][0])
                y_array.extend(self.d_dic_p["tk_"+str(i+1)][1]) 
                z_array.extend(self.d_dic_p["tk_"+str(i+1)][2])              
                gr_p = ROOT.TPolyLine3D(n,x_array,y_array,z_array)
                gr_p.SetLineColor(2)
                gr_p.SetLineStyle(1)
                gr_p.Draw("SAME")
                del x_array[:]
                del y_array[:]
                del z_array[:]
        for j in range(len(self.d_dic_n)):
            m=len(self.d_dic_n["tk_"+str(j+1)][0])
            if(m>0):
                x_array.extend(self.d_dic_n["tk_"+str(j+1)][0])
                y_array.extend(self.d_dic_n["tk_"+str(j+1)][1])
                z_array.extend(self.d_dic_n["tk_"+str(j+1)][2])                
                gr_n = ROOT.TPolyLine3D(m,x_array,y_array,z_array)
                gr_n.SetLineColor(4)
                gr_n.SetLineStyle(1)
                gr_n.Draw("SAME")
                del x_array[:]
                del y_array[:]
                del z_array[:]
        c1.SaveAs("fig/drift_path.root")
        del c1

class Amplifier:
    def CSA_amp(self,my_d,t_rise,t_fall,trans_imp):
        
        hist = ROOT.TH1F()
        my_d.sum_cu.Copy(hist)
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
            if sh_max == 0.0:
                shaper_out_V[i] = 0.0
            else:
                shaper_out_V[i] = shaper_out_Q[i] * trans_imp * 1e15 *qtot / sh_max
            hist.SetBinContent(i,shaper_out_V[i])

        charge_t=my_d.sum_cu.Integral() \
            * ((my_d.sum_cu.GetXaxis().GetXmax() \
            - my_d.sum_cu.GetXaxis().GetXmin()) \
            / my_d.sum_cu.GetNbinsX()) * 1e15
        print("total_charge:%s fc"%(qtot*1e15))
        print("total_charge:%s fc"%charge_t)
        return charge_t,qtot,hist

def draw_plot(my_detector,ele_current,qtot,my_drift,my_field,drift_path):

    
    c = ROOT.TCanvas("c", "canvas", 200,10,1000, 1000)
    c.SetLeftMargin(0.12)
    # c.SetTopMargin(0.12)
    c.SetBottomMargin(0.14)
    ROOT.gStyle.SetOptStat(0)
    my_detector.sum_cu.GetXaxis().SetTitleOffset(1.2)
    my_detector.sum_cu.GetXaxis().SetTitleSize(0.05)
    my_detector.sum_cu.GetXaxis().SetLabelSize(0.04)
    my_detector.sum_cu.GetXaxis().SetNdivisions(510)
    my_detector.sum_cu.GetYaxis().SetTitleOffset(1.1)
    my_detector.sum_cu.GetYaxis().SetTitleSize(0.05)
    my_detector.sum_cu.GetYaxis().SetLabelSize(0.04)
    my_detector.sum_cu.GetYaxis().SetNdivisions(505)
    my_detector.sum_cu.GetXaxis().CenterTitle()
    my_detector.sum_cu.GetYaxis().CenterTitle()
    my_detector.sum_cu.GetXaxis().SetTitle("Time [s]")
    my_detector.sum_cu.GetYaxis().SetTitle("Current [A]")
    my_detector.sum_cu.Draw("HIST")
    my_detector.positive_cu.Draw("SAME HIST")
    my_detector.negtive_cu.Draw("SAME HIST")
    c.Update()
    rightmax = 1.1*ele_current.GetMinimum()
    if rightmax == 0.0:
        n_scale = 0.0
    else:
        n_scale = ROOT.gPad.GetUymin() / rightmax
    ele_current.Scale(n_scale)
    ele_current.Draw("SAME HIST")
    my_detector.sum_cu.SetLineColor(3)
    my_detector.positive_cu.SetLineColor(2)
    my_detector.negtive_cu.SetLineColor(4)
    ele_current.SetLineColor(6)
    my_detector.sum_cu.SetLineWidth(2)
    my_detector.positive_cu.SetLineWidth(2)
    my_detector.negtive_cu.SetLineWidth(2)
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

    legend = ROOT.TLegend(0.5, 0.3, 0.9, 0.6)
    legend.AddEntry(my_detector.negtive_cu, "electron", "l")
    legend.AddEntry(my_detector.positive_cu, "hole", "l")
    legend.AddEntry(my_detector.sum_cu, "e+h", "l")
    legend.AddEntry(ele_current, "electronics", "l")
    legend.SetBorderSize(0)
    legend.SetTextFont(43)
    legend.SetTextSize(45)
    legend.Draw("same")

    c.Update()
    c.SaveAs("fig/basic_infor.root")
    del c

    if drift_path == 1:
        my_drift.draw_drift_path(my_detector,my_field,[100,100,100])

def draw_ele_field(my_d,my_f,plane,depth):

    c1 = ROOT.TCanvas("c", "canvas",1000, 1000)
    ROOT.gStyle.SetOptStat(ROOT.kFALSE)
    # ROOT.gStyle.SetOptFit()
    c1.SetLeftMargin(0.12)
    c1.SetRightMargin(0.2)
    c1.SetBottomMargin(0.14)
    c1.Divide(2,2)
    c1.GetPad(1).SetRightMargin(0.2)
    c1.GetPad(2).SetRightMargin(0.2)
    c1.GetPad(3).SetRightMargin(0.2)
    c1.cd(1)
    e_field=fill_his("E",depth,my_d,my_f,plane)
    # ROOT.gStyle.SetPalette(107)
    e_field.Draw("COLZ")
    # e_field.SetMaximum(1)
    # e_field.SetMinimum(0)    
    c1.Update()
    c1.cd(2)
    p_field=fill_his("P",depth,my_d,my_f,plane)
    p_field.Draw("COLZ")
    c1.SetRightMargin(0.12)
    c1.Update()
    c1.cd(3)
    w_p_field=fill_his("WP",depth,my_d,my_f,plane)
    w_p_field.Draw("COLZ")
    c1.SetRightMargin(0.12)
    c1.Update()
    c1.SaveAs("fig/ele_field"+plane+str(depth)+".root")
    del c1

def fill_his(model,depth,my_d,my_f,plane):

    if plane == "xy":
        l_x = my_d.l_x
        l_y = my_d.l_y
        t_name = plane + " at z = " + str(depth)
    elif plane == "yz":
        l_x = my_d.l_y
        l_y = my_d.l_z
        t_name = plane + " at x = " + str(depth)
    elif plane == "xz":
        l_x = my_d.l_x
        l_y = my_d.l_z
        t_name = plane + " at y = " + str(depth)
    else:
        print("the draw plane is not existing")
    nx_e =200
    ny_e =200
    e_v = ROOT.TH2F("","",nx_e,0,l_x,ny_e,0,l_y)
    if model == "E":
        v_sample = my_f.E_field
        e_v.SetTitle("electric field "+t_name)
    elif model == "P":
        v_sample = my_f.u
        e_v.SetTitle("potential "+t_name)
    elif model == "WP":
        v_sample = my_f.u_w  
        e_v.SetTitle("weigthing potential "+t_name) 
    
    # ax = plt.gca()
    for j in range (ny_e):
        for i in range(nx_e):
            x_v = (i+1)*(l_x/nx_e)
            y_v = (j+1)*(l_y/ny_e)
            try:
                if plane == "xy":
                    f_v = v_sample(x_v,y_v,depth)
                elif plane == "yz":
                    f_v = v_sample(depth,x_v,y_v)
                elif plane == "xz":
                    f_v = v_sample(x_v,depth,y_v)
                if model == "E":
                    # ax.quiver(x_v,y_v,f_v[0],f_v[1])
                    f_v = math.sqrt(math.pow(f_v[0],2)+math.pow(f_v[1],2)+math.pow(f_v[2],2))                           
            except RuntimeError:
                f_v = 0.0
                # ax.quiver(x_v,y_v,0,0)
            e_v.SetBinContent(i+1,j+1,f_v)
    # plt.savefig("fig/test.png")
    # print("save fig/test.png")
    # plt.show()
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
def save_charge(charge_t,qtot,x_v,y_v,output_path):

    with open(output_path+'.txt','a') as f:
        f.write(str(x_v)+','+str(y_v)+','+str(charge_t)+','+str(qtot)+'\n')

def threeD_time(sensor_model):
    ### define the structure of the detector
    my_detector = R3dDetector(100.0,100.0,100.0)
    my_detector.set_para(doping=10.0,voltage=-500.0,temperature=300.0) #N-type is positive and P-type is negetive, doping /um^3 
    if sensor_model == "3D":
        my_detector.set_3D_electrode(5.0,40.0)
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector,sensor_model,mesh_v=32)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)

    ### Geant4 get drift path
    # par_position = [30.,50.,100.]
    # par_out = [30.,50.,0.]
    par_position = [30.,30.,100.]
    par_out = [30.,70.,0.]
    par_direction = [par_out[0]-par_position[0],par_out[1]-par_position[1],par_out[2]-par_position[2]]
    my_g4d = MyDetectorConstruction(my_detector,maxStep=1.0*g4b.um)
    seed = random.randint(0,100000)
    my_g4v = Geant4(my_g4d,par_position,par_direction,seed,step_n=1)
    # # # ### drift of ionized particles
    my_drift = Drifts(my_g4v,0)
    my_drift.ionized_drift(my_g4v,my_field,my_detector)

    # # ### after the electronics
    my_electronics = Amplifier()
    charge_t,qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
    # # ### electric plot
    draw_ele_field(my_detector,my_field,"xz",my_detector.l_z*0.5)
    # draw_ele_field(my_detector,my_field,"xy",my_detector.l_z*0.5)
    # draw_ele_field(my_detector,my_field,"yz",my_detector.l_z*0.5)
    # # ###  current plot
    draw_plot(my_detector,ele_current,qtot,my_drift,my_field,drift_path=1)

### get the 3D time resolution
def twoD_time_scan(output,numbers,t_numbers,step_n,change_para,sensor_model):
    ## define the structure of the detector
    print(change_para)
    my_detector = R3dDetector(100.,100.,100.)
    my_detector.set_para(doping=10,voltage=-500,temperature=300)
    #N-type is  positive and P-type is negetive, doping /um^3 
    if sensor_model == "3D_scan":
        my_detector.set_3D_electrode(5.0,change_para)
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector,sensor_model,mesh_v=32)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    par_position = [50.,50.,100.]
    par_out = [50.,50.,0.]
    par_direction = [par_out[0]-par_position[0],par_out[1]-par_position[1],par_out[2]-par_position[2]]
    my_g4d = MyDetectorConstruction(my_detector,maxStep=1.0*g4b.um)
    my_g4v = Geant4(my_g4d,par_position,par_direction,numbers-step_n,step_n)   
    for i in range(numbers-step_n,numbers): 
        print("event number:%s"%i)
        # # # ### drift of ionized particles
        my_drift = Drifts(my_g4v,i-numbers+step_n)
        my_drift.ionized_drift(my_g4v,my_field,my_detector)
        x_v = 50.
        y_v = 50.
        # # ### after the electronics
        my_electronics = Amplifier()
        charge_t,qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
        ROOT.gROOT.SetBatch(1)
        c = ROOT.TCanvas("Plots", "Plots", 1000, 1000)
        ele_current.Draw("HIST")
        c.Update()
        output_path = output + "/twoD_voltage_500"
        os.system("mkdir %s -p"%(output_path))
        c.SaveAs(output_path+"/t_"+str(i)+"_vx_"+str(int(x_v))+"_vy_"+str(int(y_v))+"_events.C")
        del c
        save_charge(charge_t,qtot,x_v,y_v,output_path)

### get the 3D time resolution
def threeD_time_scan(output,numbers,t_numbers,step_n,change_para,sensor_model):
    ## define the structure of the detector
    print(change_para)
    my_detector = R3dDetector(100.,100.,100.)
    my_detector.set_para(doping=10,voltage=-500,temperature=300)
    #N-type is  positive and P-type is negetive, doping /um^3 
    if sensor_model == "3D_scan":
        my_detector.set_3D_electrode(5.0,change_para)
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector,mesh_v=32)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    ### drift path
    x_min = my_detector.e_t_2[0]-my_detector.e_t_2[2]
    x_max = my_detector.e_t_3[0]+my_detector.e_t_3[2]
    y_min = my_detector.e_t_7[1]-my_detector.e_t_7[2]
    y_max = my_detector.e_t_4[1]+my_detector.e_t_4[2]
    nx = ny = int(math.sqrt(t_numbers))
    x_step = (x_max-x_min)/(nx-1)
    y_step = (y_max-y_min)/(ny-1)
    print("number:%s"%numbers)
    n_y = int(numbers/nx)
    n_x = numbers%nx
    x_v = x_min + x_step*n_x
    y_v = y_min + y_step*n_y
    par_position = [x_v,y_v,100.]
    par_out = [x_v,y_v,0.]
    par_direction = [par_out[0]-par_position[0],par_out[1]-par_position[1],par_out[2]-par_position[2]]
    my_g4d = MyDetectorConstruction(my_detector,maxStep=1.0*g4b.um)
    my_g4v = Geant4(my_g4d,par_position,par_direction,numbers,1) 
    ### drift of ionized particles
    my_drift = Drifts(my_g4v,0)
    my_drift.ionized_drift(my_g4v,my_field,my_detector)
#     ### after the electronics
    my_electronics = Amplifier()
    charge_t,qtot,ele_current=my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10)
    ROOT.gROOT.SetBatch(1)
    c = ROOT.TCanvas("Plots", "Plots", 1000, 1000)
    ele_current.Draw("HIST")
    c.Update()
    output_path = output + "/gap_"+str(change_para)+"_voltage_350"
    os.system("mkdir %s -p"%(output_path))
    c.SaveAs(output_path+"/t_"+str(numbers)+"_vx_"+str(int(x_v))+"_vy_"+str(int(y_v))+"_events.C")
    del c
    save_charge(charge_t,qtot,x_v,y_v,output_path)

def main():
    args = sys.argv[1:]
    model = args[0]
    if model == "2D":
        threeD_time("2D")
    if model == "3D":
        threeD_time("3D")
    if model == "3D_scan":
        output = args[1]
        numbers = int(args[2])
        t_numbers = int(args[3])
        n_step = int(args[4])
        change_para = float(args[5])
        threeD_time_scan(output,numbers,t_numbers,n_step,change_para,model)
    if model == "2D_scan":
        output = args[1]
        numbers = int(args[2])
        t_numbers = int(args[3])
        step_n = int(args[4])
        change_para = float(args[5])
        twoD_time_scan(output,numbers,t_numbers,step_n,change_para,model)
    print("run end")
    
if __name__ == '__main__':
    #record run time
    starttime = time.time()
    main() 
    endtime=time.time()
    dtime = endtime-starttime
    print ("the process run time %.8s s" %dtime) 