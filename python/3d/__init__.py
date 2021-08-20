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
import gc
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
        self.t_bin = 50e-12
        self.t_end = 6.0e-9
        self.t_start = 0
        self.n_bin = int((self.t_end-self.t_start)/self.t_bin)
        self.mater = 1    # 0 is silicon, 1 is silicon carbide
        self.positive_cu = ROOT.TH1F("charge+","Positive Current",self.n_bin,self.t_start,self.t_end)
        self.negtive_cu = ROOT.TH1F("charge-","Negative Current",self.n_bin,self.t_start,self.t_end)
        self.sum_cu = ROOT.TH1F("charge","Total Current",self.n_bin,self.t_start,self.t_end)
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
        self.create_AlorSi_box(
            name = "Al",
            sidex = my_detector.l_x*g4b.um,
            sidey = my_detector.l_y*g4b.um,
            sidez = 10*g4b.um,
            translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,15000*g4b.um],
            material_type = "G4_Al",
            colour = [1,0.1,0.8],
            mother = 'world'
            )
        self.create_AlorSi_box(
            name = "Si_main",
            sidex = 1300*g4b.um,
            sidey = 1300*g4b.um,
            sidez = 33*g4b.um,
            translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,10000*g4b.um],
            material_type = "G4_Si",
            colour = [1,1,1],
            mother = 'world'
            )
        self.create_AlorSi_box(
            name = "Si_sub",
            sidex = 1300*g4b.um,
            sidey = 1300*g4b.um,
            sidez = 300*g4b.um,
            translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,9833.5*g4b.um],
            material_type = "G4_Si",
            colour = [0,0,1],
            mother = 'world'
            )
        self.create_pcb_board(
            name = "pcb1",
            sidex = 20000*g4b.um,
            sidey = 20000*g4b.um,
            sidez = 1500*g4b.um,
            translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,8933.5*g4b.um],
            tub_radius = 500*g4b.um,
            tub_depth = 1500*g4b.um,
            material_Si = "Si",
            material_O = "O",
            colour = [0,0.5,0.8],   
            mother = 'world'
            )

        self.create_sic_box(
            name = "Device",
            sidex = my_detector.l_x*g4b.um,
            sidey = my_detector.l_y*g4b.um,
            sidez = my_detector.l_z*g4b.um,
            translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,my_detector.l_z/2.0*g4b.um],
            material_Si = "Si",
            material_c = "C",
            colour = [1,0,0],
            mother = 'world'
            )

        self.create_sic_box(
            name = "SiC_sub",
            sidex = my_detector.l_x*g4b.um,
            sidey = my_detector.l_y*g4b.um,
            sidez = 350.0*g4b.um,
            translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,-175.0*g4b.um],
            material_Si = "Si",
            material_c = "C",
            colour = [0,1,1],
            mother = 'world'
            )
        self.create_pcb_board(
            name = "pcb2",
            sidex = 20000*g4b.um,
            sidey = 20000*g4b.um,
            sidez = 1500*g4b.um,
            translation = [my_detector.l_x/2.0*g4b.um,my_detector.l_y/2.0*g4b.um,-1100*g4b.um],
            tub_radius = 500*g4b.um,
            tub_depth = 1500*g4b.um,
            material_Si = "Si",
            material_O = "O",
            colour = [0,0.5,0.8],   
            mother = 'world'
            )
        self.maxStep = maxStep
        self.fStepLimit = g4b.G4UserLimits(self.maxStep)
        self.logical["Device"].SetUserLimits(self.fStepLimit)

    def create_world(self,my_d):

        self.nist = g4b.G4NistManager.Instance()
        material = self.nist.FindOrBuildMaterial("G4_AIR")  
        self.solid['world'] = g4b.G4Box("world", 25000*g4b.um,  25000*g4b.um,  25000*g4b.um)
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
        material_si = self.nist.FindOrBuildElement(kwargs['material_Si'],False)
        material_c = self.nist.FindOrBuildElement(kwargs['material_c'],False)
        sic_density = 3.2*g4b.g/g4b.cm3
        SiC = g4b.G4Material("SiC",sic_density,2) 
        SiC.AddElement(material_si,50*g4b.perCent)
        SiC.AddElement(material_c,50*g4b.perCent)
        translation = g4b.G4ThreeVector(*kwargs['translation'])
        visual = g4b.G4VisAttributes(g4b.G4Color(*kwargs['colour']))
        mother = self.physical[kwargs['mother']]
        sidex = kwargs['sidex']
        sidey = kwargs['sidey']
        sidez = kwargs['sidez']

        self.solid[name] = g4b.G4Box(name, sidex/2., sidey/2., sidez/2.)
        
        self.logical[name] = g4b.G4LogicalVolume(self.solid[name], 
                                                SiC, 
                                                name)
        self.physical[name] = g4b.G4PVPlacement(None,translation,                                                
                                            name,self.logical[name],
                                            mother, False, 0,self.checkOverlaps)
        self.logical[name].SetVisAttributes(visual)

    def create_AlorSi_box(self, **kwargs):
        name = kwargs['name']
        material_type = self.nist.FindOrBuildMaterial(kwargs['material_type'],False)

        translation = g4b.G4ThreeVector(*kwargs['translation'])
        visual = g4b.G4VisAttributes(g4b.G4Color(*kwargs['colour']))
        mother = self.physical[kwargs['mother']]
        sidex = kwargs['sidex']
        sidey = kwargs['sidey']
        sidez = kwargs['sidez']

        self.solid[name] = g4b.G4Box(name, sidex/2., sidey/2., sidez/2.)
        
        self.logical[name] = g4b.G4LogicalVolume(self.solid[name], 
                                                material_type, 
                                                name)
        self.physical[name] = g4b.G4PVPlacement(None,translation,                                                
                                            name,self.logical[name],
                                            mother, False, 0,self.checkOverlaps)
        self.logical[name].SetVisAttributes(visual)     

    def create_pcb_board(self, **kwargs):
        name = kwargs['name']
        material_si = self.nist.FindOrBuildElement(kwargs['material_Si'],False)
        material_O = self.nist.FindOrBuildElement(kwargs['material_O'],False)
        sic_density = 2.2*g4b.g/g4b.cm3
        SiO2 = g4b.G4Material("SiO2",sic_density,2) 
        SiO2.AddElement(material_si,1)
        SiO2.AddElement(material_O,2)
        translation = g4b.G4ThreeVector(*kwargs['translation'])
        visual = g4b.G4VisAttributes(g4b.G4Color(*kwargs['colour']))
        mother = self.physical[kwargs['mother']]
        sidex = kwargs['sidex']
        sidey = kwargs['sidey']
        sidez = kwargs['sidez']
        tub_radius = kwargs['tub_radius']
        tub_depth = kwargs['tub_depth']

        self.solid[name+"box"] = g4b.G4Box(name+"box", sidex/2., sidey/2., sidez/2.)
        self.solid[name+"tub"] = g4b.G4Tubs(name+"tub", 0,tub_radius,tub_depth, 0,360*g4b.deg)
        self.solid[name] = g4b.G4SubtractionSolid(name,self.solid[name+"box"],self.solid[name+"tub"])

        self.logical[name] = g4b.G4LogicalVolume(self.solid[name], 
                                                SiO2, 
                                                name)
        self.physical[name] = g4b.G4PVPlacement(None,translation,                                                
                                            name,self.logical[name],
                                            mother, False, 0,self.checkOverlaps)
        self.logical[name].SetVisAttributes(visual)

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
        # beam.SetParticleEnergy(0.546*g4b.MeV)
        beam.SetParticleMomentumDirection(g4b.G4ThreeVector(par_direction[0],par_direction[1],par_direction[2]))
        beam.SetParticleDefinition(electron)
        beam.SetParticlePosition(g4b.G4ThreeVector(par_position[0]*g4b.um,par_position[1]*g4b.um,par_position[2]*g4b.um))  

        beam2 = g4b.G4ParticleGun(1)
        beam2.SetParticleEnergy(0.546*g4b.MeV)
        beam2.SetParticleMomentumDirection(g4b.G4ThreeVector(par_direction[0],par_direction[1],par_direction[2]))
        beam2.SetParticleDefinition(electron)
        beam2.SetParticlePosition(g4b.G4ThreeVector(par_position[0]*g4b.um,par_position[1]*g4b.um,par_position[2]*g4b.um))  
        self.particleGun = beam
        self.particleGun2 = beam2

    def GeneratePrimaries(self, event):
        self.particleGun.GeneratePrimaryVertex(event)
        self.particleGun2.GeneratePrimaryVertex(event)
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
      
    def BeginOfRunAction(self, run):
        g4b.G4RunManager.GetRunManager().SetRandomNumberStore(False)
        self.eventIDs=[]
        self.edep_devices=[]
        self.p_steps=[]
        self.energy_steps=[]  
    def EndOfRunAction(self, run):
        nofEvents = run.GetNumberOfEvent()
        # print(len((self.p_steps)))    
        # print(self.p_steps)
        if nofEvents == 0:
            print("nofEvents=0")
            return
    def Record_events(self,eventID,edep_device,p_step,energy_step):
        if(len(p_step)>0):
            self.eventIDs.append(eventID)
            self.edep_devices.append(edep_device)
            self.p_steps.append(p_step)
            self.energy_steps.append(energy_step)
        else:
            self.edep_devices.append(edep_device)
            self.p_steps.append([0,0,0])
            self.energy_steps.append([0])
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
        print("eventID:%s"%eventID)
        # print(self.p_step)
        #self.fRunAction.Record_events(eventID,self.edep_device,self.p_step,self.energy_step)
        save_geant4_events(eventID,self.edep_device,self.p_step,self.energy_step)

    def RecordDevice(self, edep,point_in,point_out):
        self.edep_device += edep
        self.x_EdepC_device += edep*(point_out.getX()+point_in.getX())*0.5
        self.y_EdepC_device += edep*(point_out.getY()+point_in.getY())*0.5
        self.z_EdepC_device += edep*(point_out.getZ()+point_in.getZ())*0.5
        #print(point_in.getX()*1000,point_in.getY()*1000,point_in.getZ()*1000)
        self.p_step.append([point_in.getX()*1000,point_in.getY()*1000,point_in.getZ()*1000])
        self.energy_step.append(edep) 

def save_geant4_events(eventID,edep_device,p_step,energy_step):
    if(len(p_step)>0):
        s_eventIDs.append(eventID)
        s_edep_devices.append(edep_device)
        s_p_steps.append(p_step)
        s_energy_steps.append(energy_step)
    else:
        s_eventIDs.append(eventID)
        s_edep_devices.append(edep_device)
        s_p_steps.append([[0,0,0]])
        s_energy_steps.append([0])

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

class MyActionInitialization(g4b.G4VUserActionInitialization):
    def __init__(self,par_position,par_direction):
        g4b.G4VUserActionInitialization.__init__(self)
        self.par_position = par_position
        self.par_direction = par_direction
    # def BuildForMaster(self):
    #     self.SetUserAction(MyRunAction())
    def Build(self):
        self.SetUserAction(MyPrimaryGeneratorAction(self.par_position,self.par_direction))
        # global myRA_action
        myRA_action = MyRunAction()
        self.SetUserAction(myRA_action)
        myEA = MyEventAction(myRA_action)
        self.SetUserAction(myEA)
        self.SetUserAction(MySteppingAction(myEA))

class Geant4:
    def geant_run(self,my_g4d,par_position,par_direction,seed,step_n,geant_vis):
        if geant_vis == "1": 
            ui = None
            ui = g4b.G4UIExecutive(len(sys.argv), sys.argv)
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
        # define golbal parameter
        global s_eventIDs,s_edep_devices,s_p_steps,s_energy_steps
        s_eventIDs,s_edep_devices,s_p_steps,s_energy_steps=[],[],[],[]

        #define action
        gRunManager.SetUserInitialization(MyActionInitialization(par_position,par_direction))
        if geant_vis == "1":    
            visManager = g4b.G4VisExecutive()
            visManager.Initialize()
            UImanager = g4b.G4UImanager.GetUIpointer()
            UImanager.ApplyCommand('/control/execute init_vis.mac')
        else:
            UImanager = g4b.G4UImanager.GetUIpointer()
            UImanager.ApplyCommand('/run/initialize')
            
        gRunManager.BeamOn(step_n)
        if geant_vis == "1":  
            ui.SessionStart()
        self.p_steps=s_p_steps
        self.energy_steps=s_energy_steps
        self.edep_devices=s_edep_devices

          
# # # mobility model
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
        self.sstep=0.1/FACTOR_UNIT*FACTOR_SIZE #drift step
        self.kboltz=8.617385e-5 #eV/K
        self.max_drift_len=1e9/FACTOR_UNIT #maximum diftlenght [um]
        self.d_dic_n = {}
        self.d_dic_p = {}
        if i == 0:
            for j in range(len(my_g4v.p_steps)):
                if len(my_g4v.p_steps[j])>5:
                    self.beam_number = j
                    self.tracks_p = my_g4v.p_steps[j]
                    self.tracks_step_edep = my_g4v.energy_steps[j]
                    self.tracks_t_edep = my_g4v.edep_devices[j]            
        else:
            self.beam_number = i
            self.tracks_p = my_g4v.p_steps[i]
            self.tracks_step_edep = my_g4v.energy_steps[i]
            self.tracks_t_edep = my_g4v.edep_devices[i]
        for n in range(len(self.tracks_p)-1):
            self.d_dic_n["tk_"+str(n+1)] = [ [] for n in range(5) ]
            self.d_dic_p["tk_"+str(n+1)] = [ [] for n in range(5) ]
    def initial_parameter(self):
        self.end_cond=0
        self.d_time=1.0e-9
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
        test_p = ROOT.TH1F("test+","test+",my_d.n_bin,my_d.t_start,my_d.t_end)
        test_n = ROOT.TH1F("test-","test-",my_d.n_bin,my_d.t_start,my_d.t_end)
        test_sum = ROOT.TH1F("test sum","test sum",my_d.n_bin,my_d.t_start,my_d.t_end)
        total_pairs=0
        if (my_d.mater == 1): # silicon carbide
            sic_loss_e = 8.4 #ev
        elif (my_d.mater == 0):   # silicon
            sic_loss_e = 3.6 #ev
        for j in range(len(self.tracks_p)-1):
            for i in range(len(self.d_dic_p["tk_"+str(j+1)][2])):
                test_p.Fill(self.d_dic_p["tk_"+str(j+1)][4][i],self.d_dic_p["tk_"+str(j+1)][3][i])
            test_p = Drifts.get_current_his(self,test_p)           
            for i in range(len(self.d_dic_n["tk_"+str(j+1)][2])):
                test_n.Fill(self.d_dic_n["tk_"+str(j+1)][4][i],self.d_dic_n["tk_"+str(j+1)][3][i])
            test_n = Drifts.get_current_his(self,test_n)
            n_pairs=self.tracks_step_edep[j]*1e6/sic_loss_e
            total_pairs+=n_pairs
            test_p.Scale(n_pairs)
            test_n.Scale(n_pairs)            
            my_d.positive_cu.Add(test_p)
            my_d.negtive_cu.Add(test_n)
            test_p.Reset()
            test_n.Reset()
        laudau_t_pairs = self.tracks_t_edep*1e6/sic_loss_e
        print("laudau:%s"%laudau_t_pairs)
        if total_pairs != 0:
            n_scale = laudau_t_pairs/total_pairs
        else:
            n_scale=0
        my_d.sum_cu.Add(my_d.positive_cu)
        my_d.sum_cu.Add(my_d.negtive_cu)
        my_d.sum_cu.Scale(n_scale)

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
                    if (self.d_y<=(1.0/FACTOR_UNIT) or self.d_x<=(1.0/FACTOR_UNIT) or self.d_z<=(1.0/FACTOR_UNIT)):
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
    def __init__(self,my_d,CSA_par,BB_par,mintstep="50e-12"):
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
        IintTime = max(2.0*(self.t_rise+self.t_fall)*1e-9/self.time_unit,3.0*self.tau_BBA/self.time_unit)
        self.IMaxSh = int(max_hist_num + IintTime)

        self.charge_t=my_d.sum_cu.Integral() \
            * ((my_d.sum_cu.GetXaxis().GetXmax() \
            - my_d.sum_cu.GetXaxis().GetXmin()) \
            / my_d.sum_cu.GetNbinsX()) * 1e15
        # print("total_charge:%s fc"%(self.qtot*1e15))
        # print("total_charge:%s fc"%self.charge_t)


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
                    shaper_out_Q[i+j] += tau_fall/(tau_fall+tau_rise)*dif_shaper_Q \
                        *(math.exp(-j*self.time_unit/tau_fall)-math.exp(-j*self.time_unit/tau_rise))                
            if (abs(shaper_out_Q[i]) > abs(sh_max)):
                sh_max = shaper_out_Q[i]                
        Ci = 3.5e-11  #fF
        Qfrac = 1.0/(1.0+self.CDet*1e-12/Ci)
        self.CSA_ele = ROOT.TH1F("electronics","electronics",IMaxSh,0,IMaxSh*self.time_unit)
        for i in range(IMaxSh):
            if sh_max == 0.0:
                shaper_out_V[i] = 0.0
            elif self.CDet_j == 0:
                shaper_out_V[i] = shaper_out_Q[i] * self.trans_imp * 1e15 * self.qtot * Qfrac / sh_max            
            elif self.CDet_j ==1:
                shaper_out_V[i] = shaper_out_Q[i] * self.trans_imp * 1e15 * self.qtot / sh_max
            self.CSA_ele.SetBinContent(i,shaper_out_V[i])

        max_CSA_height = min(shaper_out_V)
        time_t = shaper_out_V.index(max_CSA_height)
        # print("CSA_time=%s"%(time_t*self.time_unit))
        

        return self.CSA_ele


    def BB_amp(self):
  
        IMaxSh = self.IMaxSh
        preamp_Q     = [0.0]*IMaxSh
        Vout_scope   = [0.0]*IMaxSh
        Iout_BB_RC   = [0.0]*IMaxSh
        Iout_C50     = [0.0]*IMaxSh   
        BBGraph      = [0.0]*IMaxSh

        step=1
        self.BB_ele = ROOT.TH1F("electronics BB","electronics BB",IMaxSh,0,IMaxSh*self.time_unit)
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
                Iout_C50[i+j]   += (dif_shaper_Q)/self.tau_scope*math.exp(-j*self.time_unit/self.tau_scope)
                Iout_BB_RC[i+j] += (dif_shaper_Q)/self.tau_BBA*math.exp(-j*self.time_unit/self.tau_BBA)
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
        
        print("BB_time=%s"%(time_t*self.time_unit))

        return self.BB_ele

    def save_ele(self,number,x_v,y_v,output_path):
        output_file = output_path + "/t_" +str(number)+"_vx_"+str(int(x_v))+"_vy_"+str(int(y_v))+"_events.csv"
        f1 = open(output_file,"w")
        f1.write("charge:%s fc\n"%(self.qtot*1e15))
        f1.write("time[ns],CSA Amplitude [mV], BB Amplitude [mV] \n")
        for i in range(self.BB_ele.GetNbinsX()):
            f1.write("%s,%s,%s \n"%(i*self.time_unit,self.CSA_ele[i],self.BB_ele[i]))
        f1.close()
        return self.charge_t,self.qtot
    def save_ele_scan_CSA(self,t_rise,t_fall,output_path):
        output_file = output_path + "/t_"+"_t_rise_"+str(t_rise)+"_t_fall_"+str(t_fall)+"_events.csv"
        f1 = open(output_file,"w")
        f1.write("charge:%s fc\n"%(self.qtot*1e15))
        f1.write("time[ns],CSA Amplitude [mV], BB Amplitude [mV] \n")
        for i in range(self.BB_ele.GetNbinsX()):
            f1.write("%s,%s,%s \n"%(i*self.time_unit,self.CSA_ele[i],self.BB_ele[i]))
        f1.close()
        return self.charge_t,self.qtot

    def save_ele_scan_BB(self,CDet_BB_imp,BB_w,output_path):
        output_file = output_path + "/t_"+"BBimpCDet_"+str(CDet_BB_imp)+"_BBw_"+str(BB_w)+"_events.csv"
        f1 = open(output_file,"w")
        f1.write("charge:%s fc\n"%(self.qtot*1e15))
        f1.write("time[ns],CSA Amplitude [mV], BB Amplitude [mV] \n")
        for i in range(self.BB_ele.GetNbinsX()):
            f1.write("%s,%s,%s \n"%(i*self.time_unit,self.CSA_ele[i],self.BB_ele[i]))
        f1.close()
        return self.charge_t,self.qtot

def draw_plot(my_detector,ele_current,my_drift,my_field,my_g4v,drift_path):
    c=ROOT.TCanvas("c","canvas1",1000,1000)
    # c = ROOT.TCanvas()
    c.cd()
    c.Update()
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
    # energy distribution in sensor 

    c1=ROOT.TCanvas("c1","canvas1",1000,1000)
    h1 = ROOT.TH1F("Edep_device", "Energy deposition in SiC", 100, 0., 0.1)
    for i in range (len(my_g4v.edep_devices)):
        h1.Fill(my_g4v.edep_devices[i])
    g1 = ROOT.TF1("m1","landau",0,0.1)
    h1.Fit(g1,"S")
    print("MPV:%s"%g1.GetParameter(1))

    h1.Draw()
    c1.SaveAs("fig/dep_SiC_energy.root")


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

def threeD_time(sensor_model,geant_vis):
    ### define the structure of the detector
    my_detector = R3dDetector(5000.0,5000.0,100.0)
    my_detector.set_para(doping=10.0,voltage=-400.0,temperature=300.0) #N-type is positive and P-type is negetive, doping /um^3 
    if sensor_model == "3D":
        my_detector.set_3D_electrode(5.0,40.0)
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector,sensor_model,mesh_v=32)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)

    ### Geant4 get drift path
    par_position = [2500.,2500.,17000.]
    par_out = [2500.,2500.,0.]
    par_direction = [par_out[0]-par_position[0],par_out[1]-par_position[1],par_out[2]-par_position[2]]
    my_g4d = MyDetectorConstruction(my_detector,maxStep=0.5*g4b.um)
    seed = random.randint(0,10000)

    my_g4v = Geant4()
    my_g4v.geant_run(my_g4d,par_position,par_direction,seed,100,geant_vis)
    # # # # ### drift of ionized particles
    my_drift = Drifts(my_g4v,0)
    my_drift.ionized_drift(my_g4v,my_field,my_detector)
    CSA_par = {'t_rise':0.7,'t_fall':1.1,'trans_imp':38,'CDet':30}  #CSA sensor parameters
    BB_par = {'BBW':0.66,'BBGain':19500,'BB_imp':10,'OscBW':2}
    my_electronics = Amplifier(my_detector,CSA_par,BB_par)
    #my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10,CDet=66)
    CSAele_current=my_electronics.CSA_amp()
    BBele_current=my_electronics.BB_amp()
    #charge_t,qtot=my_electronics.save_ele_scan_BB(CDet_BB_imp,BB_w,output_path)

    # # ### electric plot
    #draw_ele_field(my_detector,my_field,"xz",my_detector.l_z*0.5)
    draw_ele_field(my_detector,my_field,"xy",my_detector.l_z*0.5)
    draw_ele_field(my_detector,my_field,"yz",my_detector.l_z*0.5)
    # ###  current plot
    print("Test")
    draw_plot(my_detector,BBele_current,my_drift,my_field,my_g4v,drift_path=1)
    del CSAele_current
    del BBele_current
    # # ### after the electronics
    # for i in range (50):
    #     for j in range (50):
    #         t_rise=0.1 + 0.1*i
    #         t_fall=0.1 + 0.1*j
    #         # CDet_BB_imp = 100 + i*200
    #         # BB_w = 0.16 + j*0.1
    #         # BB_imp = CDet_BB_imp/30.0
    #         # print("CDet_BB_imp=%s,BBw=%s"%(CDet_BB_imp,BB_w))
    #         CSA_par = {'t_rise':t_rise,'t_fall':t_fall,'trans_imp':38,'CDet':30}  #CSA sensor parameters
    #         BB_par = {'BBW':0.66,'BBGain':19500,'BB_imp':10,'OscBW':2}
    #         my_electronics = Amplifier(my_detector,CSA_par,BB_par)
    #         #my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10,CDet=66)
    #         CSAele_current=my_electronics.CSA_amp()
    #         BBele_current=my_electronics.BB_amp()
    #         output_path = "out/ele_scan_CSA/"
    #         os.system("mkdir %s -p"%(output_path))      
    #         charge_t,qtot=my_electronics.save_ele_scan_CSA(t_rise,t_fall,output_path)
    #         #charge_t,qtot=my_electronics.save_ele_scan_BB(CDet_BB_imp,BB_w,output_path)
    #         del CSAele_current
    #         del BBele_current
    # # ### electric plot
    #draw_ele_field(my_detector,my_field,"xz",my_detector.l_z*0.5)
    #draw_ele_field(my_detector,my_field,"xy",my_detector.l_z*0.5)
    # draw_ele_field(my_detector,my_field,"yz",my_detector.l_z*0.5)
    # # ###  current plot
    #print("Test")
    #draw_plot(my_detector,BBele_current,my_drift,my_field,my_g4v,drift_path=1)
### get the 2D time resolution
def twoD_time_scan(output,numbers,t_numbers,step_n,change_para,sensor_model):
    ## define the structure of the detector
    print(change_para)
    my_detector = R3dDetector(5000.,5000.,100.)
    my_detector.set_para(doping=10,voltage=change_para,temperature=300)
    #N-type is  positive and P-type is negetive, doping /um^3 
    if sensor_model == "3D_scan":
        my_detector.set_3D_electrode(5.0,40.0)
    ### get the electric field and weighting potential
    my_field = Fenics_cal(my_detector,sensor_model,mesh_v=32)
    my_field.fenics_p_electric(my_detector) 
    my_field.fenics_p_w_electric(my_detector)
    par_position = [2500.,2500.,17000.]
    par_out = [2500.,2500.,0.]
    par_direction = [par_out[0]-par_position[0],par_out[1]-par_position[1],par_out[2]-par_position[2]]
    my_g4d = MyDetectorConstruction(my_detector,maxStep=0.5*g4b.um)
    my_g4v = Geant4()
    my_g4v.geant_run(my_g4d,par_position,par_direction,numbers-step_n,step_n,0) 
    # print(my_g4v.myRA.p_steps)
    # print(len(my_g4v.myRA.p_steps))
    for i in range(numbers-step_n,numbers): 
        if len(my_g4v.p_steps[i-numbers+step_n])>5:
            print("event number:%s"%i)
            # # # ### drift of ionized particles
            my_drift = Drifts(my_g4v,i-numbers+step_n)
            my_drift.ionized_drift(my_g4v,my_field,my_detector)
            x_v = 2500.
            y_v = 2500.
            # # ### after the electronics
            CSA_par = {'t_rise':0.4,'t_fall':0.2,'trans_imp':38,'CDet':30}  #CSA sensor parameters
            BB_par = {'BBW':0.66,'BBGain':19500,'BB_imp':10,'OscBW':2}
            my_electronics = Amplifier(my_detector,CSA_par,BB_par)
            #my_electronics.CSA_amp(my_detector,t_rise=0.4,t_fall=0.2,trans_imp=10,CDet=66)
            CSAele_current=my_electronics.CSA_amp()
            BBele_current=my_electronics.BB_amp()
            output_path = output + "/voltage_"+str(abs(change_para))
            os.system("mkdir %s -p"%(output_path))      

            charge_t,qtot=my_electronics.save_ele(i,x_v,y_v,output_path)
            save_charge(charge_t,qtot,x_v,y_v,output_path)
            del CSAele_current
            del BBele_current

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
    my_g4v = Geant4(my_g4d,par_position,par_direction,numbers,1,0) 
    ### drift of ionized particles
    my_drift = Drifts(my_g4v,0)
    my_drift.ionized_drift(my_g4v,my_field,my_detector,1)
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

def main(args):
    #args = sys.argv[1:]
    model = args[0]
    if model == "2D":
        geant_vis = args[1]
        threeD_time("2D",geant_vis)
    if model == "3D":
        geant_vis = args[1]
        threeD_time("3D",geant_vis)
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
    #main() 
    main(sys.argv[1:]) 
    endtime=time.time()
    dtime = endtime-starttime
    print ("the process run time %.8s s" %dtime)
    os._exit(0)

