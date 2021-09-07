
#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

'''
Description: Raser parameter settings      
@Date       : 2021/09/02 09:57:00
@Author     : tanyuhang
@version    : 1.0
'''
import json

# Define all input parameters used in raser main process
class Setting:
    def __init__(self,parameters):
        """
        Description:
            1.Different functions detector(), fenics(), pygeant4() define 
            different class parameters.
            2.Parameter defined according input parameters.
        Parameters:
        ---------
        det_model : str
            Define the sensor models for simulation.eg. planar-3D plugin-3D
        _pardic : dictionaries
            Storage the input parameters
        steplength : float
            The length of  each step for e-h pairs drift 
        @Modify:
        ---------
            2021/09/02
        """
        self._pardic = {}
        self.input2dic(parameters)
        self.det_model = self._pardic['det_model']
        self.read_par(self._pardic['parfile'])
        self.scan_variation()

    def input2dic(self,parameters):
        " Transfer input list to dictinary"
        for par in parameters:
            name,_,value=par.rpartition('=')
            self._pardic[name]=value

    def read_par(self,jsonfile):
        "Read the setting.json file and save the input parametersin paras"
        with open(jsonfile) as f:
            dic_pars = json.load(f)
        for dic_par in dic_pars:
            if dic_par['name'] in self.det_model:
                self.steplength = float(dic_par['steplength'])
                paras =  dic_par
        for x in paras: 
            if self.is_number(paras[x]):          
                paras[x] = float(paras[x])
            else:
                paras[x] = paras[x]
        self.paras = paras

    @property
    def detector(self):
        """
        Description:
            Define differnet types detectors parameters. Like:
            planar3D, plugin3D
        Parameters:
        ---------
        lx,ly,lz : float
            Detector length, width and height
        doping : float
            Doping concentation should times 1e12 /um^3  
            -- N-type is positive (negetive volatge applied) 
            -- P-type is negetive (positive volatge applied)
        temp : float
            Tempareture
        e_ir : float
            Radius of electrode in 3D
        e_gap : float
            Spacing between the electrodes in 3D
        @Returns:
        ---------
            A dictionary containing all parameters used in detector  
        @Modify:
        ---------
            2021/09/02
        """
        p = self.paras
        if "planar3D" in self.det_model:
            detector = {'name':'planar3D', 'lx':p['lx'], 'ly':p['ly'], 
                        'lz':p['lz'], 'doping':p['doping'], 
                        'voltage':p['voltage'], 'temp':p['temp']
                        }
            
        if "plugin3D" in self.det_model:
            detector = {'name':'plugin3D', 'lx':p['lx'], 'ly':p['ly'], 
                        'lz':p['lz'], 'doping':p['doping'], 
                        'voltage':p['voltage'], 'temp':p['temp'], 
                        'e_ir':p['e_ir'], 'e_gap':p['e_gap']
                        }
        return detector

    @property
    def fenics(self):
        """
        Description:
            Define differnet fenics parameters
        Parameters:
        ---------
        mesh : int
            Mesh precision value, the bigger the higher the accuracy
        xyscale : int
            In plane detector, scale_xy is scaling sensor 50 times at x and 
            y axis, so the precision can improve 50 times in echo distance 

        @Returns:
        ---------
            A dictionary containing all parameters used in fenics  
        @Modify:
        ---------
            2021/09/02
        """
        p = self.paras
        if "planar3D" in self.det_model:
            fenics = {'name':'planar3D', 
                      'mesh':p['mesh'], "xyscale":p['xyscale']}
        if "plugin3D" in self.det_model:
            fenics = {'name':'plugin3D', 
                      'mesh':p['mesh'], "xyscale":p['xyscale']}
        return fenics

    @property
    def pygeant4(self):
        """
        Description:
            Define differnet geant4 parameters
        Parameters:
        ---------
        maxstep : float
            Simulate the step size of each step in Geant4
        par_in : list
            Incident particle position
        par_out : list
            Theoretical position of outgoing particles
        g4_vis : bool
            False: Do not display graphical interface of geant4 partivles
            True: Do not display graphical interface of geant4 partivles
        @Returns:
        ---------
            A dictionary containing all parameters used in geant4  
        @Modify:
        ---------
            2021/09/02
        """
        p = self.paras
        if "planar3D" in self.det_model:
            pygeant4 = {'name':'planar3D',
                        'maxstep':p['maxstep'], 'g4_vis':p['g4_vis'],
                        'par_in':[p['par_inx'], p['par_iny'], p['par_inz']], 
                        "par_out":[p['par_outx'], p['par_outy'], p['par_outz']],
                        }
        if "plugin3D" in self.det_model:
            pygeant4 = {'name':'plugin3D', 
                        'maxstep':p['maxstep'], 'g4_vis':p['g4_vis'],
                        'par_in':[p['par_inx'], p['par_iny'], p['par_inz']], 
                        "par_out":[p['par_outx'], p['par_outy'], p['par_outz']],
                        }
        return pygeant4

    @property
    def amplifer(self):
        """
        Description:
            Define differnet amplifers parameters
        Parameters:
        ---------
        maxstep : float
            Simulate the step size of each step in Geant4
        par_in : list
            Incident particle position
        par_out : list
            Theoretical position of outgoing particles
        g4_vis : bool
            False: Do not display graphical interface of geant4 partivles
            True: Do not display graphical interface of geant4 partivles
        @Returns:
        ---------
            A dictionary containing all parameters used in geant4  
        @Modify:
        ---------
            2021/09/02
        """
        p = self.paras
        CSA_par = {'name':'CSA_ampl', 't_rise':p['t_rise'], 
                   't_fall':p['t_fall'], 'trans_imp':p['trans_imp'], 
                   'CDet':p['CDet']
                  }
        BB_par = {'name':'BB_ampl', 'BBW':p['BBW'], 
                  'BBGain':p['BBGain'], 'BB_imp':p['BB_imp'],'OscBW':p['OscBW']
                 }
        return [CSA_par,BB_par]

    def scan_variation(self):
        " Define parameters of batch mode"
        if "3Dscan" in self.det_model:
            self.total_events = int(self._pardic['total_e'])
            self.intance_number = int(self._pardic['instan'])
            self.g4seed = self.intance_number * self.total_events
            self.output = self._pardic["output"]
        else:
            self.total_events = 30
            self.g4seed = 0 

    def is_number(self,s):
        "Define whethre the s is a number or not"
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