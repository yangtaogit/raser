
#!/usr/bin/env python3
# -*- encoding: utf-8 -*-

'''
Description: Raser parameter settings      
@Date       : 2021/09/02 09:57:00
@Author     : tanyuhang
@version    : 1.0
'''


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
        self.steplength = 10 #um
        self.scan_variation()


    def input2dic(self,parameters):
        " Transfer input list to dictinary"
        for par in parameters:
            name,_,value=par.rpartition('=')
            self._pardic[name]=value
            
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
        self.detector_planar3D = {'name':'planar3D', 
                             'lx':5000, 'ly':5000, 'lz':100,
                             'doping':10, 'voltage':-500, 'temp':300.0}
        self.detector_plugin3D = {'name':'plugin3D', 
                             'lx':10000, 'ly':10000, 'lz':350,
                             'doping':10,'voltage':-500,'temp':300.0,
                             'e_ir':51.0, 'e_gap':150.0}
        detector_list = []
        detector_list=self.get_dics_list(self.detector_planar3D, 
                                         self.detector_plugin3D)
        detector=self.get_model_dic(detector_list)
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

        fenics_planar3D = {'name':'planar3D', 
                            'mesh':32, "xyscale":50}
        fenics_plugin3D = {'name':'plugin3D', 
                           'mesh':32, "xyscale":1}
        fenics_list = []
        fenics_list = self.get_dics_list(fenics_planar3D, fenics_plugin3D)
        fenics=self.get_model_dic(fenics_list)
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
        pla=self.detector_planar3D
        plu=self.detector_plugin3D

        pygeant4_planar3D = {'name':'planar3D',
                            'maxstep':0.5, 'g4_vis':False,
                            'par_in':[pla['lx']/2.0, pla['ly']/2.0, 17000.], 
                            "par_out":[pla['lx']/2.0, pla['ly']/2.0, 0.0],
                           }
        pygeant4_plugin3D = {'name':'plugin3D', 
                            'maxstep':0.5, 'g4_vis':False,
                            'par_in':[plu['lx']/2.0, plu['ly']/2.0, 17000.], 
                            "par_out":[plu['lx']/2.0, plu['ly']/2.0, 0.0],
                           }
        pygeant4_list = []
        pygeant4_list = self.get_dics_list(pygeant4_planar3D, pygeant4_plugin3D)
        pygeant4=self.get_model_dic(pygeant4_list)
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
        CSA_par = {'name':'CSA_ampl', 't_rise':0.7, 
                   't_fall':1.1, 'trans_imp':38, 'CDet':30
                  }
        BB_par = {'name':'BB_ampl', 'BBW':0.66, 
                  'BBGain':19500, 'BB_imp':10,'OscBW':2
                 }

        return [CSA_par,BB_par]

    def get_model_dic(self,input_list):
        """ 
        Input dictionary list
        Return a dictionary named by your detector model
        """
        for det in input_list:
            if det['name'] in self.det_model:
                out_dic = det
        return out_dic
       
    def get_dics_list(self,*input_dics):
        """ 
        Input saome dictionaries 
        Return a list contain all dictionaries
        """ 
        out_list = []
        for input_dic in input_dics:
            out_list.append(input_dic)
        return out_list

    def scan_variation(self):
        " Define parameters of batch mode"
        if "3Dscan" in self.det_model:
            self.total_events = 10
            self.intance_number = 0
            self.g4seed = self.intance_number * self.total_events
            self.output = "pyraser/unittest/"
        else:
            self.total_events = 30
            self.g4seed = 0 

