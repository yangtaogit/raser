#!/usr/bin/env python3
# -*- encoding: utf-8 -*-
'''
Description:  Define physical models for different materiales   
@Date       : 2021/09/06 18:46:00
@Author     : yangtao
@version    : 1.0
'''

""" Define Material """

import math

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



""" Define Mobility Model """

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



""" Define Avalanche Model """

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