import math
import numpy as np
import ROOT

""" Define Track """

class MIPs:

    # MIPs particle
    def __init__(self):

        track_entry = [25.0,0]
        track_exit = [25.0,50.0]
        n_div = 100

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

class TCTTracks():
    def __init__(self,my_d,laser):
        self.tech=laser["tech"]
        self.direction=laser["direction"]

        self.alpha=laser["alpha"]
        self.beta_2=laser["beta_2"]
        self.refractionIndex=laser["refractionIndex"]

        self.wavelength=laser["wavelength"]
        self.tau=laser["tau"]
        self.power=laser["power"]
        self.widthBeamWaist=laser["widthBeamWaist"]
        
        self.x_rel=laser["x_rel"]
        self.y_rel=laser["y_rel"]
        self.z_rel=laser["z_rel"]

        if "l_Reyleigh" not in laser:
            self.l_Rayleigh = np.pi*self.widthBeamWaist**2*self.refractionIndex/self.wavelength
        else:
            self.l_Rayleigh=laser["l_Rayleigh"]

        if self.tech =="SPA":
            x_min=0               #um
            x_max=my_d.det_width
            y_min=0
            y_max=my_d.det_thin
            z_min=0
            z_max=my_d.det_width

            r_step=1
            h_step=10
        elif self.tech =="TPA":
            '''
            gen_vol_rel = 5
            r_div=10
            h_div=50

            x_min = max(0,self.x_rel*my_d.det_width-gen_vol_rel*self.l_Rayleigh*1e6)
            x_max = min(my_d.det_width,self.x_rel*my_d.det_width+gen_vol_rel*self.l_Rayleigh*1e6)
            y_min = max(0,self.y_rel*my_d.det_thin-gen_vol_rel*self.widthBeamWaist*1e6)
            y_max = min(my_d.det_thin,self.y_rel*my_d.det_thin+gen_vol_rel*self.widthBeamWaist*1e6)       
            z_min = max(0,self.z_rel*my_d.det_width-gen_vol_rel*self.widthBeamWaist*1e6)
            z_max = min(my_d.det_width,self.z_rel*my_d.det_width+gen_vol_rel*self.widthBeamWaist*1e6)

            r_step = 2*gen_vol_rel*self.widthBeamWaist*1e6/r_div
            h_step = 2*gen_vol_rel*self.l_Rayleigh*1e6/h_div
            '''
            x_min=0               #um
            x_max=my_d.det_width
            y_min=0
            y_max=my_d.det_thin
            z_min=0
            z_max=my_d.det_width

            r_step=5
            h_step=100
            
        if self.direction in ("top","bottom"):
            y_step=h_step
            z_step=x_step=r_step
        elif self.direction == "edge":
            x_step=h_step
            z_step=y_step=r_step
        else:
            raise NameError(self.direction)

        xArray = np.arange(x_min+0.5*x_step,x_max-0.49*x_step,x_step)
        yArray = np.arange(y_min+0.5*y_step,y_max-0.49*y_step,y_step)
        zArray = np.arange(z_min+0.5*z_step,z_max-0.49*z_step,z_step)

        self.projGrid =np.zeros([len(xArray),len(yArray)])
        if self.direction in ("top","bottom"):
            for i in range(len(xArray)):
                for j in range(len(yArray)):
                    r2valueArray = ((zArray-self.z_rel*my_d.det_width)**2+(xArray[i]-self.x_rel*my_d.det_width)**2)*1e-12
                    if self.direction == "top":
                        h_j=yArray[j]*1e-6
                    elif self.direction == "bottom":
                        h_j=(my_d.thin-yArray[j])*1e-6
                    projArray_z_axis=self.getCarrierDensity(h_j,self.y_rel*my_d.det_thin*1e-6,r2valueArray)
                    self.projGrid[i][j]=projArray_z_axis.sum()*x_step*y_step*z_step*1e-18
        elif self.direction == "edge":
            for i in range(len(xArray)):
                for j in range(len(yArray)):
                    r2valueArray = ((zArray-self.z_rel*my_d.det_width)**2+(yArray[j]-self.y_rel*my_d.det_thin)**2)*1e-12
                    projArray_z_axis=self.getCarrierDensity(xArray[i]*1e-6,self.x_rel*my_d.det_width*1e-6,r2valueArray)
                    self.projGrid[i][j]=projArray_z_axis.sum()*x_step*y_step*z_step*1e-18

        self.track_position = [ [] for n in range (len (xArray)*len (yArray)) ]
        self.ionized_pairs = []
        for i in range(len (xArray)):
            for j in range(len (yArray)):
                x_div_point = xArray[i]
                y_div_point = yArray[j]
                self.track_position[i*len(yArray)+j].append(x_div_point)
                self.track_position[i*len(yArray)+j].append(y_div_point)
                self.ionized_pairs.append(self.projGrid[i][j])

    def getCarrierDensity(self,h,h_focus,r2):
        I = self.getIntensity(h,h_focus,r2)
        if self.tech=="SPA":
            e0 = 1.60217733e-19
            return self.alpha*I/(3.6*e0)
        elif self.tech=="TPA":
            h_Planck = 6.626*1e-34
            speedofLight = 2.998*1e8
            return self.beta_2*self.wavelength*I/(2*h_Planck*speedofLight)     

    def getIntensity(self,h,h_focus,r2):
        widthSquared= self.getWidthSquared(h-h_focus)
        if self.tech=="SPA":
            intensity = ((2*self.power)/(np.pi*widthSquared))*np.exp((-2*(r2)/(widthSquared)))*np.exp(-self.alpha*h)
            return intensity
        elif self.tech=="TPA":
            k=(self.power**2)*8*np.log(2)/(self.tau*(np.pi**2.5)*(np.log(4))**0.5)
            intensity_squared =k*np.exp(-4*r2/widthSquared)/widthSquared**2
            return intensity_squared

    def getWidthSquared(self,delta_h):
        return (self.widthBeamWaist**2)*(1+(delta_h/self.l_Rayleigh)**2)
    
    def draw_nocarrier(self,my_d,number=""):
        c1 = ROOT.TCanvas("c1","canvas2",200,10,1000,1000)
        h = ROOT.TH2D("h","pairs of carrier generation",len(self.projGrid),0,my_d.det_width,len(self.projGrid[0]),0,my_d.det_thin)
        for i in range(len(self.projGrid)):
            for j in range(len(self.projGrid[0])):
                h.Fill(self.track_position[i*len(self.projGrid[0])+j][0], self.track_position[i*len(self.projGrid[0])+j][1], self.projGrid[i][j])
        h.Draw()
        c1.SaveAs("nocarrier_"+number+".pdf")  

    def mips_ionized(self):
        self.ionized_pairs = []
        self.ionized_total_pairs = 0.

        for i in range(len(self.projGrid)):
            for j in range(len(self.projGrid[0])):
                n_pairs = self.projGrid[i][j]
                self.ionized_pairs.append(n_pairs)
                self.ionized_total_pairs += n_pairs