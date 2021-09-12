import math
import ROOT
import numpy as np
from raser.geometry import R2dDetector
from raser.geometry import R3dDetector

class MIPs:

    # MIPs particle
    def __init__(self,model,det_dic=None):
        
        if model == "NORMAL":
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

        elif model == "TCT":
            self.i_z = int(det_dic["z_nBins"])
            self.i_r = int(det_dic["r_nBins"])
            self.zlen = det_dic["det_width"]
            self.rlen = det_dic["det_thin"]
            self.z_o = det_dic["z_o"]
            self.tau = det_dic["tau"]
            self.alfa = det_dic["alfa"]
            self.power = det_dic["power"]
            self.wavelength = det_dic["wavelength"]
            self.widthBeamWaist = det_dic["widthBeamWaist"]
            self.refractionIndex = det_dic["refractionIndex"]

            self.track_position = [ [ [] for n in range(self.i_z)] for m in range(self.i_r)]
            for i in range(self.i_r):
                self.track_entry = [0,i]
                self.track_exit  = [self.zlen,i]
                track_position_x = self.track_entry[0]
                track_position_y = self.track_entry[1]
                for j in range(self.i_z):
                    x_div_point = j*(self.zlen/self.i_z)+ (self.zlen/(self.i_z*2))
                    y_div_point = i*(self.rlen/self.i_r)+ (self.rlen/(self.i_r*2))
                    self.track_position[i][j].append(x_div_point)
                    self.track_position[i][i].append(y_div_point)

        else:
            raise NameError(model)

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
        
        self.mips_ionized()
    
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


    def mips_landau(self):

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
    """
    Description:
        Transfer Carrier Distribution from Laser Coordinate System 
        to 2d Detector Coordinate System
    Parameters:
    ---------
    my_d : R2dDetector or R3dDetector
        the Detector
    laser : dict
        the Parameter List of Your Laser
    x_rel,y_rel,z_rel:
        the Normalized Coordinate for Laser Focus 
        in Detector Coordinate System
    @Modify:
    ---------
        2021/09/13
    """
    def __init__(self,my_d,laser,x_rel,y_rel,z_rel,min_carrier=0):
        self.tech=laser["tech"]
        self.direction=laser["direction"]

        self.refractionIndex=laser["refractionIndex"]

        self.wavelength=laser["wavelength"]*1e-3 #um
        self.tau=laser["tau"]
        self.power=laser["power"]
        self.widthBeamWaist=laser["widthBeamWaist"]#um
        
        self.r_step=laser["r_step"]#um
        self.h_step=laser["h_step"]#um

        if "l_Reyleigh" not in laser:
            self.l_Rayleigh = np.pi*self.widthBeamWaist**2*self.refractionIndex/self.wavelength
        else:
            self.l_Rayleigh = laser["l_Rayleigh"]#um

        if self.direction in ("top","bottom"):
            self.y_step=self.h_step
            self.z_step=self.x_step=self.r_step
            if self.direction == "top":
                def _getCarrierDensity(x,y,z):
                    return self.getCarrierDensity(y,x**2+z**2)
            if self.direction == "bottom":
                def _getCarrierDensity(x,y,z):
                    return self.getCarrierDensity(self.ly-y,x**2+z**2)

        elif self.direction == "edge":
            self.x_step=self.h_step
            self.z_step=self.y_step=self.r_step
            def _getCarrierDensity(x,y,z):
                return self.getCarrierDensity(x,y**2+z**2)
        else:
            raise NameError(self.direction)

        if isinstance(my_d,R2dDetector):
            self.lx=my_d.det_width#um
            self.ly=my_d.det_thin
            self.lz=my_d.det_width
        elif isinstance(my_d,R3dDetector):
            self.lx=my_d.l_x#um
            self.ly=my_d.l_y
            self.lz=my_d.l_z
        else:
            raise TypeError(my_d)

        if self.tech == "SPA":
            self.alpha=laser["alpha"]
            
        elif self.tech == "TPA":
            self.beta_2=laser["beta_2"]

        xArray = np.linspace(0.5*self.x_step,self.lx-0.5*self.x_step,int(self.lx/self.x_step))-self.lx*x_rel
        yArray = np.linspace(0.5*self.y_step,self.ly-0.5*self.y_step,int(self.ly/self.y_step))-self.ly*y_rel
        zArray = np.linspace(0.5*self.z_step,self.lz-0.5*self.z_step,int(self.lz/self.z_step))-self.lz*z_rel

        Y,X,Z=np.meshgrid(yArray,xArray,zArray) #Feature of numpy.meshgrid
        self.projGrid=_getCarrierDensity(X,Y,Z)\
            *self.x_step*self.y_step*self.z_step*1e-18

        if isinstance(my_d,R2dDetector):
            Y2d,X2d=np.meshgrid(yArray+self.ly*y_rel,xArray+self.lx*x_rel)
            self.projGrid=np.sum(self.projGrid,axis=2)
            self.track_position = list(np.transpose(np.array([list(np.ravel(X2d)),list(np.ravel(Y2d))])))
            self.ionized_pairs = list(np.ravel(self.projGrid))
            self.ionized_total_pairs = 0
            self.x_min,self.y_min = self.lx,self.ly
            self.x_max,self.y_max = 0,0
            for i in range(len(self.ionized_pairs)-1,-1,-1):
                if self.ionized_pairs[i]<=min_carrier:
                    del self.ionized_pairs[i]
                    del self.track_position[i]
                else:
                    self.ionized_total_pairs+=self.ionized_pairs[i]
                    self.x_min=min(self.x_min,self.track_position[i][0])
                    self.y_min=min(self.y_min,self.track_position[i][1])
                    self.x_max=max(self.x_max,self.track_position[i][0])
                    self.y_max=max(self.y_max,self.track_position[i][1])

        elif isinstance(my_d,R3dDetector):
            pass

        self.draw_nocarrier2D(x_rel,y_rel,z_rel,min_carrier)

    def getCarrierDensity(self,h,r2):
        widthSquared=(self.widthBeamWaist**2)*(1+(h/self.l_Rayleigh)**2)

        if self.tech=="SPA":
            intensity = ((2*self.power)/(np.pi*widthSquared*1e-12))*np.exp((-2*r2/(widthSquared)))*np.exp(-self.alpha*h*1e-6)
            e0 = 1.60217733e-19
            return self.alpha*intensity/(3.6*e0)
            
        elif self.tech=="TPA":
            k=(self.power**2)*8*np.log(2)/(self.tau*(np.pi**2.5)*(np.log(4))**0.5)
            intensity_squared = k*np.exp(-4*r2/widthSquared)/((widthSquared**2)*1e-24)
            h_Planck = 6.626*1e-34
            speedofLight = 2.998*1e8
            return self.beta_2*self.wavelength*1e-6*intensity_squared/(2*h_Planck*speedofLight)

    def draw_nocarrier2D(self,x_rel,y_rel,z_rel,min_carrier):
        c1 = ROOT.TCanvas("c1","canvas2",200,10,1000,1000)
        h = ROOT.TH2D("h","pairs of carrier generation",\
            int((self.x_max-self.x_min)/self.x_step)+1,self.x_min-0.5*self.x_step,self.x_max+0.5*self.x_step,\
            int((self.y_max-self.y_min)/self.y_step)+1,self.y_min-0.5*self.y_step,self.y_max+0.5*self.y_step)
        for i in range(len(self.track_position)):
            h.Fill(self.track_position[i][0], self.track_position[i][1],self.ionized_pairs[i])
        h.Draw()
        h.GetXaxis().SetTitle("Width [μm]")
        h.GetYaxis().SetTitle("Depth [μm]")
        c1.SaveAs("./fig/nocarrier_"\
            +str(round(x_rel,5))+"_"\
            +str(round(y_rel,5))+"_"\
            +str(round(z_rel,5))+"_"\
            +str(min_carrier)+".pdf")  
