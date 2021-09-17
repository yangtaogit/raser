import math
import ROOT
import numpy as np
from raser.geometry import R2dDetector
from raser.geometry import R3dDetector

""" Define Track """

'''x_rel, y_rel, z_rel : float
            the Normalized Laser Focus Position'''

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

        self.wavelength=laser["wavelength"]*1e-3 #um
        self.tau=laser["tau"]
        self.power=laser["power"]
        self.widthBeamWaist=laser["widthBeamWaist"]#um
        
        self.r_step=laser["r_step"]#um
        self.h_step=laser["h_step"]#um

        if isinstance(my_d,R2dDetector):
            self.lx=my_d.det_width#um
            self.ly=my_d.det_thin
            self.lz=my_d.det_width

        elif isinstance(my_d,R3dDetector):
            self.lx=my_d.l_x#um
            self.ly=my_d.l_y
            self.lz=my_d.l_z

        if "l_Reyleigh" not in laser:
            self.l_Rayleigh = np.pi*self.widthBeamWaist**2*self.refractionIndex/self.wavelength
        else:
            self.l_Rayleigh = laser["l_Rayleigh"]#um

        if self.direction in ("top","bottom"):
            y_step=self.h_step
            z_step=x_step=self.r_step
        elif self.direction == "edge":
            x_step=self.h_step
            z_step=y_step=self.r_step
        else:
            raise NameError(self.direction)

        self.nx=int(self.lx/x_step)
        self.ny=int(self.ly/y_step)
        self.nz=int(self.lz/z_step)

        xArray = np.arange(-1.1*self.lx+0.5*x_step,1.1*self.lx-0.49*x_step,x_step)
        yArray = np.arange(-1.1*self.ly+0.5*y_step,1.1*self.ly-0.49*y_step,y_step)
        zArray = np.arange(-1.1*self.lz+0.5*z_step,1.1*self.lz-0.49*z_step,z_step)

        self.xArray = xArray
        self.yArray = yArray
        self.zArray = zArray

        Y,X,Z=np.meshgrid(yArray,xArray,zArray) #Feature of numpy.meshgrid
        if self.direction in ("top","bottom"):
            r_squared=X**2+Z**2
            self.projGrid=self.getCarrierDensity(Y,r_squared)*x_step*y_step*z_step*1e-18
        elif self.direction == "edge":
            r_squared=Y**2+Z**2
            self.projGrid=self.getCarrierDensity(X,r_squared)*x_step*y_step*z_step*1e-18

    def getCarrierDensity(self,h,r2):
        def getIntensity(h,r2):
            widthSquared= getWidthSquared(h)
            if self.tech=="SPA":
                intensity = ((2*self.power)/(np.pi*widthSquared*1e-12))*np.exp((-2*r2/(widthSquared)))*np.exp(-self.alpha*h*1e-6)
                return intensity
            elif self.tech=="TPA":
                k=(self.power**2)*8*np.log(2)/(self.tau*(np.pi**2.5)*(np.log(4))**0.5)
                intensity_squared = k*np.exp(-4*r2/widthSquared)/((widthSquared**2)*1e-24)
                return intensity_squared

        def getWidthSquared(h):#return um^2
            return (self.widthBeamWaist**2)*(1+(h/self.l_Rayleigh)**2)

        I = getIntensity(h,r2)
        if self.tech=="SPA":
            e0 = 1.60217733e-19
            return self.alpha*I/(3.6*e0)
        elif self.tech=="TPA":
            h_Planck = 6.626*1e-34
            speedofLight = 2.998*1e8
            return self.beta_2*self.wavelength*1e-6*I/(2*h_Planck*speedofLight)

    def getTrackProfile2D(self,x_rel,y_rel,z_rel):
        """
        Description:
            Transfer Carrier Distribution from Laser Coordinate System 
            to 2d Detector Coordinate System
        Parameters:
        ---------
        x_rel,y_rel,z_rel:
            the Normalized Coordinate for Laser Focus 
            in Detector Coordinate System
        @Modify:
        ---------
            2021/09/13
        """
        self.track_position = []
        self.ionized_pairs = []
        self.ionized_total_pairs = 0
        for i in range(len(self.projGrid)):
            x_div_point = self.xArray[i]+x_rel*self.lx
            if not 0 <= x_div_point <= self.lx:
                continue
            for j in range(len (self.projGrid[0])):
                y_div_point = self.yArray[j]+y_rel*self.ly
                if not 0 <= y_div_point <= self.ly:
                    continue
                z_sum = 0
                for k in range(len (self.projGrid[0][0])):
                    z_div_point = self.zArray[k]+z_rel*self.lz
                    if 0 <= z_div_point <= self.lz:
                        z_sum += self.projGrid[i][j][k]
                self.track_position.append([x_div_point,y_div_point])
                self.ionized_pairs.append(z_sum)
                self.ionized_total_pairs += z_sum
        print(self.ionized_total_pairs)
        self.draw_nocarrier(x_rel,y_rel,z_rel)

    def draw_nocarrier(self,x_rel,y_rel,z_rel):
        c1 = ROOT.TCanvas("c1","canvas2",200,10,1000,1000)
        h = ROOT.TH2D("h","pairs of carrier generation",self.nx,0,self.lx,self.ny,0,self.ly)
        for i in range(len(self.track_position)):
            h.Fill(self.track_position[i][0], self.track_position[i][1],self.ionized_pairs[i])
        h.Draw()
        c1.SaveAs("nocarrier_"+str(round(x_rel,5))+"_"+str(round(y_rel,5))+"_"+str(round(z_rel,5))+".pdf")  
