import math
import ROOT

""" Define Track """

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

    def nocarrier(self,r):
        z_o = self.z_o
        rlen = self.rlen
        zlen = self.zlen
        r_min = -r   #um
        r_max = -r + rlen
        r_nBins = self.i_r
        z_min = 0.    #um
        z_max = zlen  #um
        z_nBins = self.i_z

        tau = self.tau
        alfa = self.alfa
        power = self.power
        wavelength = self.wavelength
        widthBeamWaist = self.widthBeamWaist
        refractionIndex = self.refractionIndex

        my_pro = ProjGrid(r_min, r_max, r_nBins, z_max, z_min, z_nBins, z_o)
        rGrid, zGrid = np.meshgrid(my_pro.rArray,my_pro.zArray)
        carriergeneration = SPAGeneration(tau,alfa,power,wavelength,widthBeamWaist,refractionIndex)
        CGrid = carriergeneration.getCarrierDensity(rGrid, zGrid, z_o*1e-6)
        xGrid = rGrid.copy()
        projGrid = CGrid.copy()

        for i_z in range(z_nBins):
            for i_r in range(r_nBins):
                x_value = xGrid[i_z, i_r]
                r_valueArray = np.sqrt(x_value*x_value+ my_pro.yArray*my_pro.yArray)
                # Get z value
                z_value = zGrid[i_z, i_r]
                z_valueArray = np.ones_like(r_valueArray)*z_value
                carr_den_projArray = carriergeneration.getCarrierDensity(r_valueArray, z_valueArray, z_o*1e-6)
                # Project sum and take into account integral step
                projGrid[i_z, i_r] = (carr_den_projArray.sum()*my_pro.y_step)/2
        self.ionized_pairs = projGrid*my_pro.x_step*my_pro.z_step
        print(self.ionized_pairs)

class SPAGeneration():
    def __init__(self,tau,alfa,power,wavelength,widthBeamWaist,refractionIndex):
        ### for Si 1064 IR ###
        self.tau = tau      #ps
        self.alfa = alfa         #Si
        self.power = power      #J/s
        self.wavelength = wavelength  #m
        self.widthBeamWaist = widthBeamWaist  #m
        self.refractionIndex = refractionIndex

    def getWidthSquared(self,z):
        return (self.widthBeamWaist**2)*(1+((self.wavelength*z)/(np.pi* (self.widthBeamWaist**2)*self.refractionIndex))**2)
    def getWidth(self,z):
        return np.sqrt(self.getWidthSquared(z))
    def getIntensity(self,r,z,z_o = 0 ):
        widthSquared= self.getWidthSquared(z-z_o)
        intensity = ((2*self.power)/(np.pi*widthSquared))*np.exp((-2*(r**2)/(widthSquared)))*np.exp(-self.alfa*z)
        return intensity
    def getCarrierDensity(self,r,z,z_o = 0 ):
        I = self.getIntensity(r,z,z_o)
        return (self.alfa*I)/(3.6*1.60217657e-19)

class ProjGrid():
    def __init__(self,r_min,r_max,r_nBins,z_max,z_min,z_nBins,z_o):
        self.r_min = r_min
        self.r_max = r_max
        self.r_nBins = r_nBins
        self.z_min = z_min
        self.z_max = z_max
        self.z_nBins = z_nBins
        self.z_o = z_o
        self.rArray = np.linspace(r_min*1e-6,r_max*1e-6,r_nBins)
        self.zArray = np.linspace(z_min*1e-6,z_max*1e-6,z_nBins)
        self.r_step = abs(r_max*1e-6-r_min*1e-6)/r_nBins
        self.z_step = abs(z_max*1e-6-z_min*1e-6)/z_nBins
        self.xArray = self.rArray
        self.yArray = self.rArray
        self.x_step = self.r_step
        self.y_step = self.r_step
