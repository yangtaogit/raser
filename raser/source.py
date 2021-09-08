import math
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
