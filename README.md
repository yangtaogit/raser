# RASER 
RAdiation SEmiconductoR

# Function
1.Simulate the time resolution of 2D and 3D SiC or SiC detector
2.Simulate the TCT and TPA scan

# Dependent packages and environments
1.python fenics
2.pybind-geant4
3.python-tk, python-ipython, tk-dev
4.build-essential libboost-all-dev qtcreator qt5-default python3-pip
5.libgl1-mesa-dev libglu1-mesa-dev libxt-dev libxmu-dev libxi-dev zlib1g-dev
libgl2ps-dev libexpat1-dev libxerces-c-dev

# Reference run setting
1.run program like:
./pyraser/raser.py det_model=planar3D parfile=setting.json
2. setting.json: define all paramters you want input
3. raser.py examples:
from setting import Setting
from geometry import R3dDetector
from pyfenics import FenicsCal
from g4particles import Particles
from calcurrent import CalCurrent

dset = Setting(args)
my_d = R3dDetector(dset)
my_f = FenicsCal(my_d, dset.fenics)
my_g4p = Particles(my_d, my_f, dset)
my_current = CalCurrent(my_d, my_f, my_g4p, dset)
