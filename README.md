RASER
======

RAdiation SEmiconductoR

https://pypi.org/project/raser

Function
======

- Simulate the time resolution of 2D and 3D SiC or SiC detector
- Simulate the TCT and TPA scan

Dependences 
======

- python fenics
- pybind-geant4
- python-tk, python-ipython, tk-dev
- build-essential libboost-all-dev qtcreator qt5-default python3-pip
- libgl1-mesa-dev libglu1-mesa-dev libxt-dev libxmu-dev libxi-dev zlib1g-dev
  libgl2ps-dev libexpat1-dev libxerces-c-dev


Example 
======

- An example python code to use raser:
``` 
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
```

> $ python example.py det_model=planar3D parfile=setting.json


Contribution 
====== 
* Xin Shi, IHEP, @xshi
* Yuhang Tan, IHEP, @tanyh2018
* Tao Yang, IHEP, @yangtaogit
* Kai Liu, IHEP, @liukaihep
* Ryuta Kiuchi, IHEP, @rkiuchi
* Jianing Lin, Jilin U, @Zombastar
* Yu Zhao, Liaoning U, @zylyz18
* Ziyi Zhao, Hunan Normal U, @zhaoziyi1
* Jia Wang, Hunan Normal U, @wangjia0203 
* Chenxi Fu, Jilin U, @heheda2001123
  
