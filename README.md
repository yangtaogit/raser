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

> $ python example.py det_model=planar3D parfile=setting.json
 
 - example.py 

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
 - setting.json 
  
```
[
{
"steplength" : "1",
"name" : "planar3D",

"lx" : "50",
"ly" : "50",
"lz" : "10",
"doping" : "10", 
"voltage" : "-100",
"temp" : "300.0",

"mesh" : "32",
"xyscale" : "1",

"maxstep" : "10",
"g4_vis" : "0",
"par_inx" : "25",
"par_iny" : "25",
"par_inz" : "17000",
"par_outx" : "25",
"par_outy" : "25",
"par_outz" : "0",

"t_rise" : "1",
"t_fall" : "1",
"trans_imp" : "1",
"CDet" : "1",
"BBW" : "1",
"BBGain" : "10",
"BB_imp" : "1",
"OscBW" : "1"
}
]

```


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
  
