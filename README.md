RASER
======

**RA**diation **SE**miconducto**R** 


[![Latest Version][pypi-image]][pypi-url] 
  

Function
======

- SiC Timing Detector Simulation 

Dependences 
======

- [ROOT](https://root.cern.ch/building-root) 
- [Geant4](https://geant4.web.cern.ch)
- [FEniCS](https://fenicsproject.org)


Example 
======

- An example python code to use raser:

> $ python example.py det_model=planar3D parfile=setting.json
 
 - example.py 

```python 
import raser
dset = raser.Setting(args)
my_d = raser.R3dDetector(dset)
my_f = raser.FenicsCal(my_d, dset.fenics)
my_g4p = raser.Particles(my_d, my_f, dset)
my_current = raser.CalCurrent(my_d, my_f, my_g4p, dset)
ele_current = raser.Amplifier(my_d, dset.amplifer)
```
 - setting.json 
  
```json
[{
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
}]

```


Contribution 
====== 
* Xin Shi, IHEP, @xshi
* Yuhang Tan, IHEP, @tanyh2018
* Tao Yang, IHEP, @yangtaogit
* Kai Liu, IHEP, @liukaihep
* Ryuta Kiuchi, IHEP, @rkiuchi
* Jianing Lin, Jilin U, @zombastar
* Yu Zhao, Liaoning U, @zylyz18
* Ziyi Zhao, Hunan Normal U, @zhaoziyi1
* Jia Wang, Hunan Normal U, @wangjia0203 
* Chenxi Fu, Jilin U, @heheda2001123
  


[pypi-image]: https://img.shields.io/pypi/v/raser.svg
[pypi-url]: https://pypi.org/project/raser

