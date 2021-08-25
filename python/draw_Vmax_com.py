#!/usr/bin/env python3

from operator import mod
import os,sys,re,traceback
from datetime import datetime
from string import Template
import sys

def draw_vmax_pin_compre():
    A = [60,80,100,120,140,160,180,200]
    for i in A:
        template_file = open(r'./python/draw/plt_Vmax_pin.c','r')
        tmpl = Template(template_file.read())
        lines = []
        lines.append(tmpl.substitute(VBIAS=i))
        python_file_name = './python/draw/temp_plt_Vmax_pin.c'
        python_file = open(python_file_name,'w')
        python_file.writelines(lines)
        python_file.close()
        os.system('root -q -b -l ./python/draw/temp_plt_Vmax_pin.c')

def draw_vmax_pin_vbias():
    python_file_name = './python/draw/draw_vmax_difV.c'
    os.system('root -q -b -l ./python/draw/draw_vmax_difV.c')

def main():
    args = sys.argv[1:]
    model = args[0]

    if model in ["pin_compre"]:
        draw_vmax_pin_compre()

    elif model in ["pin_vbias"]:
        draw_vmax_pin_vbias()

    else:
        raise NameError(model)

if __name__ == '__main__':
    main() 