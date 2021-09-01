#!/usr/bin/env python3
import os
import sys
import time
import subprocess
class Input_parameters:
    def __init__(self,args):
        self.events_each_run = int(args[1].split("_")[-1]) # events/run
        self.events_total = int(args[2].split("_")[-1]) # events/total
        self.instance_in = int(args[3].split("_")[-1]) #singularity instance start number
        self.change_p = args[4] #parameters
        self.output_path=args[5].split("=")[-1]      #output name
    @property
    def instance_number(self): #singularity instance number
        return int(self.events_total/self.events_each_run)

def main():
    args = sys.argv[1:]
    input=Input_parameters(args)
    if "3D" in args[0]:
        sub_job(input,model="3D",run_code="./run 0.2.3 3D 3D_scan")
    elif "2D" in args[0]:
        sub_job(input,model="2D",run_code="./run 0.1.3 3D 2D_scan")
    else:
        print("the scan model is wrong")

def sub_job(input,model,run_code):

    for i in range(input.instance_number):
        intance_n  = "instance"+str(input.instance_in+i)
        runcmd("singularity instance start raser.simg "+intance_n)

    for i in range(input.instance_number):
        e_number=input.events_each_run*(i+1)
        runcmd('singularity exec instance://instance%s %s %s %s %s %s &' 
        %(input.instance_in+i,run_code,e_number,input.events_each_run,input.change_p,input.output_path))
        time.sleep(1)

def runcmd(command):
    ret = subprocess.run([command],shell=True)
if __name__ == '__main__':
    main()