#!/usr/bin/env python3
import os
import sys
import time
def main():
    args = sys.argv[1:]
    args[0]
    if "3D" in args[0]:
        threeD_job(args)
    elif "2D" in args[0]:
        twoD_job(args)
    else:
        print("the scan model is wrong")

def threeD_job(args):
    n_number = int(args[1].split("_")[-1])
    t_number = int(args[2].split("_")[-1])
    change_p = args[3]
    n_in = int(args[4].split("_")[-1])
    n=int(t_number/n_number)
    
    for i in range(n):
        intance_n  = "instance"+str(n_in+i)
        os.system("singularity instance start raser.simg "+intance_n)

    for i in range(n):
        e_number=n_number*(i+1)
        print(e_number)
        os.system('singularity exec instance://instance%s ./run 0.2.3 %s %s %s &' %(i,e_number,n_number,change_p))
        time.sleep(1)

def twoD_job(args):
    n_number = int(args[1].split("_")[-1])
    t_number = int(args[2].split("_")[-1])
    change_p = args[3]
    n_in = int(args[4].split("_")[-1])
    print(n_number,t_number,change_p,n_in)

    n=int(t_number/n_number)
    for i in range(n):
        intance_n  = "instance"+str(n_in+i)
        os.system("singularity instance start raser.simg "+intance_n)

    for i in range(n):
        e_number=n_number*(i+1)
        print(e_number)
        os.system('singularity exec instance://instance%s ./run 0.1.3 %s %s %s &' %(i,e_number,n_number,change_p))
        time.sleep(1)
if __name__ == '__main__':
    main()