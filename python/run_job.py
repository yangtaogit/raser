#!/usr/bin/env python3
import os
import sys
import time
def main():
    args = sys.argv[1:]
    n_number = int(args[0])
    t_number = int(args[1])
    change_p = float(args[2])
    print(change_p)
    n_in = int(args[3])
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