import os
import sys
import time
def main():
    args = sys.argv[1:]
    model_sensor = args[0]
    output       = args[1]
    number       = int(args[2])
    t_number     = int(args[3])
    step_n		 = int(args[4])	
    change_p     = int(args[5])

    for i in range (number-step_n,number):
        os.system("python3 python/raser_3D.py %s %s %s %s %s %s &" %(model_sensor,output,i,t_number,step_n,change_p))
        time.sleep(1)
if __name__ == '__main__':
    main()