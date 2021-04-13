import numpy as np
import test2
import sys
import pandas as pd
import os
'''
Xi = np.array([62.34103, 62.98372, 63.62641, 64.2691, 64.91179, 65.55448, 66.19717, 66.83986, 67.48256, 68.12525])
Yi = np.array([-2088.29101, -2088.32331, -2088.34236, -2088.34881, -2088.34324, -2088.32624, -2088.29834, -2088.26006,
      -2088.2119, -2088.15431])
print(test2.calculate(Xi,Yi))
'''
def main(filename,storefile):
      xi = []
      yi = []
      datax = []
      datay = []
      print(filename)
      print(storefile)
      with open(filename,"r") as f:
            line = f.readline()
            num = 0
            while line:
                  xi.append(float(line.split()[0]))
                  yi.append(float(line.split()[1]))
                  num += 1
                  if(num%10==0):
                        datax .append( np.array(xi) )
                        datay .append( np.array(yi) )
                        xi = []
                        yi = []
                  line = f.readline()
      df = pd.DataFrame(columns=['V0','B0','E0','B1'])
      for i in range(len(datax)):
            df.loc[df.shape[0]+1] = test2.calculate(datax[i],datay[i])
      df.to_csv(storefile)
      print(storefile)

if __name__ == '__main__':
      filename = sys.argv[1]
      storefile = sys.argv[2]
      main(filename, storefile)