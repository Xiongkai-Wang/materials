###最小二乘法试验###
import numpy as np
from scipy.optimize import leastsq

###采样点(Xi,Yi)###


###需要拟合的函数func及误差error###
def func(p,x):
    V0,B1,E0,B0 = p
    v0v = np.power(V0 / x, 2 / 3)
    temp = np.power(v0v - 1, 3) * B1 + np.power(v0v - 1, 2) * (6 - 4 * v0v)
    y = E0 + (9.0 / 16) * V0 * B0 * temp
    return y

def error(p,x,y,s):
    print(s)
    return func(p,x)-y #x、y都是列表，故返回值也是个列表


def calculate(Xi,Yi,p0=None):
    if p0 == None:
        p0 = [min(Xi), 1, min(Yi), 5]
#TEST

#print( error(p0,Xi,Yi) )

###主函数从此开始###
    s="Test the number of iteration" #试验最小二乘法函数leastsq得调用几次error函数才能找到使得均方误差之和最小的a~c
    Para=leastsq(error,p0,args=(Xi,Yi,s)) #把error函数中除了p以外的参数打包到args中
    V0,B1,E0,B0=Para[0]
    print("V0=",V0,'\n',"B0=",B0,"E0=",E0,"B1=",B1)
    return [V0,B0,E0,B1];
###绘图，看拟合效果###
'''
import matplotlib.pyplot as plt

plt.figure(figsize=(8,6))
plt.scatter(Xi,Yi,color="red",label="Sample Point",linewidth=3) #画样本点


plt.legend()
plt.show()
'''