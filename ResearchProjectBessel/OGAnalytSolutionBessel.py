# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:07:34 2019
@author: ognyan
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, jn_zeros
import matplotlib.colors as colors
from matplotlib import cm

ax1 = plt.subplot(221)
Data1 = np.recfromcsv("ResearchProjectResultBessel1000.csv",names = "a")
Data1 = np.array(Data1)
Data2 = np.reshape(Data1,(-1,1000))
ImageData = Data2.astype("float64")
#ax = plt.subplot(211)
a = plt.imshow(ImageData,cmap = "plasma",origin = "lower")
plt.title("Jacobi Method")
plt.colorbar(a)
f1 =plt.contour(ImageData,colors='k')

#DEFINE INITIAL ARRAY DIMENSIONS:
ROWS, COLUMNS = 1000, 1000
SIZE = ROWS*COLUMNS

#DEFINE THE BOUNDS OF THE PHYSICAL SYSTEM:
w, h = 100, 100

#define step size values of the system
dr = w / (ROWS-1)
dz = h / (COLUMNS-1)

#INITIALIZE ARRAY
Array = np.zeros([ROWS,COLUMNS])
Array = Array.astype("float64")

def LaplaceSolution(Rows,width,height):
	"""Calculates the value for each point in the grid based on the analytical solution"""
	#extract the values from the grid where z = z_max and add in the C coefficient from the analytical solution
	r = np.linspace(0,width,Rows)
	z = np.linspace(0,height,Rows)
	z_max = z.max()
	r,z = np.meshgrid(r,z)
	#specify root and order, where order will be generalized
	root = jn_zeros(0,1)
	k = root / (w)
	J_0 = jv(0, k*r)
	Coeff = (np.exp((k*z_max))-np.exp((-k*z_max)))**(-1)
	soln = Coeff*J_0*(np.exp((k*z))-np.exp((-k*z)))
	return soln

V = 10
n = 1000
Height = 100*10**(-9)
width = 100*10**(-9)
Result = LaplaceSolution(n,width,Height)


x = Result[0]
AnalyticalData = Result

#plt.plot(x,y,"red")




### Visualisation
Data1 = np.recfromcsv("ResearchProjectResultBessel1000.csv",names = "a")
Data1 = np.array(Data1)
Data2 = np.reshape(Data1,(-1,1000))
ImageData = Data2.astype("float64")
#a = plt.imshow(ImageData,cmap = "plasma",origin = "lower")




plt.subplot(221)
p1 = plt.imshow(ImageData,origin = "lower")
cb1 = plt.colorbar(p1,shrink = 0.9)
cb1.set_label( "Voltage (V)")
plt.xlabel("Radius (m)")
plt.ylabel("Height (m)")


plt.title("Simulation")
plt.subplot(222)
p2 = plt.imshow(AnalyticalData,origin = "lower")
cb2 = plt.colorbar(p2,shrink = 0.9)
cb2.set_label("Voltage (V)")
plt.xlabel("Radius (m)")
plt.ylabel("Height (m)")
plt.title("Analytical Solution")



Convergence = abs(ImageData - AnalyticalData)/abs(AnalyticalData)

Convergence = Convergence[1:-1,1:-1]
plt.subplot(223)
p3 = plt.imshow(Convergence,origin = "lower")
plt.title("Convergence")
plt.xlabel("Radius (m)")
plt.ylabel("Height (m)")
cb3 = plt.colorbar(p3,shrink = 0.9)
cb3.set_label("Percentage error")
Convergence = 100*np.sum(Convergence)/(1000**2)
print("Percentage Error",Convergence)




NList = np.array([100,200,300,400,500,600,700,800,900,1000])
AccuracyList = []
for i in NList:
    Data1 = np.recfromcsv("ResearchProjectResultBessel"+str(i)+".csv",names = "a")
    Data1 = np.array(Data1)
    Data2 = np.reshape(Data1,(-1,i))
    ROWS, COLUMNS = i,i
    r = np.linspace(0,Height,ROWS)
    z = np.linspace(0,Height,COLUMNS)
    ImageData = Data2.astype("float64")
    Solution_final =  LaplaceSolution(i,width,Height)
    Convergence = abs(ImageData - Solution_final)/abs(Solution_final)
    Convergence = Convergence[1:-1,1:-1]
    Convergence = 100*np.sum(Convergence)/(i*i)
    AccuracyList.append(Convergence)
    print("Percentage Error",Convergence)



plt.subplot(224)
plt.semilogy(NList,AccuracyList, "x")
plt.title("Comparative Error")
plt.xlabel("Number of points along one axis, J")
plt.ylabel("Percentage Error")
plt.subplots_adjust(wspace = 0.7,hspace = 0.4)
plt.savefig("BesselGeneral.jpg",quality = 95)
plt.close()



ax = plt.subplot(111)
ax.plot(x,AnalyticalData[50,:],ImageData[50:],linewidth = 1)
ax.legend(("Analytical Solution","Simulation result"))
plt.xlabel("Number of points along one axis, J")
plt.ylabel("Percentage Error of total matrix")
plt.savefig("BesselPlot")
plt.close()