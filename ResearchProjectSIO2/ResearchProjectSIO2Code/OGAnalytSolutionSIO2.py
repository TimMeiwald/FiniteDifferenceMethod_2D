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

ax1 = plt.subplot(321)
Data1 = np.recfromcsv("ResearchProjectResultSIO2_100.csv",names = "a")
Data1 = np.array(Data1)
Data2 = np.reshape(Data1,(-1,100))
ImageData = Data2.astype("float64")
#ax = plt.subplot(211)
a = plt.imshow(ImageData,cmap = "plasma",origin = "lower")
plt.title("Jacobi Method")
plt.colorbar(a)
f1 =plt.contour(ImageData,colors='k')

#DEFINE INITIAL ARRAY DIMENSIONS:
ROWS, COLUMNS = 100, 100
SIZE = ROWS*COLUMNS

#DEFINE THE BOUNDS OF THE PHYSICAL SYSTEM:
w, h = 100, 100
r = np.linspace(-w/2,w/2,ROWS)
z = np.linspace(0,h,COLUMNS)
#define step size values of the system
dr = w / (ROWS)
dz = h / (COLUMNS)

#INITIALIZE ARRAY
Array = np.zeros([ROWS,COLUMNS])
Array = Array.astype("float64")


VOLT_CONST = 10


def NoCharge_Dielectric(VOLT_CONST,HEIGHT,WIDTH,Rows,COLS): 
	"""Calculates the potential inside the cylinder with a Silicon Dioxide layer"""
	#Defining sides and boundary
	H_STEP_SIZE = HEIGHT / (COLS-1)
	Z = np.linspace(0,HEIGHT-H_STEP_SIZE,COLS)
	TOPSIDEBESSEL = np.linspace(-WIDTH/2,WIDTH/2,Rows)
	TOPSIDEBESSEL,Z = np.meshgrid(TOPSIDEBESSEL,Z)
	Z_MAX = Z.max()
	BOUNDARY = int((Rows)/2) #this is where the boundary is
	#Z_BOUND = HEIGHT - BOUNDARY*H_STEP_SIZE
	Z_BOUND = BOUNDARY
	#Differentiating material regions
	SILICON_REGION = Z[Z <= Z_BOUND]
	VACUUM_REGION = Z[Z > Z_BOUND]
	#Defining permitivities
	EPSILON_VACUUM = 1
	EPSILON_SILICON = 4.1
	#Defining analytical solution constant
	A_CONST = VOLT_CONST/(Z_MAX+Z_BOUND*(EPSILON_VACUUM/EPSILON_SILICON))
	#Defining analytical solutions depending on regions
	POTENTIAL_IN_SILICON = EPSILON_VACUUM/EPSILON_SILICON*A_CONST*SILICON_REGION
	POTENTIAL_IN_VACUUM = A_CONST*VACUUM_REGION + Z_BOUND*(EPSILON_VACUUM/EPSILON_SILICON-1)*A_CONST
	#Calculating total potential in system
	TOTAL_POTENTIAL = np.append(POTENTIAL_IN_SILICON,POTENTIAL_IN_VACUUM)
	#Z = np.append(SILICON_REGION,VACUUM_REGION)
	TOTAL_POTENTIAL = np.reshape(TOTAL_POTENTIAL,(-1,Rows))
	TOTAL_POTENTIAL = TOTAL_POTENTIAL.astype("float64")
	TOTAL_POTENTIAL = TOTAL_POTENTIAL[:,::-1]
	return TOTAL_POTENTIAL

Solution_final = NoCharge_Dielectric(VOLT_CONST,h,w,ROWS,COLUMNS)

#2D VISUALISATION
ax = plt.subplot(322)
cp = plt.imshow(Solution_final,extent=[-w/2,w/2,h,0],cmap=cm.plasma)
ax.xaxis.tick_top()
ax.set_xlabel('X LABEL')    
ax.xaxis.set_label_position('top') 
plt.gca().invert_yaxis()
plt.xlabel('r')
plt.ylabel('z')
f2 =plt.contour(r, z, Solution_final,colors='k',origin = "lower")
plt.colorbar(cp)
plt.show()


plt.subplot(323)
Convergence = abs(ImageData - Solution_final)/abs(Solution_final)

Convergence = Convergence[1:-1,1:-1]
cp = plt.imshow(Convergence,origin = "lower")
plt.colorbar(cp)
Convergence = 100*np.sum(Convergence)/(ROWS*COLUMNS)
print("Percentage Error",Convergence)


plt.subplot(324)
x = np.arange(0,COLUMNS)
plt.plot(x,ImageData[:,50])


"""
NList = np.array([100,200,300,400,500,600,700,800,900,1000])
AccuracyList = []
for i in NList:
    Data1 = np.recfromcsv("ResearchProjectResultSIO2_"+str(i)+".csv",names = "a")
    Data1 = np.array(Data1)
    Data2 = np.reshape(Data1,(-1,i))
    ROWS, COLUMNS = i,i
    w, h = 100, 100

    ImageData = Data2.astype("float64")
    Solution_final =  NoCharge_Dielectric(VOLT_CONST,h,w,i,i)
    print(Solution_final)
    Convergence = abs(ImageData - Solution_final)/abs(Solution_final)
    Convergence = Convergence[1:-1,1:-1]
    Convergence = 100*np.sum(Convergence)/(i*i)
    AccuracyList.append(Convergence)
    print("Percentage Error",Convergence)




plt.subplot(325)
plt.semilogy(NList,AccuracyList, "x")
plt.title("Error when comparing simulation to analytical solution")
plt.xlabel("N, for NxN matrix")
plt.ylabel("% Error")
"""





