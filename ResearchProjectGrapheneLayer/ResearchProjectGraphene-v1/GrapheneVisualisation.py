# -*- coding: utf-8 -*-
"""
Created on Fri Feb 22 17:17:48 2019

@author: timme
"""
import numpy as np
import matplotlib.pyplot as plt
from matplotlib import cm

Dim = 100
f = plt.figure(figsize = (10,10))
ax1 = plt.subplot(211)
Data1 = np.recfromcsv("ResearchProjectResultGraphene_100.csv",names = "a")
Data1 = np.array(Data1)
Data2 = np.reshape(Data1,(-1,Dim))
ImageData = Data2.astype("float64")
a = plt.imshow(ImageData,cmap = "plasma",origin = "lower")
plt.title("Simulation")
plt.colorbar(a)






"""Parameters"""
ROWS, COLS = 100, 100
SIZE = ROWS*COLS
WIDTH, HEIGHT = 100, 100
TOPSIDEBESSEL = np.linspace(-WIDTH/2,WIDTH/2,ROWS)
R_STEP_SIZE = WIDTH / (ROWS-1)
H_STEP_SIZE = HEIGHT / (COLS-1)
Z = np.linspace(0,HEIGHT-H_STEP_SIZE,COLS)

"""Initialization"""
Array = np.zeros([ROWS,COLS], dtype = np.float64)

def NoCharge_Dielectric(TOPSIDEBESSEL,Z,VOLT_CONST): 
	"""Calculates the potential inside the cylinder with a Silicon Dioxide layer"""
	#Defining sides and boundary
	TOPSIDEBESSEL,Z = np.meshgrid(TOPSIDEBESSEL,Z)
	Z_MAX = Z.max()
	BOUNDARY = int((ROWS)/2) #this is where the boundary is
	#Z_BOUND = HEIGHT - BOUNDARY*H_STEP_SIZE
	Z_BOUND = BOUNDARY
	#Differentiating material regions
	SILICON_REGION = Z[Z <= Z_BOUND]
	VACUUM_REGION = Z[Z > Z_BOUND]
	#Defining permitivities
	EPSILON_VACUUM = 1
	EPSILON_SILICON = 4.1
	#Defining analytical solution constant
	A_CONST = VOLT_CONST/(Z_MAX+Z_BOUND*(EPSILON_VACUUM/EPSILON_SILICON-1))
	#Defining analytical solutions depending on regions
	POTENTIAL_IN_SILICON = EPSILON_VACUUM/EPSILON_SILICON*A_CONST*SILICON_REGION
	POTENTIAL_IN_VACUUM = A_CONST*VACUUM_REGION + Z_BOUND*(EPSILON_VACUUM/EPSILON_SILICON-1)*A_CONST
	#Calculating total potential in system
	TOTAL_POTENTIAL = np.append(POTENTIAL_IN_SILICON,POTENTIAL_IN_VACUUM)
	Z = np.append(SILICON_REGION,VACUUM_REGION)
	
	return TOTAL_POTENTIAL, Z


TOTAL_POT = NoCharge_Dielectric(TOPSIDEBESSEL,Z,10)[0]
Z = NoCharge_Dielectric(TOPSIDEBESSEL,Z,10)[1]

#Plotting analytical solution
plt.subplot(212)
y = ImageData[::-1,int(Dim/2)]
x = np.arange(Dim)
ax = plt.plot(x,y,"-.")
plt.plot(Z[::-1],TOTAL_POT)
plt.legend(("Graphene","SIO2"))
plt.xlabel("Z(Height)")
plt.ylabel("V(Potential)")

plt.close(f)


y = ImageData[:,int(Dim/2)]
x = np.arange(Dim)
ax2 = plt.plot(x,y,"-")
plt.plot(Z[:],TOTAL_POT)
plt.legend(("Graphene","SIO2"))
plt.xlabel("Z(Height)")
plt.ylabel("V(Potential)")
plt.axvline(50,ls = "-.",color = "r")