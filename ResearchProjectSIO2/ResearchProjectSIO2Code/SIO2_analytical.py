# -*- coding: utf-8 -*-
"""
Created on Tue Jan 29 15:07:34 2019
@author: ognyan
note: CLEAN UP THIS SHIT
"""
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import jv, jn_zeros
from matplotlib import cm

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

def NoCharge_Dielectric(TOPSIDEBESSEL,Z,VOLT_CONST,Rows): 
	"""Calculates the potential inside the cylinder with a Silicon Dioxide layer"""
	#Defining sides and boundary
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
	A_CONST = VOLT_CONST/(Z_MAX+Z_BOUND*(EPSILON_VACUUM/EPSILON_SILICON-1))
	#Defining analytical solutions depending on regions
	POTENTIAL_IN_SILICON = EPSILON_VACUUM/EPSILON_SILICON*A_CONST*SILICON_REGION
	POTENTIAL_IN_VACUUM = A_CONST*VACUUM_REGION + Z_BOUND*(EPSILON_VACUUM/EPSILON_SILICON-1)*A_CONST
	#Calculating total potential in system
	TOTAL_POTENTIAL = np.append(POTENTIAL_IN_SILICON,POTENTIAL_IN_VACUUM)
	Z = np.append(SILICON_REGION,VACUUM_REGION)
	
	return TOTAL_POTENTIAL, Z


TOTAL_POT = NoCharge_Dielectric(TOPSIDEBESSEL,Z,10,ROWS)[0]
Z = NoCharge_Dielectric(TOPSIDEBESSEL,Z,10,ROWS)[1]
print(TOTAL_POT)
#Plotting analytical solution
fig = plt.figure(figsize = (20,20))
ax = plt.subplot(524)
plt.plot(Z[::-1],TOTAL_POT)


#Plotting analytical solution -> colormap/colorbar
ax = plt.subplot(521)
RESHAPED = np.reshape(TOTAL_POT,(-1,ROWS))
RESHAPED = RESHAPED.astype("float64")
RESHAPED = RESHAPED[:,::-1]
plt.title("Analytical")
cp = plt.imshow(RESHAPED,extent=[-WIDTH/2,WIDTH/2,HEIGHT,0],cmap=cm.plasma)

ax.xaxis.tick_top()
ax.set_xlabel('X LABEL')    
ax.xaxis.set_label_position('bottom') 
plt.gca().invert_yaxis()
plt.xlabel('r')
plt.ylabel('z')
plt.colorbar(cp)
plt.show()


ax1 = plt.subplot(522)
Data1 = np.recfromcsv("ResearchProjectResultSIO2_100.csv",names = "a")
Data1 = np.array(Data1)
Data2 = np.reshape(Data1,(-1,100))
ImageData = Data2.astype("float64")
a = plt.imshow(ImageData,cmap = "plasma",origin = "lower")
plt.title("Simulation")
plt.colorbar(a)
#f1 =plt.contour(ImageData,colors='k')



#CONVERGENCE/ERROR
plt.subplot(523)
Convergence = abs(ImageData[1:-1,1:-1] - RESHAPED[1:-1,1:-1])/abs(RESHAPED[1:-1,1:-1])
Convergence = Convergence[1:-2,1:-2]
cp = plt.imshow(Convergence)
plt.colorbar(cp)
Convergence = 100*np.sum(Convergence)/(ROWS*COLS)
print("Percentage Error",Convergence)


#Compute the simulation for every NXN from the array below
#Store that in separate csv files
#Use the code below to compute % error against dimension size of matrix

NList = np.array([100,200,300,400,500,600,700,800,900,1000])
AccuracyList = []
for i in NList:
    Data1 = np.recfromcsv("ResearchProjectResultSIO2_"+str(i)+".csv",names = "a")
    Data1 = np.array(Data1)
    Data2 = np.reshape(Data1,(-1,i))
    ROWS, COLS = i,i
    w, h = 100, 100
    TOP = np.linspace(-w/2,w/2,i)
    H_STEP_SIZE = HEIGHT / (i-1)
    z = np.linspace(0,HEIGHT-H_STEP_SIZE,i)

    ImageData = Data2.astype("float64")
    
    TOTAL_POT = NoCharge_Dielectric(TOP,z,10,i)[0]
    RESHAPED = np.reshape(TOTAL_POT,(-1,i))
    RESHAPED = RESHAPED.astype("float64")
    RESHAPED = RESHAPED[:,::-1]


    
    Convergence = abs(ImageData[1:-1,1:-1] - RESHAPED[1:-1,1:-1])/abs(RESHAPED[1:-1,1:-1])
    Convergence = Convergence[1:-2,1:-2]
    Convergence = 100*np.sum(Convergence)/(i*i)
    AccuracyList.append(Convergence)
    print("Percentage Error",Convergence)

plt.subplot(525)
plt.semilogy(NList,AccuracyList, "x")
plt.title("Error when comparing simulation to analytical solution")
plt.xlabel("N, for NxN matrix")
plt.ylabel("% Error")










"""3D Visualization"""
#3D VISUALISATION: UNCOMMENT BELOW TO PLOT IN 3D
#figure = plt.figure(3,figsize=(10,8))
#ax = figure.add_subplot(111,projection='3d')
#plt.subplots_adjust(bottom = 0, top = 1, left = 0,right = 1)

#r = np.linspace(-0.5*w,0.5*w,ROWS)
#angular = np.linspace(0,2*np.pi,50)
#R, angular = np.meshgrid(r, angular)
#ord = 0
#root = jn_zeros(ord,1)[-1]
#k = root / (0.5*w)
#Z = jv(ord,k*R)

#X,Y = R*np.cos(angular), R*np.sin(angular)
#ax.plot_surface(X,Y,Z, cmap = cm.plasma, alpha = 0.8)

#plt.show()