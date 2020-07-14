# -*- coding: utf-8 -*-
"""
Created on Mon Mar 18 15:47:09 2019

@author: ognyan
"""
import numpy as np
import matplotlib.pyplot as plt
import time
from numba import jit, float64, int32, int64, prange
from matplotlib import cm
import matplotlib.colors as colors
###############################################################################
# Define the System
###############################################################################

#ARRAY DIMENSIONS: N ROWS, M COLUMNS
Dim = 1000
ROWS, COLS = Dim, Dim
size = ROWS * COLS
#DEFINE PHYSICAL SYSTEM: w - WIDTH; h - HEIGHT; dx & dy POINT SPACING
boundaryindex = (ROWS-1)//2
w, h = 10000,10000 # nm
dr = w / (COLS - 1)
dz = h / (ROWS - 1)
r = np.linspace(0, w, ROWS)
z = np.linspace(0, h, COLS)
eps1 = 1.0
eps2 = 4.1
lw = 2 # LINE WIDTH

###############################################################################
# Work Out Analytical Soln
###############################################################################

r,z = np.meshgrid(r,z)

zmax = z.max()
z0 = 0 # only solved analytically for bottom potential = 0
zb = z[boundaryindex]
zb = h-boundaryindex*dz
z2 = z[z <= zb]
z1 = z[z > zb]

e = -1.602e-19
gamma = 0.539
eps0 = 8.85418782e-12*(1e-9) # converted to units of F/nm
k = -e/(eps0*np.pi*gamma**2) # gamma left in units of eV nm

eps1 = 1
eps2 = 4.1

VB = 0 # BOTTOM POTENTIAL
VT = 1 # TOP POTENTIAL

# QUADRATIC COEFFICIENTS
A = (k*(zb-zmax)**2)
B = ((2*k*VT + eps2/zb)*(zb-zmax)-eps1)
C = (k*VT**2 + eps2*(VT-VB)/zb)

# CONSTANTS IN ANALYTICAL SOLN
a = (-B - np.sqrt((B**2-4*A*C)))/(2*A)
b = VT - a*zmax
c = (a*zb-a*zmax+VT-VB)/zb
d  = VB

# ANALYTICAL SOLNS
v1 = a*z1+b
v2 = c*z2+d


#plt.plot(z,v,'-', label='Analytical Solution', linewidth = lw+2)
#plt.xlabel('z (nm)')
#plt.ylabel('Potential (V)')



"""COMPARISON BETWEEN THIS ANALYTICAL SOLUTION AND THE NUMERICAL SOLUTION"""

TOTAL_POT = np.append(v2,v1)
Z = np.append(z2,z1)

"""Add as a separate plot -> this is just showing the linearity of the potential."""
#Plotting potential against Z
"""
plt.figure(1)
Data1 = np.recfromcsv("ChebyshevSiO2_simple.csv",names = "a")
Data1 = np.array(Data1)
Data2 = np.reshape(Data1,(-1,COLS))
plt.plot(Z,TOTAL_POT, color = "green",label = "Graphene")
plt.title('Potential magnitude against height \n Comparison between scenarios with and without a graphene layer')
plt.xlabel('Height (nm)')
plt.ylabel('Potential (V)')
boundaryforplot = (h-1)//2
Middle = boundaryforplot
line = plt.axvline(x=Middle, label = "Boundary")
plt.legend()
"""

"""Visualization"""
#Plotting numerical solution -> colormap/colorbar
plt.figure(2)
ax1 = plt.subplot(221)
Data1 = np.recfromcsv("ResearchProjectResultGraphene_1000.csv",names = "a")
Data1 = np.array(Data1)
Data2 = np.reshape(Data1,(-1,COLS))
ImageData = Data2.astype("float64")
a = plt.imshow(ImageData,extent=[-w/2,w/2,0,h],cmap=cm.inferno,origin="lower",norm=colors.PowerNorm(gamma=1./2.5))
plt.title("Numerical Solution -> Gauss-Seidel+Chebyshev")
plt.xlabel('Radius (nm)')
plt.ylabel('Height (nm)')
plt.colorbar(a)

#Plotting analytical solution -> colormap/colorbar
ax = plt.subplot(222)
RESHAPED = np.reshape(TOTAL_POT,(COLS,-1))
RESHAPED = RESHAPED.astype("float64")
RESHAPED = RESHAPED[:,::-1]
cp = plt.imshow(RESHAPED,extent=[-w/2,w/2,h,0],cmap=cm.inferno,norm=colors.PowerNorm(gamma=1./2.5))
ax.set_xlabel('X LABEL')    
plt.gca().invert_yaxis()
plt.title("Analytical solution")
plt.xlabel('Radius (nm)')
plt.ylabel('Height (nm)')
plt.colorbar(cp)
plt.show()

#Plotting the Error colourmap
plt.subplot(223)
Convergence = abs(ImageData[1:-1,1:-1] - RESHAPED[1:-1,1:-1])/abs(RESHAPED[1:-1,1:-1])
Convergence = Convergence[1:-2,1:-2]
cp = plt.imshow(Convergence, cmap = cm.inferno, norm=colors.PowerNorm(gamma=1./1.5), origin="lower")
plt.title("Error")
plt.colorbar(cp)
Convergence = 100*np.sum(Convergence)/(ROWS*COLS)
print("Percentage Error compared to analytical solution",Convergence)



#Numerical values of charge density (sigma) for both numerical and analytical and error associated with them
#Here ive done analytical solution again because the "z" needed in order to calculate sigma analytical is in a 
#different form to what i have it up in the code. 
#TO CONDENSE THE CODE: just give the z needed to calculate sigma a different name and use that
#up in the code so as
#not to be required to put the whole solution down here again
#however, not really needed as this is merely a comparison that is not very much computationally demanding 

z = np.linspace(0, h, ROWS)
zmax = z.max()
z0 = 0 # only solved analytically for bottom potential = 0
zb = z[boundaryindex]
z2 = z[z <= zb]
z1 = z[z > zb]
e = -1.60217662*10**(-19)
gamma = 0.539
eps0 = 8.85418782e-12*(1e-9) # converted to units of F/nm
k = -e/(eps0*np.pi*gamma**2) # gamma left in units of eV nm
eps1 = 1
eps2 = 4.1
VB = 0 # BOTTOM POTENTIAL
VT = 1 # TOP POTENTIAL
# QUADRATIC COEFFICIENTS
A = (k*(zb-zmax)**2)
B = ((2*k*VT + eps2/zb)*(zb-zmax)-eps1)
C = (k*VT**2 + eps2*(VT-VB)/zb)
# CONSTANTS IN ANALYTICAL SOLN
a = (-B - np.sqrt((B**2-4*A*C)))/(2*A)
b = VT - a*zmax
c = (a*zb-a*zmax+VT-VB)/zb
d  = VB
# ANALYTICAL SOLNS
v1 = a*z1+b
v2 = c*z2+d
z = np.append(z2,z1)



# SIGMA: CHARGE DENSITY, UNITS:
phi = ImageData[::1,1] # POTENTIAL,arbitrary r index chosen (as parallel plate), all z indecies
v = np.append(v2,v1)

sigmaanalytical = v[boundaryindex]**2*e/(np.pi*gamma**2)
sigmanumerical = phi[boundaryindex]**2*e/(np.pi*gamma**2)
print('\nCharge Density from analytical solution:',sigmaanalytical,'[C/nm^2]','\n\nCharge Density from numerical solution:',sigmanumerical,'[C/nm^2]')
diffsig = np.abs(sigmaanalytical-sigmanumerical)
print('\nDifference in charge densities from both methods',diffsig,'[C/nm^2]')
percenterrorsig = np.abs(diffsig/sigmaanalytical * 100)
print('\n% Error in charge densities from both methods', percenterrorsig)


# N: NUMBER DENSITY OF CHARGE CARRIERS
numberdensity_analytical = sigmaanalytical/e
numberdensity_numerical = sigmanumerical/e

print('\nNumber density of charge carriers from analytical solution:',numberdensity_analytical,'[1/nm^2]','\n\nNumber density of charge carriers from numerical solution:',numberdensity_numerical,'[1/nm^2]')
diffnum = np.abs(numberdensity_analytical-numberdensity_numerical)
print('\nDifference in number densities from both methods:',diffnum,'[1/nm^2]')
error_in_numbers_percent = np.abs(diffnum/numberdensity_analytical * 100)
print('\n% Error in number densities from both methods:',error_in_numbers_percent)