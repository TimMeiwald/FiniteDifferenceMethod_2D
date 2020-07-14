import numpy as np
import matplotlib.pyplot as plt


z = np.linspace(0,100,101)
zmax = 100
zb = z[50]


e = -1.602e-19
gamma = 0.539
eps0 = 8.85418782e-12*(1e-9) # converted to units of F/nm
k = -e/(eps0*np.pi*gamma**2) # gamma left in units of eV nm
print(k)
#sigma = 0.2**2/(np.pi*gamma**2)
#print(sigma)

eps1 = 1
eps2 = 4.1

VB = 0
VT = 10

A = (k*(zb-zmax)**2)
B = ((2*k*VT + eps2/zb)*(zb-zmax)-eps1)
C = (k*VT**2 + eps2*(VT-VB)/zb)

print(A,B,C)

a = (-B - np.sqrt((B**2-4*A*C)))/(2*A)
b = VT - a*zmax
c = (a*zb-a*zmax+VT-VB)/zb
d  = VB

print(a,b,c,d)

plt.plot(z[50:],a*z[50:]+b)
#plt.plot(z[50:],aa*z[50:]+bb)
plt.plot(z[0:51],c*z[0:51]+d)
#plt.plot(z[0:51],cc*z[0:51]+d)
#plt.xlim(0,zmax)
#plt.ylim(0,VT)
plt.show()
