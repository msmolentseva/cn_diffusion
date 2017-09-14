#diffusion equation du/dt=c*d2u/dx2 in 1D
#Crank-Nicolson discretization scheme in time, centered discretization scheme in space
#periodic boundary conditions
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sl
from mpl_toolkits.mplot3d.axes3d import Axes3D
plt.ion()

#define initial function
def u(x,y):
  x0=10.5
  y0=10.5
  sigmax=10.0
  sigmay=10.0
  if (np.abs(x-x0)<=sigmax):
    res=0.5*(1.0+np.cos(np.pi*(x-x0)/sigmax))
  #  res=0.5*np.exp(1.0-4*np.log(2)*((x-x0)**2/(2.0*sigmax**2)))
  else:
    res=0.0
  return res

dx=0.4
dy=0.4
L=np.arange(0.0,50.0,dx)
H=np.arange(0.0,50.0,dy)
dt=0.5
c=0.1
n=len(L)
m=len(H)
du=np.zeros((n,m))
du1=np.zeros((n,m))

print 2.0*c*dt/np.min(dx**2,dy**2)

mu=0.5*dt*c/dx**2
la=0.5*dt*c/dy**2

for i in range(0,n):
  for j in range(0,m):
    du[i,j]=u(L[i],H[j])
#print du
#fig = plt.figure(figsize=(14,6))
#ax = fig.add_subplot(1, 1, 1, projection='3d')
#X,Y=np.meshgrid(L,H)
#p=ax.plot_surface(X,Y, du, rstride=4, cstride=4, linewidth=0)
#ax.set_zlim(0, 1.0)
#cb = fig.colorbar(p, shrink=0.5)
#levels=np.linspace(np.min(du),np.max(du),40)
#plt.contourf(L, H, du,levels=levels)
#plt.show()
c1=0.5
for t in np.arange(0.0,10.1,dt):
  #predictor du1
  #mu=0.5*dt*c/dx**2
  for j in range(1,m-1):
    du1[0,j]=mu*(du[1,j]) + (1.0-2.0*mu-2.0*la)*du[0,j] + la*(du[0,j+1] + du[0,j-1])
    du1[n-1,j]=mu*(du[n-2,j]) + (1.0-2.0*mu-2.0*la)*du[n-1,j] + la*(du[n-1,j+1] + du[n-1,j-1])
  for i in range(1,n-1):
    du1[i,0]=mu*(du[i+1,0] + du[i-1,0]) + (1.0-2.0*mu-2.0*la)*du[i,0] + la*du[i,1]
    du1[i,m-1]=mu*(du[i+1,m-1] + du[i-1,m-1]) + (1.0-2.0*mu-2.0*la)*du[i,m-1] + la*du[i,m-2]
#    du1[i,0]=mu*(du[i+1,0] + du[i-1,0]) + (1.0-2.0*mu-2.0*la)*du[i,0] + la*(du[i,1] + du[i,n-1])
#    du1[i,m-1]=mu*(du[i+1,m-1] + du[i-1,m-1]) + (1.0-2.0*mu-2.0*la)*du[i,m-1] + la*(du[i,0] + du[i,m-2])
    for j in range(1,m-1):
      du1[i]=mu*(du[i+1,j] + du[i-1,j]) + (1.0-2.0*mu-2.0*la)*du[i,j] + la*(du[i,j+1] + du[i,j-1])
  du1[0,0]=mu*(du[1,0]) + (1.0-2.0*mu-2.0*la)*du[0,0] + la*(du[0,1])
  du1[n-1,m-1]=mu*(du[n-2,m-1]) + (1.0-2.0*mu-2.0*la)*du[n-1,m-1] + la*(du[n-1,m-2])
  du=du1.copy()

  #plt.clf()
#plt.ylim((0,5.0e-9))
fig = plt.figure(figsize=(14,6))
ax = fig.add_subplot(1, 1, 1, projection='3d')
X,Y=np.meshgrid(L,H)
p=ax.plot_surface(X,Y, du, rstride=4, cstride=4, linewidth=0)
ax.set_zlim(0, 1.0)
 # plt.pause(0.1)
#plt.show()
raw_input()
