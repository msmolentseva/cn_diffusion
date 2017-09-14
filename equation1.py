#diffusion equation du/dt=c*d2u/dx2 in 1D
#Crank-Nicolson discretization scheme in time, centered discretization scheme in space
#periodic boundary conditions
import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as sl
from matplotlib.animation import FuncAnimation
plt.ion()

def D(x):
  res = x**2 + 1
  return res

#define initial function
def u(x):
  x0=10.5
  sigma=10.0
  if (np.abs(x-x0)<=sigma):
    res=0.5*(1.0+np.cos(np.pi*(x-x0)/sigma))
  else:
    res=0.0
  return res

dx=0.4
L=np.arange(0.0,50.0,dx)
dt=10
c=0.01
du=0.0*np.ones(len(L))
du1=0.0*np.ones(len(L))
du2=0.0*np.ones(len(L))
n=len(L)-1

mu=0.5*dt*c/dx**2
if mu<=0.5:
  print "Courant number ", mu, " -> stable"
else:
  print "Courant number ", mu, " -> unstable"

for i in range(0,n):
  du[i]=u(L[i])

plt.plot(L, du)
#print du
#print mu

for t in np.arange(0.0,1000.1,dt):
  #predictor du1
  #mu=0.5*dt*c/dx**2
  du1[0]=mu*du[1] + (1.0-2.0*mu)*du[0] + mu*du[n]
  du1[n]=(1.0-2.0*mu)*du[n] + mu*du[n-1] + mu*du[0]
  for i in range(1,n-1):
    du1[i]=mu*du[i+1] + (1.0-2.0*mu)*du[i] + mu*du[i-1]
    #print du1[i]
  #corrector du
#solving three-diagonal matrix |C B 0 0 0 ... 0|
#                              |A C B 0 0 ... 0|
#                              |0 A C B 0 ... 0|
#                              |0     ...     0|
#                              |0 ... 0 0 A C B|
#                              |0 ... 0 0 0 A C|
# A=B=-mu, C=1+2*mu
#  beta=np.zeros(n+1)
#  alfa=np.zeros(n+1)
  #steps forward
#  beta[1]=-du1[0]/mu
#  alfa[1]=(1.0+2.0*mu)/mu
#  for i in range(0,n):
#    a=-mu*alfa[i]+1.0+2.0*mu
#    alfa[i+1]=mu/a
#    beta[i+1]=(du1[i]+mu*beta[i])/a
  #steps backward
#  du[n]=(du1[n]+mu*beta[n])/(1.0+2.0*mu-mu*alfa[n])
#  for i in range(n,0,-1):
#    du[i-1]=alfa[i]*du[i]+beta[i]

  M=np.zeros((n+1,n+1))
  for i in range(1,n):
      M[i,i-1]=-mu
      M[i,i]=1.0+2.0*mu
      M[i,i+1]=-mu
  M[0,0]=1.0+2.0*mu
  M[n,n]=1.0+2.0*mu
  M[0,1]=-mu
  M[n,n-1]=-mu
#  du = sl.solve_triangular(M,du1)
  du=du1.copy()

#  plt.clf()
  plt.ylim((0,1.0))
  plt.plot(L, du)
  plt.pause(0.1)
raw_input()
#print du
    
