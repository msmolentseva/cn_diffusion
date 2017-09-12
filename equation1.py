#diffusion equation du/dt=c*d2u/dx2 in 1D
#Crank-Nicolson discretization scheme in time, centered discretization scheme in space
#periodic boundary conditions
import numpy as np
import matplotlib.pyplot as plt
plt.ion()

def D(x):
  res = x**2 + 1
  return res

#define initial function
def u(x):
  x0=0.7
  sigma=0.5
  if (np.abs(x-x0)<=sigma):
    res=0.5*(1.0+np.cos(np.pi*(x-x0)/sigma))
  else:
    res=0.0
  return res

dx=0.5
L=np.arange(0.0,50.0,dx)
dt=0.5
c=0.2
du=0.0*np.ones(len(L))
du1=0.0*np.ones(len(L))
n=len(L)-1

mu=0.5*dt*c/dx**2
if mu<=0.5:
  print "Courant number ", mu, " -> stable"
else:
  print "Courant number ", mu, " -> unstable"

for i in range(0,n-1):
  du[i]=u(L[i])

#plt.plot(L, du)


for t in np.arange(0.0,30.0,dt):
  #predictor du1
  du1[0]=mu*du[1] + (1.0-2.0*mu)*du[0] + mu*du[n-1]
  du1[n]=(1.0-2.0*mu)*du[n] + mu*du[n-1] + mu*du[0]
  for i in range(1,n-1):
    x=dx*(i+1)
    du1[i]=mu*du1[i+1] + (1.0-2.0*mu)*du1[i] + mu*du1[i-1]
  
  #corrector du
#solving three-diagonal matrix |C B 0 0 0 ... 0|
#                              |A C B 0 0 ... 0|
#                              |0 A C B 0 ... 0|
#                              |0     ...     0|
#                              |0 ... 0 0 A C B|
#                              |0 ... 0 0 0 A C|
# A=B=-mu, C=1+2*mu
  beta=np.zeros(n+1)
  alfa=np.zeros(n+1)
  #steps forward
  beta[1]=-du1[0]/mu
  alfa[1]=(1+2*mu)/mu
  for i in range(0,n-1):
    a=-mu*alfa[i]+1+2*mu
    alfa[i+1]=mu/a
    beta[i+1]=(du1[i]+mu*beta[i])/a
  #steps backward
  du[n]=(du1[n]+mu*beta[n])/(1+2*mu-mu*alfa[n])
  for i in range(n,0,-1):
    du[i-1]=alfa[i]*du[i]+beta[i]

  plt.clf()
  plt.ylim((0,0.05))
  plt.plot(L, du)
  plt.pause(0.1)
raw_input()
    
