#diffusion advection equation du/dt=c*d2u/dx2 + k*du/dx in 1D
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
  x0=10.5
  sigma=10.0
  if (np.abs(x-x0)<=sigma):
    res=0.5*(1.0+np.cos(np.pi*(x-x0)/sigma))
  else:
    res=0.0
  return res

dx=0.4
L=np.arange(0.0,50.0,dx)
dt=0.5
c=0.5
k=-0.3
du=0.0*np.ones(len(L))
du1=0.0*np.ones(len(L))
du2=0.0*np.ones(len(L))
n=len(L)-1

la=0.25*dt*k/dx
mu=0.5*dt*c/dx**2
print np.max(np.abs(la),np.abs(mu))

for i in range(0,n-1):
  du[i]=u(L[i])

plt.plot(L, du)
#print du
#print mu

for t in np.arange(0.0,30.1,dt):
  #predictor du1
  #mu=0.5*dt*c/dx**2
  du1[0]=(mu + la)*du[1] + (1.0-2.0*mu)*du[0] + (mu + la)*du[n]
  du1[n]=(1.0-2.0*mu)*du[n] + (mu + la)*du[n-1] + (mu + la)*du[0]
  for i in range(1,n-1):
    du1[i]=(mu + la)*du[i+1] + (1.0-2.0*mu)*du[i] + (mu + la)*du[i-1]
    #print du1[i]
  #corrector du
#solving three-diagonal matrix |C B 0 0 0 ... 0|
#                              |A C B 0 0 ... 0|
#                              |0 A C B 0 ... 0|
#                              |0     ...     0|
#                              |0 ... 0 0 A C B|
#                              |0 ... 0 0 0 A C|
  a=-mu-la
  c=1.0+2.0*mu
  b=-mu+la
  beta=np.zeros(n+1)
  alfa=np.zeros(n+1)
  #steps forward
  beta[1]=du1[0]/c
  alfa[1]=-b/c
  for i in range(0,n):
    dev=a*alfa[i]+c
    alfa[i+1]=-b/dev
    beta[i+1]=(du1[i]-a*beta[i])/dev
  #steps backward
  du[n]=(du1[n]-a*beta[n])/(c+a*alfa[n])
  for i in range(n,0,-1):
    du[i-1]=alfa[i]*du[i]+beta[i]
  plt.clf()
  plt.ylim((0,1.0))
  plt.plot(L, du)
  plt.pause(0.1)
raw_input()
