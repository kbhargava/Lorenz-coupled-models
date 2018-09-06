#!/usr/bin/env python

################################################################################

#  Author            : Kriti Bhargava
#  Created on        : Thu May 10 00:00:22 EDT 2018
#  Last Modified on  : Thu Sep  6 12:24:51 EDT 2018


################################################################################
# An implementation of Lorenz, Edward Norton (1963). "Deterministic nonperiodic 
# flow".
# dx/dt = sigma * (y-x)
# dy/dt = x * ( rho - z ) - y
# dz/dt = x * y - beta * z
# The model can be integrated using the double approx. method described in L63 
# or Runge-Kutta 4 (RK4).

# This is a simplified version without any save_state or get_history feature
################################################################################

import numpy as np

class Lorenz63:

  def __init__(self,x0,par,dt,solver="RK4",F=[0,0,0],mem=1):
    self.x = x0+0.1*np.random.normal(size=(mem,np.size(x0)))
    self.dt = dt
    self.par = par
    self.t=0
    self.F=F
    self.mem =mem
    if solver == "RK2":
      self._solver = self._rk2
    else:
      self._solver = self._rk4

  def advance(self,tmax=1,tw=1,F=[0,0,0]):
    """
    Integrate the model state by tmax timesteps
    Saving every tw timestep
    """
    for i in range(tmax):
      self._solver(F)
    return self.x,self.t  
  
  def _rk4(self,F):
    """
    Runge-Kutta 4-stage integration scheme
    """
    x_old = self.x
    k1 = self.dxdt(F)
    self.x = x_old +k1*self.dt/2.0
    k2 = self.dxdt(F)
    self.x = x_old+k2*self.dt/2.0
    k3 = self.dxdt(F)
    self.x = x_old+k3*self.dt
    k4 = self.dxdt(F)
    x_new = x_old + (self.dt/6)*(k1+2.0*k2+2.0*k3+k4)
    self.x = x_new
    self.t+=self.dt

  def _rk2(self,F):
    """
    Runge-Kutta 2-stage integration scheme
    """
    x_old = self.x
    k1 = self.dxdt(F)
    self.x = x_old +k1*self.dt
    k2 = self.dxdt(F)
    x_new = x_old + (self.dt/2)*(k1+k2)
    self.x = x_new
    self.t+=self.dt
  
  def dxdt(self,F):
    k = np.empty_like(self.x)
    F=F+self.F
    x,y,z = self.x.T
    sigma, rho, beta = self.par
    k[:,0] = sigma*(y-x)+F[0]
    k[:,1] = rho*x-y-x*z+F[1]
    k[:,2] = x*y -beta*z+F[2]
    return k

################################################################################    

################################################################################    
if __name__ == '__main__':
  par = [10.0,28.0,8.0/3.0]
  x0  = [0.0,1.0,0.0]
  
  dt = 0.01
  timesteps = 200
  x=np.zeros((timesteps,3))
  t=np.zeros((timesteps))
  model = Lorenz63(x0,par,dt,solver="RK4",F=[0,0,0],mem=1)
  for i in range(timesteps):
    x[i,:],t[i]= model.advance()
  for t,s in zip(t,x):
      print("%04d % 05d % 05d % 05d" %(round(t*1e2), s[0], s[1], s[2]))
