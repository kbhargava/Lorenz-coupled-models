#!/usr/bin/env python

################################################################################

#	Author				: Kriti Bhargava
#	Created on			: Thu May 10 00:00:22 EDT 2018
#	Last Modified on	:

################################################################################
# An implementation of Lorenz, Edward Norton (1963). "Deterministic nonperiodic 
# flow".
# dx/dt = sigma * (y-x)
# dy/dt = x * ( rho - z ) - y
# dz/dt = x * y - beta * z
# The model can be integrated using the double approx. method described in L63 
# or Runge-Kutta 4 (RK4).
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
		self._history = {'time':[self.t],'state': np.array(self.x)[np.newaxis,:,:]}

	def advance(self,tmax=1,tw=1,writeout=True):
		"""
		Integrate the model state by tmax timesteps
		Saving every tw timestep
		"""
		for i in range(tmax):
			self._solver()
			if (writeout):
				if((i+1)%(tw) == 0):
					self.save_state()
	
	def save_state(self):
		"""
		Save time and state to history
		"""
		self._history['time'].append(self.t)
		state = np.array(self.x[np.newaxis,:,:])
		self._history['state'] = np.vstack([self._history['state'],state])

	def get_history(self):
		"""
		Return the dictionary containing list of times and corresponding states
		"""
		#state = (np.array(self._history['state'])).rehsape(
		self._history['state'] = (np.squeeze(self._history['state']))
		return self._history

	def get_state(self,settime=False):
		"""
		Return the current state
		"""
		if (settime):
			self.t=0
			self._history.clear()
			self._history = {'time':[self.t],'state': np.array(self.x)}
		return self.x	
	
	def _rk4(self):
		"""
		Runge-Kutta 4-stage integration scheme
		"""
		x_old = self.x
		k1 = self.dxdt()
		self.x = x_old +k1*self.dt/2.0
		k2 = self.dxdt()
		self.x = x_old+k2*self.dt/2.0
		k3 = self.dxdt()
		self.x = x_old+k3*self.dt
		k4 = self.dxdt()
		x_new = x_old + (self.dt/6)*(k1+2.0*k2+2.0*k3+k4)
		self.x = x_new
		self.t+=self.dt

	def _rk2(self):
		"""
		Runge-Kutta 2-stage integration scheme
		"""
		x_old = self.x
		k1 = self.dxdt()
		self.x = x_old +k1*self.dt
		k2 = self.dxdt()
		x_new = x_old + (self.dt/2)*(k1+k2)
		self.x = x_new
		self.t+=self.dt
	
	def dxdt(self):
		k = np.empty_like(self.x)
		F=self.F
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
	x0	= [0.0,1.0,0.0]
	
	dt = 0.01
	timesteps = 200

	model = Lorenz63(x0,par,dt,solver="RK4",F=[0,0,0],mem=1)

	model.advance(tmax=timesteps,tw=10,writeout=True)
	h = model.get_history()
	for t,s in zip(h['time'],h['state']):
			print("%04d % 05d % 05d % 05d" %(round(t*1e2), 10*s[0], 10*s[1], 10*s[2]))
