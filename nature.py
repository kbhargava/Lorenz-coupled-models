#!/usr/bin/env python

import sys
from Lorenz63 import *
import plot_data 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation
"""
Setting Experiment details
"""
exp_name_list = ["nature", "syserr", "parerr", "numerr"]
"""
Setting Parameters
"""
f       = 10
theta   = 3*np.pi/4

par     = [10.0,28.0,8.0/3.0] # values of sigma, rho and beta
par_err = [10.0, 29.4,8.0/3.0]

F0      = [0,0,0]
F_err   = [f*np.cos(theta), f*np.sin(theta),0]

"""
Initial Value
"""
x0 = np.load("IC.npy")#[0.0,1.0,0.0]
"""
Timesteps and max timestep information
"""
dt            = 0.01				 # timestep
dt_err        = 0.02				 # 
t_spinup      = 0 #100000			 # timesteps for spinning up
tmax          = 100		 # max time in Model Time units
timesteps     = int(tmax/dt)		 # max timesteps 1 time unit = 5 days
timesteps_err = int(tmax/dt_err)      
"""
Model initialization and spin up
"""
model_nature = Lorenz63(x0,par,dt,solver="RK4",F=F0)

model_nature.advance(tmax=t_spinup,tw=1,writeout=False) #spinup without creating history
x0 = model_nature.get_state(settime=True)

model_syserr	 = Lorenz63(x0,par    ,dt    ,solver="RK4",F=F_err) 
model_parerr	 = Lorenz63(x0,par_err,dt    ,solver="RK4",F=F0   ) 
model_numerr	 = Lorenz63(x0,par    ,dt_err,solver="RK4",F=F0   ) 

"""
Running nature run and rk2 model
"""
model_nature.advance(tmax=timesteps,tw=2) 
model_syserr.advance(tmax=timesteps,tw=2) 
model_parerr.advance(tmax=timesteps,tw=2) 
model_numerr.advance(tmax=timesteps_err,tw=1) 

"""
Get the model states
"""
h_nature = model_nature.get_history()
h_syserr = model_syserr.get_history()
h_parerr = model_parerr.get_history()
h_numerr = model_numerr.get_history()

states_nature = np.array(h_nature['state'])
states_syserr = np.array(h_syserr['state'])
states_parerr = np.array(h_parerr['state'])
states_numerr = np.array(h_numerr['state'])
time          = np.array(h_numerr['time'])

#x0 = model_nature.get_state(settime=False)
#np.save("IC",x0)
data = [states_nature, states_syserr, states_parerr, states_numerr] 
"""
Save state
"""
for idx, exp_name in enumerate(exp_name_list):
	np.save("nature_run_%s"%(exp_name), data[idx] )
#	np.save("clim_std_%s"%(exp_name), np.std(data[idx],axis=0,ddof=1 ))
#	np.save("clim_mn_%s"%(exp_name), np.mean(data[idx],axis=0 ))

"""
Calculate Errors
"""

plot_data.anim_L63(data,exp_name_list)
print(time)
rmse = np.zeros(np.shape(states_nature))
rmse[:,0] = np.sqrt(np.mean((states_nature-states_syserr)**2,axis=1))
rmse[:,1] = np.sqrt(np.mean((states_nature-states_parerr)**2,axis=1))
rmse[:,2] = np.sqrt(np.mean((states_nature-states_numerr)**2,axis=1))

plot_data.plot2d(rmse, time,exp_name_list[1:],"RMSE with perfect IC")

