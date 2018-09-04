#!/usr/bin/env python

from Lorenz63 import *
from makeobs import *
import plot_data 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation
"""
Setting Experiment details
"""
exp_name_list = ["nature", "syserr", "parerr", "numerr", "obs"]

nens          = 3 #nens
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
Timesteps and max timestep information
"""
dt            = 0.01				 # timestep
dt_err        = 0.02				 # 
t_spinup      = 0 #100000			 # timesteps for spinning up
tmax          = 100			 # max time in Model Time units
timesteps     = int(tmax/dt)		 # max timesteps 1 time unit = 5 days
timesteps_err = int(tmax/dt_err)      
"""
Loading truth and Initial Value
"""
x0      = np.load("IC.npy")#[0.0,1.0,0.0]
std     ={}
for exp_name in exp_name_list:
	std[exp_name] = std["nature"]*0.1 if exp_name=="obs" else np.load("clim_std_%s.npy" %exp_name)
x0 = x0+np.random.normal(size=3)
"""
Make observations
"""
x_truth = np.load("nature_run_nature.npy")
obs_data = obs_net(x_truth,std["obs"],3)
yo,H,R = obs_data.makeobs()
"""
Model initialization and spin up
"""
model_nature = Lorenz63(x0,par,dt,solver="RK4",F=F0)
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
np.save("IC_perturb",x0)
data = [states_nature, states_syserr, states_parerr, states_numerr, yo] 
"""
Save state
"""
for idx, exp_name in enumerate(exp_name_list):
	np.save("obs_10" if exp_name == "obs" else "free_run_%s"%(exp_name), data[idx] )
plot_data.anim_L63(data,exp_name_list)

"""
Calculate Errors
"""
exp_num =np.size(exp_name_list)
rmse = np.zeros((timesteps_err+1,exp_num))
rmse[:,0] = np.sqrt(np.mean((x_truth-states_nature)**2,axis=1))
rmse[:,1] = np.sqrt(np.mean((x_truth-states_syserr)**2,axis=1))
rmse[:,2] = np.sqrt(np.mean((x_truth-states_parerr)**2,axis=1))
rmse[:,3] = np.sqrt(np.mean((x_truth-states_numerr)**2,axis=1))
rmse[:,4] = np.sqrt(np.mean((x_truth-yo           )**2,axis=1))

plot_data.plot2d(rmse,time,exp_name_list,"RMSE with perturbed IC")
