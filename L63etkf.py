#!/usr/bin/env python

from Lorenz63 import *
from makeobs import *
from DA import *
import plot_data 
import matplotlib.pyplot as plt
import sys
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation
"""
Setting Experiment details
"""
exp_name_list = ["nature", "syserr", "parerr", "numerr", "obs"]

nens		  = 3 #nens
"""
Setting Parameters
"""
f	   = 10
theta   = 3*np.pi/4

par	 = [10.0,28.0,8.0/3.0] # values of sigma, rho and beta
par_err = [10.0, 29.4,8.0/3.0]

F0	  = [0,0,0]
F_err   = [f*np.cos(theta), f*np.sin(theta),0]

"""
Timesteps and max timestep information
"""
dt		      = 0.01				 # timestep
dt_err	      = 0.02				 # 
t_spinup	  = 0 #100000			 # timesteps for spinning up
tmax		  = 100			 # max time in Model Time units
timesteps	  = int(tmax/dt)		 # max timesteps 1 time unit = 5 days
timesteps_err = int(tmax/dt_err)	  
taw           = 2
"""
Loading truth and Initial Value
"""
x0	  = np.load("IC_perturb.npy")#[0.0,1.0,0.0]
std	 ={}
for exp_name in exp_name_list:
	std[exp_name] = std["nature"]*0.1 if exp_name=="obs" else np.load("clim_std_%s.npy" %exp_name)

x_truth = np.load("nature_run_nature.npy")
"""
Make observations
"""
obs_data = obs_net(x_truth,std["obs"],3)
yo,H,R = obs_data.makeobs()
"""
Model initialization 
"""
model_nature	 = Lorenz63(x0,par	,dt	,solver="RK4",F=F0   ,mem=nens)
model_syserr	 = Lorenz63(x0,par	,dt	,solver="RK4",F=F_err,mem=nens) 
model_parerr	 = Lorenz63(x0,par_err,dt	,solver="RK4",F=F0   ,mem=nens) 
model_numerr	 = Lorenz63(x0,par	,dt_err,solver="RK4",F=F0   ,mem=nens) 

"""
Running DA system

"""
time=[]
x_fcst_nature= np.zeros((int(timesteps_err+1/taw),3))
x_fcst_syserr= np.zeros((int(timesteps_err+1/taw),3))
x_fcst_parerr= np.zeros((int(timesteps_err+1/taw),3))
x_fcst_numerr= np.zeros((int(timesteps_err+1/taw),3))
x_anal_nature= np.zeros((int(timesteps_err+1/taw),3))
x_anal_syserr= np.zeros((int(timesteps_err+1/taw),3))
x_anal_parerr= np.zeros((int(timesteps_err+1/taw),3))
x_anal_numerr= np.zeros((int(timesteps_err+1/taw),3))
obs          = np.zeros((int(timesteps_err+1/taw),3))
das = DA(x0,obs=yo,nens=nens,H=H,obserr=R,method="ETKF")
for n in range(timesteps_err):
	model_nature.advance(tmax=2,writeout=False) 
	model_syserr.advance(tmax=2,writeout=False) 
	model_parerr.advance(tmax=2,writeout=False) 
	model_numerr.advance(writeout=False) 
	if ( n%(taw) ==0 ):
		"""
		Saving truth, forecast and obs at analysis step
		"""
		print (n)
		time.append(model_nature.t)
		x_fcst_nature[n,:] = np.mean(model_nature.x,axis=0)
		x_fcst_syserr[n,:] = np.mean(model_syserr.x,axis=0)
		x_fcst_parerr[n,:] = np.mean(model_parerr.x,axis=0)
		x_fcst_numerr[n,:] = np.mean(model_numerr.x,axis=0)
		obs[n]=yo[2*n]
		x_truth_anal= x_truth[2*n]
		"""
		DA STEP
		"""
		model_nature.x = np.transpose(das.etkf(Xb=model_nature.x.T,obs=yo[2*n]))
		model_syserr.x = np.transpose(das.etkf(Xb=model_syserr.x.T,obs=yo[2*n]))
		model_parerr.x = np.transpose(das.etkf(Xb=model_parerr.x.T,obs=yo[2*n]))
		model_numerr.x = np.transpose(das.etkf(Xb=model_numerr.x.T,obs=yo[2*n]))
		x_anal_nature[n,:] = np.mean(model_nature.x,axis=0)
		x_anal_syserr[n,:] = np.mean(model_syserr.x,axis=0)
		x_anal_parerr[n,:] = np.mean(model_parerr.x,axis=0)
		x_anal_numerr[n,:] = np.mean(model_numerr.x,axis=0)
	model_nature.save_state()
	model_syserr.save_state()
	model_parerr.save_state()
	model_numerr.save_state()



	
data = [x_anal_nature, x_anal_syserr, x_anal_parerr, x_anal_numerr, obs] 
data_fcst = [x_fcst_nature, x_fcst_syserr, x_fcst_parerr, x_fcst_numerr, obs] 
"""
Save state
"""
for idx, exp_name in enumerate(exp_name_list):
	np.save("obs_da" if exp_name == "obs" else "anal_run_%s"%(exp_name), data[idx] )
plot_data.anim_L63(data,exp_name_list)

"""
Calculate Errors
"""
exp_num =np.size(exp_name_list)
rmse = np.zeros((timesteps_err+1,exp_num))
rmse[:,0] = np.sqrt(np.mean((x_truth_anal-x_anal_nature)**2,axis=1))
rmse[:,1] = np.sqrt(np.mean((x_truth_anal-x_anal_syserr)**2,axis=1))
rmse[:,2] = np.sqrt(np.mean((x_truth_anal-x_anal_parerr)**2,axis=1))
rmse[:,3] = np.sqrt(np.mean((x_truth_anal-x_anal_numerr)**2,axis=1))
rmse[:,4] = np.sqrt(np.mean((x_truth_anal-obs		  )**2,axis=1))

plot_data.plot2d(rmse,time,exp_name_list,"RMSE of analysis with \n models using ETKF",figname="RMSE_anal.png")

label_list = ["Forceast", "Analysis", "Observatios"]
exp_num =np.size(label_list)
rmse = np.zeros((timesteps_err+1,exp_num))
rmse[:,0] = np.sqrt(np.mean((x_truth_anal-x_fcst_nature)**2,axis=1))
rmse[:,1] = np.sqrt(np.mean((x_truth_anal-x_anal_nature)**2,axis=1))
rmse[:,2] = np.sqrt(np.mean((x_truth_anal-obs		  )**2,axis=1))

plot_data.plot2d(rmse,time,exp_name_list,"RMSE of analysis with perfect \nmodel using ETKF",figname="RMSE_nature_DA.png")
