#!/usr/bin/env python

from Lorenz63 import *
import plot_data 
import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation

model = Lorenz63(x0=[0.0,1.0,0.0],par=[10.0,28.0,8.0/3.0],dt=0.1, mem = 10)
for n in range(10):
	model.advance(tmax=2,tw=2)
	print model.t
	print model.x
data = model.get_history()
data_array = np.array(data['state'])
time = data['time']
print np.shape(data_array)
print "time=", time
print "data_array", data_array















"""
Setting Experiment details
exp_name_list = ["nature", "syserr", "parerr", "numerr"]

clim = np.empty((4,3))
mn   = np.empty_like(clim)
#Save state

for idx, exp_name in enumerate(exp_name_list):
	clim[idx,:] = np.load("clim_std_%s.npy"%(exp_name))
	mn[idx,:]   = np.load("clim_mn_%s.npy"%(exp_name))
	print exp_name
	print "Climatology = " , clim[idx]
	print "Mean        = " , mn[idx]

plt.subplot(121)
plt.imshow(clim)
plt.subplot(122)
plt.imshow(mn)

plt.colorbar()
"""
