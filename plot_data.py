#!/usr/bin/env python

import matplotlib.pyplot as plt
import numpy as np
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation

def anim_L63(data,label,figname="Lorenz.mp4"):

	N_trajectories = np.shape(data)[0]#3 # 3 trajectories for nature,syserr and error
	xi =data
	# Set up figure & 3D axis for animation
	fig = plt.figure()
	ax = fig.add_axes([0,0,1,1], projection='3d')
	ax.axis('off')
	
	# Choose color for each trajectory
	colors = plt.cm.jet(np.linspace(0,1,N_trajectories))
	
	# Set up lines and points
	lines = sum([ax.plot([], [], [], '-', c=c,label=label[idx]) for idx,c in enumerate(colors)], [])
	pts = sum([ax.plot([], [], [], 'o', c=c) for c in colors], [])
	
	# prepare the axes limits
	ax.set_xlim((-25, 25))
	ax.set_ylim((-35, 35))
	ax.set_zlim((5, 55))
	
	# set point-of-view: specified by (altitude degrees, azimuth degrees)
	ax.view_init(30, 0)
	legend=plt.legend()	
	# initialization function: plot the background of each frame
	def init():
		for line, pt in zip(lines, pts):
			line.set_data([], [])
			line.set_3d_properties([])
	
			pt.set_data([], [])
			pt.set_3d_properties([])
		return lines + pts
	
	# animation function.  This will be called sequentially with the frame number
	def animate(i):
	# we'll step two time-steps per frame.  This leads to nice results.
		i = (2 * i) #% states_nature.shape[1]
		
		for line, pt, xi in zip(lines, pts, data):
			#x, y, z = xi[:i,0]#.T
			x = xi[:i,0]#.T
			y = xi[:i,1]#.T
			z = xi[:i,2]#.T
			line.set_data(x, y)
			line.set_3d_properties(z)
			
			pt.set_data(x[-1:], y[-1:])
			pt.set_3d_properties(z[-1:])
		
		ax.view_init(30, 0.3 * i)
		fig.canvas.draw()
		fig.suptitle("Timesteps =%s" %(i))
		return lines + pts
	
	# instantiate the animator.
	anim = animation.FuncAnimation(fig, animate, init_func=init,frames=500, interval=30, blit=True)
																
	# Save as mp4. This requires mplayer or ffmpeg to be installed
	anim.save(figname, fps=15)#, extra_args=['-vcodec', 'libx264'])
	plt.close()
	#plt.show()
	#"""

def plot2d(data,time, exp_name,title,figname="RMSE_nature.png"):
	nt, exp_num = np.shape(data)
	colors = plt.cm.brg(np.linspace(0,1,exp_num))
	fig = plt.figure(figsize=(10,20))
	for n,c,exp in zip(range(exp_num),colors,exp_name):
		plt.semilogy(time,data[:,n],color=c, label=exp)
	plt.legend()
	plt.ylabel("RMSE")
	plt.xlabel(" Time in MTU ")
	fig.suptitle(title)
#	fig.savefig("RMSE_nature.png", dpi=300)
	fig.savefig(figname, dpi=300)
	plt.show()


