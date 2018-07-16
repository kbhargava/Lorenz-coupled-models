#!/usr/bin/env python

from Lorenz63 import *

import matplotlib.pyplot as plt
from mpl_toolkits.mplot3d import Axes3D
from matplotlib.colors import cnames
from matplotlib import animation
"""
Setting Parameters
"""
par = [10.0,28.0,8.0/3.0] # values of sigma, rho and beta

"""
Initial Value
"""
x0 = [0.0,1.0,0.0]

"""
Timesteps and max timestep information
"""
dt =0.01				 # timestep
t_spinup = 1000			 # timesteps for spinning up
timesteps = 50			 # max timesteps 1 time unit = 5 days

"""
Model initialization and spin up
"""
model_nature = Lorenz63(x0,par,dt,solver="RK4",F=0,mem=1)

model_nature.advance(tmax=t_spinup,tw=1,writeout=False) #spinup without creating history
x0 = model_nature.get_state(settime=True)
print (model_nature.t)

model_rk2	 = Lorenz63(x0,par,dt,solver="RK2",F=0,mem=1) 

"""
Running nature run and rk2 model
"""
model_nature.advance(tmax=timesteps) 
model_rk2.advance(tmax=timesteps) 

"""
Get the model states
"""
h_nature = model_nature.get_history()
h_rk2	 = model_rk2.get_history()
states_nature = np.array(h_nature['state'])
states_rk2 = np.array(h_rk2['state'])

"""
Calculate Errors
"""
state_error = states_rk2-states_nature

#plt.plot(h_nature['time'],states_nature[:,0], 'b-', label =" x_nature")
#plt.plot(h_nature['time'],states_rk2[:,0], 'r--', label =  " x_rk2")
#plt.plot(h_nature['time'],states_nature[:,1], 'r-', label =" y_nature")
#plt.plot(h_nature['time'],states_rk2[:,1], 'r--', label =  " y_rk2")
#plt.plot(h_nature['time'],states_nature[:,2], 'm-', label =" z_nature")
#plt.plot(h_nature['time'],states_rk2[:,2], 'm--', label =  " z_rk2")

plt.plot(h_nature['time'],state_error[:,0], 'black', label =" x_error")
plt.plot(h_nature['time'],state_error[:,1], 'r-', label =" y_nature")
plt.plot(h_nature['time'],state_error[:,2], 'm-', label =" z_nature")
plt.legend()
plt.savefig("err.png",dpi=300)
plt.show()
plt.close()
"""
fig = plt.figure()
ax = fig.add_subplot(111, projection='3d')

#ax.plot(states_nature[:,0], states_nature[:,1], states_nature[:,2],'b',zdir ='z', label ='Nature run Lorenz63')
#ax.plot(states_rk2[:,0], states_rk2[:,1], states_rk2[:,2],'r',zdir='z', label ='RK2 run Lorenz63')
ax.plot(state_error[:,0], state_error[:,1], state_error[:,2],'g',zdir='y', label ='Error in Lorenz63')
ax.set_xlabel('X axis')
ax.set_ylabel('Y axis')
ax.set_zlabel('Z axis')
ax.legend()
plt.show()
"""
"""
N_trajectories = 1#3 # 3 trajectories for nature,rk2 and error

# Set up figure & 3D axis for animation
fig = plt.figure()
ax = fig.add_axes([0,0,1,1], projection='3d')
ax.axis('off')

# Choose color for each trajectory
colors = plt.cm.jet(np.linspace(0,1,N_trajectories))

# Set up lines and points
lines = sum([ax.plot([], [], [], '-', c=c) for c in colors], [])
pts = sum([ax.plot([], [], [], 'o', c=c) for c in colors], [])

# prepare the axes limits
ax.set_xlim((-25, 25))
ax.set_ylim((-35, 35))
ax.set_zlim((5, 55))

# set point-of-view: specified by (altitude degrees, azimuth degrees)
ax.view_init(30, 0)

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
	i = (2 * i) % x_t.shape[1]
	
	for line, pt, xi in zip(lines, pts, ):
		x, y, z = xi[:i].T
		line.set_data(x, y)
		line.set_3d_properties(z)
		
		pt.set_data(x[-1:], y[-1:])
		pt.set_3d_properties(z[-1:])
	
	ax.view_init(30, 0.3 * i)
	fig.canvas.draw()
	return lines + pts

# instantiate the animator.
anim = animation.FuncAnimation(fig, animate, init_func=init,frames=500, interval=30, blit=True)
															
# Save as mp4. This requires mplayer or ffmpeg to be installed
#anim.save('lorentz_attractor.mp4', fps=15, extra_args=['-vcodec', 'libx264'])
plt.show()
"""
