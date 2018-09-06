#!/usr/bin/env python
################################################################################
#                         EXPERIMENT PARERR SETTING
################################################################################
"""
Setting Experiment details
"""
exp_name = "parerr"
nens      = 3 #nens

"""
Setting Parameters
"""
par   = [10.0, 29.4,8.0/3.0] # values of sigma, rho and beta
F0    = [0,0,0]

"""
Timesteps and max timestep information
"""
dt            = 0.01         # timestep
t_spinup      = 0 #100000       # timesteps for spinning up
tmax          = 100       # max time in Model Time units
timesteps     = int(tmax/dt)     # max timesteps 1 time unit = 5 days
taw           = 1
"""

Path to Initial Value
"""
IC_file  = "data/IC/IC.npy"

"""
Make observations
"""
makeobs  = True

"""
Setting DA parameters
"""
do_DA    =  False
taw      =  1 # Assimilation window


"""
