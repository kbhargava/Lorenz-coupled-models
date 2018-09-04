#!/usr/bin/env python

################################################################################

#	Author				: Kriti Bhargava
#	Created on			: Wed Jul 25 22:53:50 EDT 2018 
#	Last Modified on	:

################################################################################
# An implementation of Different Data Assimilation Schemes
# 
# 1. ETKF based on Hunt et al. 2007 
# 
# 
# 
# 
################################################################################

import numpy as np

class obs_net:

	def __init__(self,x_truth,obserr,nobs):
		self.nt     = np.shape(x_truth)[0]
		self.ndim   = np.shape(x_truth)[1]
		self.nobs   = nobs
		self.x_truth= x_truth
		self.obs    = np.zeros((self.nt, nobs))
		self.R      = np.diag(obserr)
		self.H      = np.zeros((nobs,self.ndim))

	def makeobs(self):
		"""
		Makeobs
		"""
		err       = np.random.normal(size=np.shape(self.x_truth))
		self.getH()	
		H         = self.H
		self.R    = np.sqrt(np.dot(np.dot(np.square(self.R),H),np.transpose(H))) 
		for i in range(self.nobs):
			self.obs[:,i] = self.x_truth[:,i] +err[:,i]*self.R[i,i]
		return self.obs, self.H, self.R
	
	def getH(self):
		"""
		Calculate H
		"""
		
		for n in range(self.nobs):
			self.H[n,n] = 1
