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
import numpy.matlib
class DA:

	def __init__(self,x0,Xb=[],obs=[],H=[],nens=1,obserr=[],method="ETKF"):
		self.x0     =x0
		self.ndim   = np.size(x0)
		self.nens   = nens
		self.Xb     = np.empty((self.ndim,nens))
		self.obs    = obs
		self.R      = np.matrix(obserr)
		self.method = method
		self.H      = H

	def init_ens(self):
		"""
		Initialize the ensemble members
		"""
		self.Xb = self.x0+np.matmul(np.random.normal(size=np.shape(self.Xb)),R)
		return self.Xb
	
	def etkf(self, Xb,obs):
		"""
		Apply ETKF
		"""
		nens = self.nens
		R    = self.R

		Xb = np.matrix(Xb)
		Yo = np.matrix(obs)

		Yb = np.matrix(np.empty_like(Xb))
		for i in range(nens):
			Yb[:,i] = np.dot(self.H,Xb[:,i])
		xmean = np.mean(Xb,axis=1)
		ymean = np.mean(Yb,axis=1)

		Xb    = Xb - np.matlib.repmat(xmean, 1, self.nens)
		Yb    = Yb - np.matlib.repmat(ymean, 1, self.nens)

		# Calculate C= Yb^T*R^(-1)
		Rinv = np.linalg.inv(R)
		C    = np.dot(np.transpose(Yb),Rinv)
		
		# Calculate Pa = [(k-1)I/rho +C*Yb]^(-1)
		I         = np.eye(nens)
		rho       = 1.0 # no inflation
		eigArg    = (nens-1)*I/rho + np.dot(C,Yb)
		lamda, P  = np.linalg.eigh(eigArg)
		Linv      = np.diag(1.0/lamda)
		PLinv     = np.dot(P,Linv)
		Pa        = np.dot(PLinv, np.transpose(P))

		# Calculate Wa = [(k-1)Pa]^(1/2)
		Linvsqrt  = np.diag(1/np.sqrt(lamda))
		PLinvsqrt = np.dot(P, Linvsqrt)
		Wa        = np.sqrt((nens-1))*np.dot(PLinvsqrt,np.transpose(P))

		# Compute the mean update wma = Pa*C*yo-ym
		d         = Yo -ymean
		Cd        = np.dot(C,d)
		wma       = np.dot(Pa,Cd)
		print(np.shape(wma))
		Wa        = Wa + np.matlib.repmat(wma, 1, nens)
		Xa        = np.dot(Xb,Wa)+ np.matlib.repmat(xmean, 1, nens)
		return Xa

