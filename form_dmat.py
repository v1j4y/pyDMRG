#!/usr/bin/env python

import numpy

def form_dmat(evec, dim12):
	'''
	For the reduced density matrix on the A* sub-block
	'''
	
	projvec = numpy.reshape(evec, (dim12, dim12))
	dmat = projvec.dot(projvec.transpose())
	return dmat

if __name__=="__main__":

	evec 	= numpy.array([0.5, 0.5, 0.5, 0.5])
	dim12 	= 2
	dmat 	= numpy.zeros((dim12, dim12))

	dmat = form_dmat(evec, dim12)

	print dmat
