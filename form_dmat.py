#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron, linalg

def form_dmat(evec, dim12):
	'''
	For the reduced density matrix on the A* sub-block

	Returns:
		dmat	: denaity matrix
		matO	: rotation matrix
	'''
	
	projvec = numpy.reshape(evec, (dim12, dim12))
	dmat = projvec.dot(projvec.transpose())

	neig = len(dmat)-2

	evals, evecs = linalg.eigs(dmat, neig, which="LM", tol=0)

	if dim12 <= 16:
		matO = evecs
	else:
		matO = evecs[::,0::31]

	return dmat, matO

if __name__=="__main__":

	evec 	= numpy.array([-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
	evec 	= evec/numpy.linalg.norm(evec)
	dim12 	= 4
	dmat 	= numpy.zeros((dim12, dim12))

	dmat, matO = form_dmat(evec, dim12)

	print matO
