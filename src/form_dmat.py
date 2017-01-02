#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron, linalg

def is_Hermitian(matrix):
	return (matrix.transpose() == matrix).all()

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

	print is_Hermitian(dmat)
	print dmat
#evals, evecs = linalg.eigs(dmat, k=neig, which="LM", tol=0)
	evals, evecs = numpy.linalg.eigh(dmat)
	print evals
	print evecs[0].dot(evecs[0])
	print evecs[1]

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
