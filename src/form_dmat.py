#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron, linalg

def is_Hermitian(matrix):
	return numpy.allclose(matrix.T, matrix)

def pprint(matrix):
	print numpy.array_str(matrix, precision=2, suppress_small=True)

def form_dmat(evec, dim12, dim34, m):
	'''
	For the reduced density matrix on the A* sub-block

	Returns:
		dmat	: density matrix
		matO	: rotation matrix
	'''
	
	projvec = numpy.reshape(evec, (dim12, dim34))
	dmat = projvec.dot(projvec.T)


	if dim12 > m:
		evals, evecs = numpy.linalg.eigh(dmat)
#	neig = m
#	evals, evecs = linalg.eigsh(dmat, k=neig, which="SA", tol=1e-08, maxiter=100000)
#	idx = evals.argsort()[::-1]   
#	evals = evals[idx]
#	evecs = evecs[:,idx]
	else:
		evals, evecs = numpy.linalg.eigh(dmat)

	
	if dim12 > m:
		matO = evecs[:,dim12-1-m+1:dim12-1]
	else:
		matO = evecs

	return dmat, matO

if __name__=="__main__":

	evec 	= numpy.array([-0.5, -0.5, -0.5, -0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5])
	evec 	= evec/numpy.linalg.norm(evec)
	dim12 	= 4
	dmat 	= numpy.zeros((dim12, dim12))

	dmat, matO = form_dmat(evec, dim12, dim12, 32)
