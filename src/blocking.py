#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron, linalg
from sites import Site
from form_H12 import form_H12
from form_dmat import form_dmat

def is_Hermitian(matrix):
	return (matrix.transpose() == matrix).all()

def pprint(matrix):
	print numpy.array_str(matrix, precision=2, suppress_small=True)

if __name__=="__main__":
	
	site1 	= Site()
	site2 	= Site()

	nsites 	= 4

	for i in range(1,nsites):

		dim1 	= site1.H.shape[0]
		dim2 	= site2.H.shape[0]

		dim12	= dim1*dim2
		dim34	= dim1*dim2
		dim1234	= dim12*dim34

		if dim12 < 54:
			m 	= dim12
		else:
			m	= 54

		dimO 	= m

		sites = []

		sites.append(site1)
		sites.append(site2)

		'''
		Forming block H12
		'''

		site12 = Site()
		form_H12(sites[0], sites[1], site12)
		dim12 = site12.H.shape[0]

		'''
		Forming block H34
		'''

		site34 = Site()
		form_H12(sites[1], sites[0], site34)

		'''
		Forming block H1234
		'''

		site1234 = Site()
		form_H12(site12, site34, site1234)
		dim1234 = site1234.H.shape[0]

		'''
		diagonalize H1234
		'''

		neig = 1

		print is_Hermitian(site1234.H.toarray())
#	evals, evec = linalg.eigs(site1234.H, neig, which="SM", tol=1e-09, maxiter=400)
		evals, evecs = numpy.linalg.eigh(site1234.H.toarray())
		evec = evecs[::,0]

		'''
		calculate dmat
		'''
		
		dmat, matO = form_dmat(evec, dim12)
		print evals, site12.H.shape, dim12, matO.shape


		matOT 	= matO.transpose()
		site12.H = matOT.dot(site12.H.dot(matO))
		site12.Jz[0] = matOT.dot(site12.Jz[0].dot(matO))
		site12.Sp[0] = matOT.dot(site12.Sp[0].dot(matO))
		site12.Sm[0] = matOT.dot(site12.Sm[0].dot(matO))
		site12.Jz[1] = matOT.dot(site12.Jz[1].dot(matO))
		site12.Sp[1] = matOT.dot(site12.Sp[1].dot(matO))
		site12.Sm[1] = matOT.dot(site12.Sm[1].dot(matO))
		
		site1 = site12
		pprint(matO.dot(matOT))
		pprint(site1.H)
