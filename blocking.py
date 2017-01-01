#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron, linalg
from sites import Site
from form_H12 import form_H12
from form_dmat import form_dmat

if __name__=="__main__":
	
	site1 	= Site()
	site2 	= Site()

	nsites 	= 3

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

		evals, evec = linalg.eigs(site1234.H, neig, which="SR", tol=0)

		'''
		calculate dmat
		'''
		
		dmat, matO = form_dmat(evec, dim12)

		print  evals, dim1

		matOT 	= matO.transpose()
		site12.H = matOT.dot(site12.H.dot(matO))
		site12.Jz[0] = matOT.dot(site12.Jz[0].dot(matO))
		site12.Sp[0] = matOT.dot(site12.Sp[0].dot(matO))
		site12.Sm[0] = matOT.dot(site12.Sm[0].dot(matO))
		site12.Jz[1] = matOT.dot(site12.Jz[1].dot(matO))
		site12.Sp[1] = matOT.dot(site12.Sp[1].dot(matO))
		site12.Sm[1] = matOT.dot(site12.Sm[1].dot(matO))
		
		site1 = site12
