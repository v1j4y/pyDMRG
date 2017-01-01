#!/usr/bin/env python

import numpy
from sites import Site
from scipy.sparse import csr_matrix, kron

def form_H12(site1, site2, site12):
	'''
	Blocking of two sites into one 
	'''

	dim1 = site1.H.shape[0]
	dim2 = site2.H.shape[0]
	print dim1, dim2

	row 	= numpy.array([0])
	column 	= numpy.array([0])
	data 	= numpy.array([0.])
	H12 = csr_matrix((data ,(row ,column )),shape=(dim1*dim2,dim1*dim2))
	H12 = 	    (1.0)*kron(site1.Jz, site2.Jz)
	H12 = H12 + (0.5)*kron(site1.Sp, site2.Sm)
	H12 = H12 + (0.5)*kron(site1.Sm, site2.Sp)
	H12 = H12 + (1.0)*kron(site1.H, numpy.eye(dim2))
	H12 = H12 + (1.0)*kron(numpy.eye(dim1), site2.H)
	site12.H = H12

if __name__=="__main__":

	sites = []
	site12 = Site()
	nsites = 2
	for i in range(nsites):
		sites.append(Site())

	form_H12(sites[0], sites[1], site12)

	print site12.H.toarray()
