#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron
from sites import Site
from form_H12 import form_H12

if __name__=="__main__":
	
	site1 	= Site()
	site2 	= Site()

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

	'''
	Forming block H34
	'''

	site34 = Site()
	form_H12(sites[0], sites[1], site34)

	'''
	Forming block H1234
	'''

	site1234 = Site()
	form_H12(site12, site34, site1234)

	print site1234.H.toarray()
