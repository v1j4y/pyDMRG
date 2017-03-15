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

	row 	= numpy.array([0])
	column 	= numpy.array([0])
	data 	= numpy.array([0.])
	H12 	= csr_matrix((data ,(row ,column )),shape=(dim1*dim2,dim1*dim2))
	H12 	= 	    (1.0)*kron(site1.Jz[1], site2.Jz[0])
	H12 	= H12 + (0.5)*kron(site1.Sp[1], site2.Sm[0])
	H12 	= H12 + (0.5)*kron(site1.Sm[1], site2.Sp[0])
	H12 	= H12 + (1.0)*kron(site1.H, numpy.eye(dim2))
	H12 	= H12 + (1.0)*kron(numpy.eye(dim1), site2.H)

	Jz12	= 	    (1.0)*kron(site1.Jz[0], numpy.eye(dim2))
	site12.Jz[0] = Jz12
	Jz12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Jz[1])
	site12.Jz[1] = Jz12

	Sp12	= 	    (1.0)*kron(site1.Sp[0], numpy.eye(dim2))
	site12.Sp[0] = Sp12
	Sp12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Sp[1])
	site12.Sp[1] = Sp12

	Sm12	= 	    (1.0)*kron(site1.Sm[0], numpy.eye(dim2))
	site12.Sm[0] = Sm12
	Sm12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Sm[1])
	site12.Sm[1] = Sm12

	site12.H = H12

def form_H12_Hubb(site1, site2, site12):
	'''
	Blocking of two sites into one for the Hubbard Hamiltonian
	Hubb = Sum_i,j[ t_ij (c^dag_a c_b + c^dag_b c_a) ] + Sum_i [ U (n_i,a n_i,b) ]
	'''

	dim1 = site1.H.shape[0]
	dim2 = site2.H.shape[0]

	row 	= numpy.array([0])
	column 	= numpy.array([0])
	data 	= numpy.array([0.])
	H12 	= csr_matrix((data ,(row ,column )),shape=(dim1*dim2,dim1*dim2))
	H12 	= 	    (1.0)*kron(site1.Jz[1], site2.Jz[0])
	H12 	= H12 + (0.5)*kron(site1.Sp[1], site2.Sm[0])
	H12 	= H12 + (0.5)*kron(site1.Sm[1], site2.Sp[0])
	H12 	= H12 + (1.0)*kron(site1.H, numpy.eye(dim2))
	H12 	= H12 + (1.0)*kron(numpy.eye(dim1), site2.H)

	Jz12	= 	    (1.0)*kron(site1.Jz[0], numpy.eye(dim2))
	site12.Jz[0] = Jz12
	Jz12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Jz[1])
	site12.Jz[1] = Jz12

	Sp12	= 	    (1.0)*kron(site1.Sp[0], numpy.eye(dim2))
	site12.Sp[0] = Sp12
	Sp12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Sp[1])
	site12.Sp[1] = Sp12

	Sm12	= 	    (1.0)*kron(site1.Sm[0], numpy.eye(dim2))
	site12.Sm[0] = Sm12
	Sm12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Sm[1])
	site12.Sm[1] = Sm12

	site12.H = H12

if __name__=="__main__":

	sites = []
	site12 = Site()
	nsites = 2
	for i in range(nsites):
		sites.append(Site())

	form_H12(sites[0], sites[1], site12)

	print site12.H.toarray()
