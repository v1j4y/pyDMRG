#!/usr/bin/env python

import numpy
from sites import Site, SiteHubb
from scipy.sparse import csr_matrix, kron
from hermitian import is_Hermitian

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
	Hubb = Sum_i,j[ t_ij Sum_sigma (c^dag_sigma c_sigma + c^dag_sigma c_sigma) ] + Sum_i [ U (n_i,a n_i,b) ]
	'''

	dim1 = site1.H.shape[0]
	dim2 = site2.H.shape[0]

	row 	= numpy.array([0])
	column 	= numpy.array([0])
	data 	= numpy.array([0.])
	H12 	= csr_matrix((data ,(row ,column )),shape=(dim1*dim2,dim1*dim2))
	H12 	= 	    (1.0)*kron(site1.Ca[1], site2.Da[0])
	H12 	= H12 + (1.0)*kron(site1.Da[1], site2.Ca[0])
	H12 	= H12 + (1.0)*kron(site1.Cb[1], site2.Db[0])
	H12 	= H12 + (1.0)*kron(site1.Db[1], site2.Cb[0])
	H12 	= H12 + (2.0)*kron(site1.Na[1], site2.Nb[0])
	H12 	= H12 + (1.0)*kron(site1.H, numpy.eye(dim2))
	H12 	= H12 + (1.0)*kron(numpy.eye(dim1), site2.H)

	Na12	= 	    (1.0)*kron(site1.Na[0], numpy.eye(dim2))
	site12.Na[0] = Na12
	Na12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Na[1])
	site12.Na[1] = Na12

	Nb12	= 	    (1.0)*kron(site1.Nb[0], numpy.eye(dim2))
	site12.Nb[0] = Nb12
	Nb12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Nb[1])
	site12.Nb[1] = Nb12

	Ca12	= 	    (1.0)*kron(site1.Ca[0], numpy.eye(dim2))
	site12.Ca[0] = Ca12
	Ca12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Ca[1])
	site12.Ca[1] = Ca12

	Cb12	= 	    (1.0)*kron(site1.Cb[0], numpy.eye(dim2))
	site12.Cb[0] = Cb12
	Cb12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Cb[1])
	site12.Cb[1] = Cb12

	Da12	= 	    (1.0)*kron(site1.Da[0], numpy.eye(dim2))
	site12.Da[0] = Da12
	Da12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Da[1])
	site12.Da[1] = Da12

	Db12	= 	    (1.0)*kron(site1.Db[0], numpy.eye(dim2))
	site12.Db[0] = Db12
	Db12	= 	    (1.0)*kron(numpy.eye(dim1), site2.Db[1])
	site12.Db[1] = Db12

	site12.H = H12

if __name__=="__main__":

	sites = []
	site12 = SiteHubb()
	nsites = 2
	for i in range(nsites):
		sites.append(SiteHubb())

	form_H12_Hubb(sites[0], sites[1], site12)

	print (site12.Na[0].toarray()  + site12.Na[1].toarray() )
	print site12.H.toarray()
