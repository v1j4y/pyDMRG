#!/usr/bin/env python

import numpy
import copy
from scipy.sparse import csr_matrix, kron, linalg
from sites import Site, SiteHubb
from form_H12 import form_H12
from form_dmat import form_dmat

def is_Hermitian(matrix):
	return numpy.allclose(matrix.T, matrix)

def pprint(matrix):
	print numpy.array_str(matrix, precision=2, suppress_small=True)

def blocking(nsites):
	
	site0 	= Site()
	site1 	= Site()
	site2 	= Site()

	for i in range(1,nsites):

		dim1 	= site1.H.shape[0]
		dim0 	= site0.H.shape[0]
		dim2 	= site2.H.shape[0]

		dim12	= dim1*dim0
		dim34	= dim2*dim0
		dim1234	= dim12*dim34

		m		= 32
		dimO 	= m

		'''
		Forming block H12
		'''

		site12 = Site()
		form_H12(site1, site0, site12)
		dim12 = site12.H.shape[0]

		'''
		Forming block H34
		'''

		site34 = Site()
		form_H12(site2, site0, site34)

		'''
		Forming block H1234
		'''

		site1234 = Site()
		form_H12(site12, site34, site1234)
		dim1234 = site1234.H.shape[0]

		'''
		diagonalize H1234
		'''

		neig = 2

		evals, evec = linalg.eigsh(site1234.H.toarray(), k=neig, which="SA", tol=1e-09, maxiter=100000)
#	evals, evecs = numpy.linalg.eigh(site1234.H.toarray())
#	evec = evecs[:,0]
#	print(site1234.H.toarray())

		'''
		calculate dmat
		'''
		
		dmat, matO, discardw = form_dmat(evec[:,0], dim12, dim12, m)


		matOT 	= matO.transpose()
		site12.H = matOT.dot(site12.H.dot(matO))
		site12.Jz[0] = matOT.dot(site12.Jz[0].dot(matO))
		site12.Sp[0] = matOT.dot(site12.Sp[0].dot(matO))
		site12.Sm[0] = matOT.dot(site12.Sm[0].dot(matO))
		site12.Jz[1] = matOT.dot(site12.Jz[1].dot(matO))
		site12.Sp[1] = matOT.dot(site12.Sp[1].dot(matO))
		site12.Sm[1] = matOT.dot(site12.Sm[1].dot(matO))
	
		site1 = site12
		site2 = site12
		sitefilename = "site_"+str(i+1)
		numpy.save(sitefilename, site1)
#	sitefilename = "sz0"
#	numpy.savetxt(sitefilename, sztot)

def blockingHubb(nsites):
	
	site0 	= SiteHubb()
	site1 	= SiteHubb()
	site2 	= SiteHubb()

	for i in range(1,nsites):

		dim1 	= site1.H.shape[0]
		dim0 	= site0.H.shape[0]
		dim2 	= site2.H.shape[0]

		dim12	= dim1*dim0
		dim34	= dim2*dim0
		dim1234	= dim12*dim34

		m		= 32
		dimO 	= m

		'''
		Forming block H12
		'''

		site12 = SiteHubb()
		form_H12(site1, site0, site12)
		dim12 = site12.H.shape[0]

		'''
		Forming block H34
		'''

		site34 = SiteHubb()
		form_H12(site2, site0, site34)

		'''
		Forming block H1234
		'''

		site1234 = SiteHubb()
		form_H12(site12, site34, site1234)
		dim1234 = site1234.H.shape[0]

		'''
		diagonalize H1234
		'''

		neig = 2

		evals, evec = linalg.eigsh(site1234.H.toarray(), k=neig, which="SA", tol=1e-09, maxiter=100000)
#	evals, evecs = numpy.linalg.eigh(site1234.H.toarray())
#	evec = evecs[:,0]
#	print(site1234.H.toarray())

		'''
		calculate dmat
		'''
		
		dmat, matO, discardw = form_dmat(evec[:,0], dim12, dim12, m)


		matOT 	= matO.transpose()
		site12.H = matOT.dot(site12.H.dot(matO))
		site12.Jz[0] = matOT.dot(site12.Jz[0].dot(matO))
		site12.Sp[0] = matOT.dot(site12.Sp[0].dot(matO))
		site12.Sm[0] = matOT.dot(site12.Sm[0].dot(matO))
		site12.Jz[1] = matOT.dot(site12.Jz[1].dot(matO))
		site12.Sp[1] = matOT.dot(site12.Sp[1].dot(matO))
		site12.Sm[1] = matOT.dot(site12.Sm[1].dot(matO))
	
		site1 = site12
		site2 = site12
		sitefilename = "site_"+str(i+1)
		numpy.save(sitefilename, site1)
#	sitefilename = "sz0"
#	numpy.savetxt(sitefilename, sztot)

if __name__=="__main__":
	
	nsites 	= 4
	sites = []
	for i in range(0,4):
		sites.append(Site())
	'''
		mattrans = numpy.array(
							[[0/numpy.sqrt(2),	0,0,	 0/numpy.sqrt(2)],
							 [0,	1/numpy.sqrt(2), 1/numpy.sqrt(2),	0],
							 [0,	1/numpy.sqrt(2),-1/numpy.sqrt(2),	0],
							 [0/numpy.sqrt(2),	0,0,	-0/numpy.sqrt(2)]]
							)
	#mattrans  = numpy.eye(4)
		mattransT = mattrans.T
		pprint(mattransT.dot(mattrans))
	#pprint(sites[0].Sp[0].toarray())
	#	for i in range(0,3):
	#		sites[i].Jz[0] = mattransT.dot(sites[i].Jz[0].dot(mattrans))
	#		sites[i].Jz[1] = mattransT.dot(sites[i].Jz[1].dot(mattrans))
	#		sites[i].Sp[0] = mattransT.dot(sites[i].Sp[0].dot(mattrans))
	#		sites[i].Sp[1] = mattransT.dot(sites[i].Sp[1].dot(mattrans))
	#		sites[i].Sm[0] = mattransT.dot(sites[i].Sm[0].dot(mattrans))
	#		sites[i].Sm[1] = mattransT.dot(sites[i].Sm[1].dot(mattrans))
		
	#pprint(sites[0].Sp[0])
	
		H = numpy.zeros(shape=(16,16))
		site = copy.copy(sites[0])
		for i in range(0,1):
			site0 = copy.copy(site)
			dim2 = sites[i].H.shape[0]
			dim1 = site0.H.shape[0]
			site.H 	= 	   	   (1.0)*kron(site0.Jz[1], sites[i].Jz[0])
			site.H 	= site.H + (0.5)*kron(site0.Sp[1], sites[i].Sm[0])
			site.H 	= site.H + (0.5)*kron(site0.Sm[1], sites[i].Sp[0])
			site.H 	= site.H + (1.0)*kron(site0.H, numpy.eye(dim2))
			site.H 	= site.H + (1.0)*kron(numpy.eye(dim1), sites[i].H)
			Jz12	= 	    (1.0)*kron(site0.Jz[0], numpy.eye(dim2))
			site.Jz[0] = Jz12
			Jz12	= 	    (1.0)*kron(numpy.eye(dim1), sites[i].Jz[1])
			site.Jz[1] = Jz12
	
			Sp12	= 	    (1.0)*kron(site0.Sp[0], numpy.eye(dim2))
			site.Sp[0] = Sp12
			Sp12	= 	    (1.0)*kron(numpy.eye(dim1), sites[i].Sp[1])
			site.Sp[1] = Sp12
	
			Sm12	= 	    (1.0)*kron(site0.Sm[0], numpy.eye(dim2))
			site.Sm[0] = Sm12
			Sm12	= 	    (1.0)*kron(numpy.eye(dim1), sites[i].Sm[1])
			site.Sm[1] = Sm12
	
		pprint(site.H.toarray())
		site.Jz[0] = mattransT.dot(site.Jz[0].dot(mattrans))
		site.Jz[1] = mattransT.dot(site.Jz[1].dot(mattrans))
		site.Sp[0] = mattransT.dot(site.Sp[0].dot(mattrans))
		site.Sp[1] = mattransT.dot(site.Sp[1].dot(mattrans))
		site.Sm[0] = mattransT.dot(site.Sm[0].dot(mattrans))
		site.Sm[1] = mattransT.dot(site.Sm[1].dot(mattrans))
		site0 = copy.copy(site)
		dim2 = site.H.shape[0]
		dim1 = site.H.shape[0]
		site.H 	= 	   	   (1.0)*kron(site0.Jz[1], site0.Jz[0])
		site.H 	= site.H + (0.5)*kron(site0.Sp[1], site0.Sm[0])
		site.H 	= site.H + (0.5)*kron(site0.Sm[1], site0.Sp[0])
		site.H 	= site.H + (1.0)*kron(site0.H, numpy.eye(dim2))
		site.H 	= site.H + (1.0)*kron(numpy.eye(dim1), site0.H)
		pprint(site.H.toarray())
	'''
	blockingHubb(nsites)
