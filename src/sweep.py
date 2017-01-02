#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron, linalg
from sites import Site
from form_H12 import form_H12
from form_dmat import form_dmat, pprint, is_Hermitian
from blocking import blocking

def sweep(nsites):

	'''
	now do first sweep
	'''

	site0 	= Site()

	for i in range(nsites/2, 1, -1):

			sitefilename = "site_"+str(i)+".npy"
			site2 	= numpy.load(sitefilename).item()

			if i == nsites/2:
				site1 	= site2

			dim1 	= site1.H.shape[0]
			dim0 	= site0.H.shape[0]
			dim2 	= site2.H.shape[0]

			dim12	= dim1*dim0
			dim34	= dim2*dim0
			dim1234	= dim12*dim34

			if dim12 < 54:
				m 	= dim12
			else:
				m	= 54

			dimO 	= m

			sites = []

			sites.append(site1)
			sites.append(site0)
			sites.append(site0)
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
			form_H12(sites[2], sites[3], site34)

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

			evals, evecs = linalg.eigsh(site1234.H, k=neig, which="SA", tol=1e-08, maxiter=100000)
#	eval	s, evecs = numpy.linalg.eigh(site1234.H.toarray())
			idx = evals.argsort()[::-1]   
			evals = evals[idx]
			evecs = evecs[:,idx]
			evec	 = evecs[:,1]
			print evals[1], evals[0]

			'''
			calculate dmat
			'''
		
			dmat, matO = form_dmat(evec, dim12, dim34)


			matOT 	= matO.transpose()
			site12.H = matOT.dot(site12.H.dot(matO))
			site12.Jz[0] = matOT.dot(site12.Jz[0].dot(matO))
			site12.Sp[0] = matOT.dot(site12.Sp[0].dot(matO))
			site12.Sm[0] = matOT.dot(site12.Sm[0].dot(matO))
			site12.Jz[1] = matOT.dot(site12.Jz[1].dot(matO))
			site12.Sp[1] = matOT.dot(site12.Sp[1].dot(matO))
			site12.Sm[1] = matOT.dot(site12.Sm[1].dot(matO))
			
			site1 = site12
			sitefilename = "site_"+str((nsites+1)-i)
			numpy.save(sitefilename, site1)

	'''
	now do second sweep
	'''

	site0 	= Site()

	for i in range(nsites-1, nsites/2-1, -1):

			sitefilename = "site_"+str(i)+".npy"
			site2 	= numpy.load(sitefilename).item()

			if i == nsites-1:
				site1 	= Site()

			dim1 	= site1.H.shape[0]
			dim0 	= site0.H.shape[0]
			dim2 	= site2.H.shape[0]

			dim12	= dim1*dim0
			dim34	= dim2*dim0
			dim1234	= dim12*dim34

			if dim12 < 54:
				m 	= dim12
			else:
				m	= 54

			dimO 	= m

			sites = []

			sites.append(site1)
			sites.append(site0)
			sites.append(site0)
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
			form_H12(sites[2], sites[3], site34)

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

			evals, evecs = linalg.eigsh(site1234.H, k=neig, which="SA", tol=1e-08, maxiter=100000)
#	eval	s, evecs = numpy.linalg.eigh(site1234.H.toarray())
			idx = evals.argsort()[::-1]   
			evals = evals[idx]
			evecs = evecs[:,idx]
			evec	 = evecs[:,1]
			print evals[1], evals[0]-evals[1]

			'''
			calculate dmat
			'''
			
			dmat, matO = form_dmat(evec, dim12, dim34)


			matOT 	= matO.transpose()
			site12.H = matOT.dot(site12.H.dot(matO))
			site12.Jz[0] = matOT.dot(site12.Jz[0].dot(matO))
			site12.Sp[0] = matOT.dot(site12.Sp[0].dot(matO))
			site12.Sm[0] = matOT.dot(site12.Sm[0].dot(matO))
			site12.Jz[1] = matOT.dot(site12.Jz[1].dot(matO))
			site12.Sp[1] = matOT.dot(site12.Sp[1].dot(matO))
			site12.Sm[1] = matOT.dot(site12.Sm[1].dot(matO))
			
			site1 = site12
			sitefilename = "site_"+str((nsites+1)-i)
			numpy.save(sitefilename, site1)

if __name__ == "__main__":

	nsites = 14

	'''
	do blocking first
	'''

	blocking(nsites)
	print "done blocking"

	sweep(nsites)
	print "done sweep"
