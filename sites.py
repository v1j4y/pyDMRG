#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix

class Site(object):
	'''Class that contains the MPO's for each site
	
	Attributes:
		Jz:	The Jz matrix
		Sp: The step up matrix
		Sm:	The step down matrix
		H:	The Hamiltonian which is the 0 matrix
	'''
		

	def __init__(self):
		rowindicesJz 	= numpy.array([0,1])
		columnindicesJz = numpy.array([0,1])
		dataJz 			= numpy.array([1./2.,-1./2.])
		self.Jz = csr_matrix((dataJz,(rowindicesJz,columnindicesJz)),shape=(2,2))

		rowindicesSp 	= numpy.array([0])
		columnindicesSp = numpy.array([1])
		dataSp 			= numpy.array([1.])
		self.Sp = csr_matrix((dataSp,(rowindicesSp,columnindicesSp)),shape=(2,2))

		rowindicesSm 	= numpy.array([1])
		columnindicesSm = numpy.array([0])
		dataSm 			= numpy.array([1.])
		self.Sm = csr_matrix((dataSm,(rowindicesSm,columnindicesSm)),shape=(2,2))

		rowindicesH  	= numpy.array([0])
		columnindicesH  = numpy.array([0])
		dataH  			= numpy.array([0.])
		self.H  = csr_matrix((dataH ,(rowindicesH ,columnindicesH )),shape=(2,2))

if __name__ == "__main__":

	sites = []

	nsites = 1
	for i in range(nsites):
		sites.append(Site())

	print sites[0].Jz.toarray()
