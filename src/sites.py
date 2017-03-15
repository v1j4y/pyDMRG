#!/usr/bin/env python

import numpy
from scipy.sparse import csr_matrix, kron

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
		Jz 				= csr_matrix((dataJz,(rowindicesJz,columnindicesJz)),shape=(2,2))
		self.Jz			= []
		self.Jz.append(Jz)
		self.Jz.append(Jz)

		rowindicesSp 	= numpy.array([1])
		columnindicesSp = numpy.array([0])
		dataSp 			= numpy.array([1.])
		Sp 				= csr_matrix((dataSp,(rowindicesSp,columnindicesSp)),shape=(2,2))
		self.Sp			= []
		self.Sp.append(Sp)
		self.Sp.append(Sp)

		rowindicesSm 	= numpy.array([0])
		columnindicesSm = numpy.array([1])
		dataSm 			= numpy.array([1.])
		Sm 				= csr_matrix((dataSm,(rowindicesSm,columnindicesSm)),shape=(2,2))
		self.Sm			= []
		self.Sm.append(Sm)
		self.Sm.append(Sm)

		rowindicesH  	= numpy.array([0])
		columnindicesH  = numpy.array([0])
		dataH  			= numpy.array([0.])
		self.H  = csr_matrix((dataH ,(rowindicesH ,columnindicesH )),shape=(2,2))

class SiteHubb(Site):
	'''
	Hubbard model

	Class that contains the MPO's for each site

			0 a  b  ab
	basis: [0,01,10,11]
	
	Attributes:
		Jz:	The Jz matrix
		Sp: The step up matrix
		Sm:	The step down matrix
		Ca:	The alpha creation operator
		Cb:	The beta creation operator
		Da:	The alpha destruction operator
		Db:	The beta destruction operator
		H:	The Hamiltonian which is the 0 matrix
	'''
		

	def __init__(self):
		obj = Site()
		rowindicesJz 	= numpy.array([1,2])
		columnindicesJz = numpy.array([1,2])
		dataJz 			= numpy.array([1./2.,-1./2.])
		Jz 				= csr_matrix((dataJz,(rowindicesJz,columnindicesJz)),shape=(4,4))
		self.Jz			= []
		self.Jz.append(Jz)
		self.Jz.append(Jz)

		Sp 				= kron(obj.Sm[0], obj.Sp[0])
		self.Sp			= []
		self.Sp.append(Sp)
		self.Sp.append(Sp)

		Sm 				= kron(obj.Sp[0], obj.Sm[0])
		self.Sm			= []
		self.Sm.append(Sm)
		self.Sm.append(Sm)

		Ca 				= kron(numpy.eye(2), obj.Sp[0])
		self.Ca			= []
		self.Ca.append(Ca)
		self.Ca.append(Ca)

		Cb 				= kron(obj.Sp[0], numpy.eye(2))
		self.Cb			= []
		self.Cb.append(Cb)
		self.Cb.append(Cb)

		Da 				= kron(numpy.eye(2), obj.Sm[0])
		self.Da			= []
		self.Da.append(Da)
		self.Da.append(Da)

		Db 				= kron(obj.Sm[0], numpy.eye(2))
		self.Db			= []
		self.Db.append(Db)
		self.Db.append(Db)

		Na 				= Ca.dot(Da)
		self.Na			= []
		self.Na.append(Na)
		self.Na.append(Na)

		Nb 				= Cb.dot(Db)
		self.Nb			= []
		self.Nb.append(Nb)
		self.Nb.append(Nb)

		rowindicesH  	= numpy.array([0])
		columnindicesH  = numpy.array([0])
		dataH  			= numpy.array([0.])
		self.H  = csr_matrix((dataH ,(rowindicesH ,columnindicesH )),shape=(4,4))

if __name__ == "__main__":

	sites = []

	nsites = 2
	for i in range(nsites):
		sites.append(SiteHubb())

	print sites[0].Na[0].toarray()
