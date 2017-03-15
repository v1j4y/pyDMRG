#!/usr/bin/env python

from blocking import blocking
from sweep import sweep


if __name__ == "__main__":

	nsites = 6

	'''
	do blocking first
	'''

	blocking(nsites)
	print "done blocking"

	'''
	next do the weep iterations
	'''

	sweep(nsites)
	print "done sweep"
