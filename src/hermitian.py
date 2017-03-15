#!/usr/bin/env python

import numpy

def is_Hermitian(matrix):
	return numpy.allclose(matrix.T, matrix)
