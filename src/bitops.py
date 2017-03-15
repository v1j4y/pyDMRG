#!/usr/bin/env python

def popcnt(n):
	res = 0
	value = n
	while(value):
		value &= (value-1)
		res += 1
	return res

def int2str(n):
	return "{0:b}".format(n)

if __name__ == "__main__":
	for i in range(10):
		print popcnt(i), int2str(i)
