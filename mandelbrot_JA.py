#!/usr/bin/python

from math import sqrt
from time import sleep

grid = [[((float(x)/120 + 1j*(float(y)/60))) for x in range(-160,171)] for y in range(120,-81,-1)]

def magnitude(n):
	return sqrt(n.real**2 + n.imag**2)

def m_iterate(grid1, cycles):
	new_grid = []
	for row in grid1:
		new_row = []
		for number in row:
			z = 0
			for cycle in range(cycles):
				new_number = z**2 + number
				if magnitude(new_number) < 3.0:
					z = new_number
				else:
					z = (3.0 + 0.0j)
			new_row.append(z)
		new_grid.append(new_row)
	return new_grid

for time in range(45):
	for row in m_iterate(grid,time):
		print "".join(["X" if magnitude(i) <= (2.0) else " " for i in row])
	sleep(0.25)
