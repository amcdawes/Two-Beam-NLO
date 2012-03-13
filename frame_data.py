#!/usr/bin/env python
# encoding: utf-8
"""
frame_data.py

Created by Andrew M.C. Dawes on 2008-03-15.
Copyright (c) 2008 Andrew M.C. Dawes.
Some rights reserved. License: Creative Commons GNU GPL:
http://creativecommons.org/licenses/GPL/2.0/
"""

import sys
from optparse import OptionParser
from math import pi, sin, cos

VERSION = "0.4.0"

CHANGES = """
0.1.0 added pinhole function to determine start and stop coordinates based on center and width of pinhole
0.2.0 added K_to_fft function, this is useful for other things too
0.3.0 added spots function to generalize the code
0.3.1 changed the azimuth by pi so I don't measure the switch beam.
0.4.0 changed so the azimuth is correct for 0 = far right.
"""

def K_to_fft(K, N, xstep):
	""" convert coordinate in k to coordinate in fft output grid """
	return K*N*xstep/(2*pi)
	
def spots(K, N, xstep, on_azimuth, off_azimuth):
 	ring_radius = K_to_fft(K, N, xstep)
	
	#on_azimuth = 0
	#off_azimuth = pi/6
	
	onspot = [int(ring_radius*cos(on_azimuth) + N/2), int(ring_radius*sin(on_azimuth) + N/2)]
	offspot = [int(ring_radius*cos(off_azimuth) + N/2), int(ring_radius*sin(off_azimuth) + N/2)]
	return [onspot, offspot]

def pinhole(coordinates, D):
	""" calculate the row and column start and stop points for a square pinhole D by D centered at grid index [i, j] """
	if len(coordinates)==2:
		x = coordinates[0]
		y = coordinates[1]
		row_start = y - D/2
		row_stop = y + D/2
		col_start = x - D/2
		col_stop = x + D/2
		#print row_start, row_stop, col_start, col_stop
		return [row_start, row_stop, col_start, col_stop]
	else:
		return 0

def main(argv=None):

	global VERSION
	if argv is None: argv = sys.argv
	parser = OptionParser(usage="usage: %prog [options] ...",
						  version=VERSION)
	parser.set_defaults(verbose=False,gridsize=256,K=14.49)
	parser.add_option("-o", "--output",
					  action="store",
					  type="string",
					  help="send output to FILE",
					  metavar="FILE",
					  dest="outputfile")
	parser.add_option("-v", "--verbose",
					  action="store_true",
					  help="generate verbose output",
					  dest="verbose")
	parser.add_option("-k", action="store", type="float", help="K vector of pattern", dest="K")
	parser.add_option("-n", action="store", type="int", help="grid size N x N", dest="gridsize") 

	(options, args) = parser.parse_args()
	
		
	#[onspot, offspot] = spots(14.49, 256, 4.0/256, 0, pi/6)
	[onspot_coord, offspot_coord] = spots(options.K, options.gridsize, 4.0/options.gridsize, pi, pi/2)
	
	[on_row_start, on_row_stop, on_col_start, on_col_stop] = pinhole(onspot_coord, 6)
	[off_row_start, off_row_stop, off_col_start, off_col_stop] = pinhole(offspot_coord, 6)

	outputfile = open('response.dat','w')
	# main loop:
	for file in args:
		# image = read_array(file)
		data = open(file)
		
		row = 0
		offspot = 0
		onspot = 0
		for line in data:
			# take the first index between off_row_start and off_row_stop
			if row >= off_row_start and row < off_row_stop:
				offspot += sum(map(float,line.split()[off_col_start:off_col_stop]))
				
			# take the first index between on_row_start and on_row_stop
			if row >= on_row_start and row < on_row_stop:
				onspot += sum(map(float,line.split()[on_col_start:on_col_stop]))
			row += 1
		data.close()

		
		frame = file.split("_")[0]
		
		outputfile.write("%s\t%10f\t%10f\n" % (frame, onspot, offspot))
		
		print "%s\t%10.8f\t%10.8f\n" % (frame, onspot, offspot)
		
	print "onspot pinhole ", pinhole(onspot_coord, 6), "\n"
	print "offspot pinhole ", pinhole(offspot_coord, 6), "\n"
	print "onspot coord ", onspot_coord, "\n"
	print "offspot coord ", offspot_coord, "\n"


if __name__ == "__main__":
    sys.exit(main())
