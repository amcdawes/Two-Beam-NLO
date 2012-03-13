#!/usr/bin/env python
# encoding: utf-8
"""
traverse.py

Created by Andrew M. C. Dawes on 2009-03-19.
Copyright (c) 2009 Pacific University. All rights reserved.
"""

import sys, os
import getopt
from subprocess import Popen


help_message = '''
traverses a set of directories and runs the same command in each. It assumes 
names like runA, runB, runC etc. The command is used as:

traverse.py 

or

traverse.py -n 3 (to indicate 3 runs, or runA, runB, and runC)

So far, the command is hardwired to be: frame_data.py -k 15.49 -n 50 *_ff_int.dat
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def main(argv=None):
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ho:vn:", ["help", "output=","numruns="])
		except getopt.error, msg:
			raise Usage(msg)
		
		number_of_runs = 10 # note, this an input so it may be overridden
		
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-o", "--output"):
				output = value
			if option in ("-n", "--numruns"):
				number_of_runs = int(value)
		
		
		
		# create a list containing the run directory names ["runA", "runB", ...]
		tree = ["run"+letter for letter in map(chr, range(65,65+number_of_runs))]
		
		command = "frame_data.py -k 15.49 -n 50 *_ff_int.dat"
		
		startingdir = os.getcwd()
		
		for folder in tree:
			os.chdir(startingdir+"/"+folder)
			currentdir = os.getcwd()
			print "running ", command, " in ", currentdir
			framedata = Popen(command,shell=True).wait()
			print ("\n\n\tdone with "+command+", returned: %i\n" % framedata)
	
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2


if __name__ == "__main__":
	sys.exit(main())
