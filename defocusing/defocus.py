#!/usr/bin/env python
# encoding: utf-8
"""
switching.py

Created by  on 2007-12-06.
Copyright (c) 2007 Andrew M.C. Dawes.
Some rights reserved. License: Creative Commons GNU GPL:
http://creativecommons.org/licenses/GPL/2.0/
"""

VERSION = 0.6

import sys
import os
import getopt
from subprocess import *


help_message = '''
Execute the command ./switching with parameters specified in a script file of the following format:

workdir	sigma Tr IL NZ pumpx pumpy kin bkinx bkiny noiseamp	seed probeamp probewidth switchon

i.e.:

newdirA 400 90 0.565 20 0.0 0.0 49 0.0 0.0 0.0000001 12345678 0.000002 0.2 60
newdirB 400 90 0.565 20 0.0 0.0 49 0.0 0.0 0.0000001 12345678 0.000002 0.2 60
'''


class Usage(Exception):
	def __init__(self, msg):
		self.msg = msg


def simulate(scriptfilename):
	startingpath = os.getcwd() # save the current path
	scriptfile = file(scriptfilename,'r')

	print 'Using %s as script file' % scriptfilename

	for line in scriptfile:
		entry = line.split()
		if len(entry)!=15:
			print "Script line doesn't have 15 entries"
			break

		workingdir = entry[0]
		if not os.path.exists("./" + workingdir + "/"):
			os.mkdir(workingdir) # create the directory specified in the script
			os.chdir(workingdir) # change to this dir
			
			logfile = file("log.dat", 'w')
			logfile.write(line)
			logfile.close()
			
			print "Running Simulation in: ", workingdir
			# trapping stderr to a PIPE keeps it from being shown on the terminal
			# to do this, add stderr=PIPE as well:
			p = Popen(["~/code/cpp/nonlinearoptics/defocusing/defocus"], shell=True, stdin=PIPE, stdout=PIPE, close_fds=True)
			stdout_text, stderr_text = p.communicate(
													entry[1]+"\n"
													+entry[2]+"\n"
													+entry[3]+"\n"
													+entry[4]+"\n"
													+entry[5]+"\n"
													+entry[6]+"\n"
													+entry[7]+"\n"
													+entry[8]+"\n"
													+entry[9]+"\n"
													+entry[10]+"\n"
													+entry[11]+"\n"
													+entry[12]+"\n"
													+entry[13]+"\n"
													+entry[14]+"\n")
			
			# as a log, print the stdout:
			print stdout_text
													
			print "cleaning up in " + workingdir
			
			retcode1 = Popen("mogrify -format jpg -quality 90 *.png",shell=True).wait()
			print ("\n\n\tdone mogrify, returned: %i\n" % retcode1)
			
			retcode2 = Popen("ffmpeg -i %03d_farfield.jpg " + workingdir + "_farfield.mp4",shell=True).wait()
			print ("\n\n\tdone ffmpeg round 1, returned: %i\n" % retcode2)

			retcode3 = Popen("ffmpeg -i %03d_frame.jpg " + workingdir + "_frames.mp4",shell=True).wait()
			print ("\n\n\tdone ffmpeg round 2, returned: %i\n" % retcode3)
				
			if (retcode3 == 0 and retcode2 == 0):	
				retcode4 = Popen("rm -f *.jpg *.png",shell=True).wait()
				print ("\n\n\tdone rm, returned: %i\n" % retcode4)
			
			
		else:
			print "Directory " + workingdir + " already exists, skipping to next script entry."
		
		os.chdir(startingpath) # set back to original working dir


def main(argv=None):
	if argv is None:
		argv = sys.argv
	try:
		try:
			opts, args = getopt.getopt(argv[1:], "ho:v", ["help", "output="])
		except getopt.error, msg:
			raise Usage(msg)
	
		# option processing
		for option, value in opts:
			if option == "-v":
				verbose = True
			if option in ("-h", "--help"):
				raise Usage(help_message)
			if option in ("-o", "--output"):
				output = value
				
		for scriptfile in args:
			simulate(scriptfile)
	
	except Usage, err:
		print >> sys.stderr, sys.argv[0].split("/")[-1] + ": " + str(err.msg)
		print >> sys.stderr, "\t for help use --help"
		return 2

if __name__ == "__main__":
	sys.exit(main())
