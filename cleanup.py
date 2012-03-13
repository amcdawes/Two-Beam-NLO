#!/usr/bin/env python
# encoding: utf-8
"""
cleanup.py

A script for converting from png to jpg files, and then making a movie. Cleans up after itself also.

Copyright (c) 2007 Andrew M.C. Dawes.
Some rights reserved. License: Creative Commons GNU GPL:
http://creativecommons.org/licenses/GPL/2.0/
"""

VERSION = 0.6

import sys
import os
import getopt
from subprocess import *

def convert(workingdir):
	convertcode = Popen("mogrify -format jpg -quality 90 *.png",shell=True).wait()
	print ("\n\n\tdone mogrify, returned: %i\n" % convertcode)

	ffmpegcode = Popen("ffmpeg -i %03d_farfield.jpg " + workingdir + "_farfield.mp4",shell=True).wait()
	print ("\n\n\tdone ffmpeg round 1, returned: %i\n" % ffmpegcode)

	framempegcode = Popen("ffmpeg -i %03d_frame.jpg " + workingdir + "_frames.mp4",shell=True).wait()
	print ("\n\n\tdone ffmpeg round 2, returned: %i\n" % framempegcode)
	
	if (ffmpegcode == 0 and framempegcode == 0):	
		deletecode = Popen("rm -f *.jpg *.png",shell=True).wait()
		print ("\n\n\tdone rm, returned: %i\n" % deletecode)

def main(argv=None):
	workingdir = os.path.basename(os.getcwd())
	convert(workingdir)

if __name__ == "__main__":
	sys.exit(main())