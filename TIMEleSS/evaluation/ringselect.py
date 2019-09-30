#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This is part of the TIMEleSS tools
http://timeless.texture.rocks/
Copyright (C) S. Merkel, Universite de Lille, France
This program is free software; you can redistribute it and/or
modify it under the terms of the GNU General Public License
as published by the Free Software Foundation; either version 2
of the License, or (at your option) any later version.
This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.
You should have received a copy of the GNU General Public License
along with this program; if not, write to the Free Software
Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA  02110-1301, USA.
"""

# This script is based on ringselect.py from the ImageD11 project (J. Wright) and was modified for the purposes of the TIMEleSS project.

# System functions, to manipulate command line arguments
import sys
import argparse
import os.path
#from __future__ import print_function

from ImageD11.columnfile import columnfile
from ImageD11.unitcell import unitcell_from_parameters, unitcell
from ImageD11.parameters import read_par_file
from ImageD11.transform import compute_tth_eta
import polyxsim.reflections as rf

# Maths stuff
import numpy as np
import pylab as pl

# TIMEleSS parsing utilities
from TIMEleSS.general import multigrainOutputParser

# Will use to crystallography functions in xfab.symmetry
#import xfab.symmetry
from xfab import tools,structure,sg,symmetry

#######################################################################################################################

# This function is copied from ImageD11.unitcell and then modified.
def makerings(limit,tol=0.001):
	"""
	Makes a list of computed powder rings
	The tolerance is the difference in d* to decide
	if two peaks overlap
	"""
	"""
	peaks=unitcell.gethkls(self,limit+tol) # [ ds, [hkl] ]
	
	# reflections.py (rf) needs some param file which we don't have so we create a dummy dictionary which serves as it.
	param = {'theta_min': 0.0,
		  'theta_max': 9.0, 
		  'wavelength': 0.2894, 
		  'unit_cell_phase_0': [3.32,3.32,3.32,90.00,90.00,90.00], 
		  'sgname_phase_0': 'pm-3m', 
		  'cell_choice_phase_0': 'standard', 
		  'structure_phase_0': 'KCl-B2_20GPa.cif', 
		  'sgno_phase_0': 221, 
		  'structure_int': 5}
	hkl = rf.gen_miller(param, 0)
	struct = rf.open_structure(param, 0)
	hkl_intensity = rf.calc_intensity(hkl, struct)
	
	# If hkl and peaks are equal: leave them. Else: Delete them from the peaks list.
	for item in hkl_intensity:
		if item[3] = 0:
			peaks.remove()
	
	
	ringds=[]   # a list of floats
	ringhkls={} # a dict of lists of integer hkl
	# Append first peak
	peak = peaks[0]
	ringds.append(peak[0])
	ringhkls[peak[0]] = [peak[1]]
	for peak in peaks[1:]:
		if abs(peak[0] - ringds[-1]) < tol:
			ringhkls[ringds[-1]].append(peak[1])
		else:
			ringds.append(peak[0])
			ringhkls[ringds[-1]]= [peak[1]]
	ringtol = tol
"""
#################################################################################################

# Get the peaks of a certain phase
def ringselect(inputfile, par, tol, output_phase_only, tthrange, pixels, top, phase):
	c = columnfile(inputfile)
	p = read_par_file(par)
	c.filter( c.Number_of_pixels > pixels )
	u = unitcell_from_parameters( p )
	w = p.get("wavelength")

	tth, eta = compute_tth_eta( (c.sc, c.fc), **p.parameters)

	dsmax = 2*np.sin(1.03*tth.max()*np.pi/360)/w
	
	u.makerings(dsmax)		# Actual program
	#makerings(dsmax)		# Attempt to insert the cif file
	
	#print u.ringds
	#print u.ringhkls
	
	mask = np.zeros( c.nrows, dtype=np.bool )
	for d in u.ringds:
		tthc = np.arcsin(w*d/2)*360/np.pi
		M = len(u.ringhkls[d])
		#print ("hkl",u.ringhkls[d][-1],M, end=" " )
		sel =  abs( tth - tthc ) < tol
		nring = sel.sum()
		#print (nring, end=" " )
		if nring == 0:
			print()
			continue
		intensities = c.IMax_int[sel]
		idxsel = np.argsort( intensities )
		npkstokeep = top * M
		idxkeep = np.arange(c.nrows)[sel] 
		#print( npkstokeep, len(idxkeep), end=" ")
		if npkstokeep < nring:
			idxs = np.take( idxkeep, idxsel[-npkstokeep:] )
		else:
			idxs = idxkeep 
		#print( len(idxkeep),len(idxsel[-npkstokeep:] ), end=" ")
		np.put( mask, idxs, True)
		#print( mask.sum() )
	c.filter( mask )
	pl.show()
	if output_phase_only == "inputfile_phase_only.flt":
		alt_name = "%s_%s_only.flt" % (inputfile, phase)
		c.writefile(alt_name)
		print "\nPeaks belonging to %s were extracted from >%s< and saved in >%s<." % (phase, inputfile, alt_name)
	else:
		c.writefile(output_phase_only)
		print "\nPeaks belonging to %s were extracted from >%s< and saved in >%s<." % (phase, inputfile, output_phase_only)

##########################################################################################################

# Get the leftover peaks
def ringselect_reverse(inputfile, output_phase_removed, output_phase_only, phase):
	alt_name = "%s_%s_only.flt" % (inputfile, phase)
	with open(inputfile, 'r') as file1:
		if output_phase_only == "inputfile_phase_only.flt":
			with open(alt_name, 'r') as file2:
				same = set(file1).difference(file2)
		else:
			with open(output_phase_only, 'r') as file2:
				same = set(file1).difference(file2)

	same.discard('\n')

	if output_phase_removed == "inputfile_phase_removed.flt":
		alt_name2 = "%s_%s_removed.flt" % (inputfile, phase)
		with open(alt_name2, 'w') as file_out:
			file_out.write("#  sc  fc  omega  Number_of_pixels  avg_intensity  s_raw  f_raw  sigs  sigf  covsf  sigo  covso  covfo  sum_intensity  sum_intensity^2  IMax_int  IMax_s  IMax_f  IMax_o  Min_s  Max_s  Min_f  Max_f  Min_o  Max_o  dety  detz  onfirst  onlast  spot3d_id\n")
			for line in same:
				file_out.write(line)
		print "Afterwards, the peak list >%s< was subtracted by the peaks from >%s<. Result was saved in >%s<.\n" % (inputfile, alt_name, alt_name2)
	else:
		with open(output_phase_removed, 'w') as file_out:
			file_out.write("#  sc  fc  omega  Number_of_pixels  avg_intensity  s_raw  f_raw  sigs  sigf  covsf  sigo  covso  covfo  sum_intensity  sum_intensity^2  IMax_int  IMax_s  IMax_f  IMax_o  Min_s  Max_s  Min_f  Max_f  Min_o  Max_o  dety  detz  onfirst  onlast  spot3d_id\n")
			for line in same:
				file_out.write(line)
		print "Afterwards, the peak list >%s< was subtracted by the peaks from >%s<. Result was saved in >%s<.\n" % (inputfile, output_phase_only, output_phase_removed)
	
##########################################################################################################################

#Main subroutines
class MyParser(argparse.ArgumentParser):
	"""
	Extend the regular argument parser to show the full help in case of error
	"""
	def error(self, message):
		
		sys.stderr.write('\nError : %s\n\n' % message)
		self.print_help()
		sys.exit(2)

def main(argv):
	
	parser = MyParser(usage='%(prog)s -i inputfile -p parameterfile -t tolerance [options]', description="Takes a list of peaks and filters out those which belong to a certain phase. As a result, two files are created: One containing the extracted peaks only and one containing the leftover peaks.\nThis is part of the TIMEleSS project\nhttp://timeless.texture.rocks\n")
	
	# Required arguments
	parser.add_argument('-i', '--inputfile', required=True, help="Input peaks file (flt) from PeakSearch (required).")
	parser.add_argument('-p', '--parameter_file', required=True, help="Parameter file from ImageD11 (required).")
	parser.add_argument('-t', '--twotheta_tolerance', required=True, help="Tolerance in 2theta. (Example: Tolerance of 0.15 means: 2theta +- 0.15 degrees) (required).", type=float)
	parser.add_argument('-phase', '--phase', required=True, help="Comment on which part of the data (e.g. which phase) is removed. Default is %(default)s", default="phase", type=str)
	
	# Optional arguments
	parser.add_argument('-r', '--twotheta_range', required=False, help="Maximum 2theta angle to be considered. Default is %(default)s", default=20.0, type=float)
	parser.add_argument('-s', '--peaksize_min', required=False, help="Minimum peak size (in pixels) to be considered. Default is %(default)s", default=0, type=int)
	parser.add_argument('-n', '--peaknumber_max', required=False, help="Maximum number of peaks per ring (put a really large number if you want to be sure to have them all). Default is %(default)s", default=100000, type=int)
	parser.add_argument('-only', '--output_phase_only', required=False, help="Output peaks file (flt) of the removed peaks (Please add the ending .flt). Default is %(default)s", default="inputfile_phase_only.flt", type=str)
	parser.add_argument('-removed', '--output_phase_removed', required=False, help="Output peaks file (flt) of the leftover after removing some peaks (Please add the ending .flt). Default is %(default)s", default="inputfile_phase_removed.flt", type=str)
	
	args = vars(parser.parse_args())
	
	inputfile = args['inputfile']
	par = args['parameter_file']
	tol = args['twotheta_tolerance']
	output_phase_only = args['output_phase_only']
	output_phase_removed = args['output_phase_removed']
	tthrange = args['twotheta_range']/2
	pixels = args['peaksize_min']
	top = args['peaknumber_max']
	phase = args['phase']
	
	ringselect(inputfile, par, tol, output_phase_only, tthrange, pixels, top, phase);
	ringselect_reverse(inputfile, output_phase_removed, output_phase_only, phase);
	print "Ringselect created two files:\n1. %s_%s_only.flt, which contains only the peaks that belong to %s.\n2. %s_%s_removed.flt, which contains the leftover (the peaks that don't belong to to %s)." % (inputfile, phase, phase, inputfile, phase, phase)
		
##########################################################################################################

# Calling method 1 (used when generating a binary in setup.py)
def run():
	main(sys.argv[1:])

# Calling method 2 (if run from the command line)
if __name__ == "__main__":
	main(sys.argv[1:])
