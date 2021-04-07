#!/usr/bin/env python
# -*- coding: utf-8 -*-


"""
This is part of the TIMEleSS tools
http://timeless.texture.rocks/

Copyright (C) M. Krug, Münster Univ. Germany
Copyright (C) S. Merkel, Univ. Lille France

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

# Python 2 to python 3 migration tools
from __future__ import absolute_import
from __future__ import print_function

# System functions, to manipulate command line arguments
import sys
import argparse
import os.path
from argparse import RawTextHelpFormatter

# Maths stuff
import matplotlib
import matplotlib.pyplot as plt
matplotlib.use("Qt5Agg")
import numpy

# TIMEleSS parsing utilities
from TIMEleSS.general import multigrainOutputParser
from TIMEleSS.general import indexedPeak3DXRD


def absolute_grainsizes(grainsizelist, beamsize_H, beamsize_V, rotationangle, gasketthickness, indexquality, volume):
    with open(grainsizelist) as g:
        grainsizes = g.readlines()
    total = 0
    for grain in grainsizes:
        grain = float(grain)
        if volume == False:
            grain = 4/3*numpy.pi()*grain^3 # Turn grain radii into grain volumes
        total += grain
    total = total * indexquality / 100 # Account for the indexing quality
    
    # Calculate the sample chamber volume. For more info on the formula ask M. Krug.
    samplechambervol = gasketthickness * beamsize_H * beamsize_V * numpy.arccos(rotationangle*numpy.pi/180/2)
    
    ratio_V = total / samplechambervol # How many µm^3 equals one relative grain size unit
    ratio_V = float(ratio_V)
    print (ratio_V)
    ratio_R = ratio_V**(1/3)
    
    # Create a new file that contains the absolute grain size 
    newfile = grainsizelist[:-4] + "_abs.txt"
    grainsizes_new = [] # Make a list of the new grain sizes
    string = ""
    for grain in grainsizes:
        grain = float(grain)
        if volume == False:
            grain = grain * ratio_R
        else:
            grain = grain * ratio_V
        grainsizes_new.append(grain)
        string += str(grain) + "\n"
    f= open(newfile,"w+")
    f.write(string)
    f.close()
    
    if volume == True:
        print ("\nA volumetric grain size of 1.0 in your list of grainsizes corresponds to %0.3f µm^3." % (ratio_V))
        print ("\nSaved new list of grain volumes (in µm^3) as %s." % (newfile))
    else:
        print ("\nA grain radius of 1.0 in your list of grainsizes corresponds to %0.3f µm." % (ratio_R))
        print ("\nSaved new list of grain radii (in µm) as %s." % (newfile))
    return grainsizes_new



################################################################
#
# Main subroutines
#
#################################################################

class MyParser(argparse.ArgumentParser):
    """
    Extend the regular argument parser to show the full help in case of error
    """
    def error(self, message):
        
        sys.stderr.write('\nError : %s\n\n' % message)
        self.print_help()
        sys.exit(2)


def main(argv):
    """
    Main subroutine
    """
    
    parser = MyParser(usage='%(prog)s  [options] -w wavelength logfile.log fltfile.flt ciffile.cif', description="""Estimation of relative grain volumes based on diffraction intensities

Example:
    %(prog)s -w 0.3738 mylogfile.log peaks-t100.flt crystal.cif

This is part of the TIMEleSS project\nhttp://timeless.texture.rocks
""", formatter_class=RawTextHelpFormatter)
    
    # Required arguments
    parser.add_argument('grainsizelist', help="Path and name for the GrainSize output file")
    parser.add_argument('-H', '--beamsize_H', required=True, help="Beamsize perpendicular to rotation axis (µm), usually horizontal (required)", type=float)
    parser.add_argument('-V', '--beamsize_V', required=True, help="Beamsize parallel to rotation axis (µm), usually vertical (required)", type=float)
    parser.add_argument('-r', '--rotationangle', required=True, help="Full rotation angle. Example: [-28,+28] rotation = 56 degrees (required)", type=float)
    parser.add_argument('-t', '--gasketthickness', required=True, help="Thickness of your gasket indentation (µm, required)", type=float)
    parser.add_argument('-i', '--indexquality', required=True, help="Percentage of indexed g-vectors (in %, required)", type=float)
    
    # Optionnal arguments
    parser.add_argument('-vol', '--volume', required=False, help="If True, treats grainsizelist as list of grain volumes. If False, treats grainsizelist as list of grain radii. Default is %(default)s", default=True, type=bool)
    parser.add_argument('-hist', '--histogram_bins', required=False, help="If a histogram shall be plotted, give the number of histogram bins here. Default is %(default)s", default=None, type=int)
    
    # Parse arguments
    args = vars(parser.parse_args())
    grainsizelist = args['grainsizelist']
    beamsize_H = args['beamsize_H']
    beamsize_V = args['beamsize_V']
    rotationangle = args['rotationangle']
    gasketthickness = args['gasketthickness']
    indexquality = args['indexquality']
    volume = args['volume']
    histogram_bins = args['histogram_bins']

    grainsizes_new = absolute_grainsizes(grainsizelist, beamsize_H, beamsize_V, rotationangle, gasketthickness, indexquality, volume)

    # Make a histogram
    if histogram_bins != None:
        if volume == True:
            print ("Plotting histogram ...\n")
            plt.hist(grainsizes_new, bins = histogram_bins)
            plt.xlabel("Grain volume ($\mu$m^3)")
            plt.ylabel("Number of grains")
            plt.title("n = %s" % len(grainsizes_new), fontsize = 20)
            plt.show()
        else:
            radii = []
            for item in grainsizes_new:
                radius = (3*item/4/numpy.pi())^(1/3)
                radii.append(radius)
            print ("Plotting histogram ...\n")
            plt.hist(radii, bins = histogram_bins)
            plt.xlabel("Grain radii ($\mu$m)")
            plt.ylabel("Number of grains")
            plt.title("n = %s" % len(grainsizes_new), fontsize = 20)
            plt.show()
            
        
# Calling method 1 (used when generating a binary in setup.py)
def run():
    main(sys.argv[1:])

# Calling method 2 (if run from the command line)
if __name__ == "__main__":
    main(sys.argv[1:])
