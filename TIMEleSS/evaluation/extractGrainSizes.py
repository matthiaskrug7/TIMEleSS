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
# from TIMEleSS.evaluation import testGrainsPeaks  # Not necessary
from TIMEleSS.general import indexedPeak3DXRD
from TIMEleSS.general import cifTools


# This script determines an average relative grain size of a multigrain dataset.

# S. Merkel: Not sure what this is for
# g = testGrainsPeaks.grainPlotData()
# Moved histogram out of the function. I do not like to mix calculations and gui in a single function

def grainSizeEstimate(logfile, fltfile, ciffile, wavelength, ttheta_min=None, ttheta_max=None, output=None, kickoutfactor=20):
    """
    Determination of grain size statistics based on diffraction intensities
    Diffracting intensities for each peak is normalized by the theoretical intensity (structure factor and Lorentz correction) for that peak
    
    Params:
     - logfile
     - fltfile
     - ciffile
     - wavelength (in angstroms)
     - ttheta_min: min value for 2theta is simulation of peaks from cif, in degrees, guessed from GS output file if set to None)
     - ttheta_max: max value for 2theta is simulation of peaks from cif, in degrees, guessed from GS output file if set to None)
     - output: name of output file to save results in a file
     - kickoutfactor: remove grains for which the averageI is >= kickoutfactor*medianI
    
    Returns:
        a list of grain volumes, assuming that grain volume is proportionnal to diffraction intensity

    Created: 3/2021, M. Krug, Univ. Münster, S. Merkel, Univ. Lille, France
    """

    # S. Merkel: Not sure what this is for
    # g.parseInputFiles(logfile,fltfile,parfile)
    
    # Parsing indexing output files
    grains = multigrainOutputParser.parseGrains(logfile)
    ngrains = len(grains)
    print("Parsed grains from %s, found %d grains" % (logfile, ngrains))
    [peaksflt,idlist,header] = multigrainOutputParser.parseFLT(fltfile)
    npeaks = len(idlist)
    
    #with open(ciffile) as f: # Only if you saved a list of theoretical intensities first
    #    cifpeaks = f.readlines()
    #cifpeaks = cifpeaks[4:]
    
    # Preparing a list of theoretical intensities based on CIF
    if (ttheta_min is None):
        mintt = []
        for grain in grains:
            mintt.append(grain.getMinTwoTheta())
        ttheta_min = min(mintt)-1.
    if (ttheta_max is None):
        maxtt = []
        for grain in grains:
            maxtt.append(grain.getMaxTwoTheta())
        ttheta_max = max(maxtt)+1.
    cifpeaks = cifTools.peaksFromCIF(ciffile, ttheta_min,  ttheta_max, wavelength) 
    print("Calculated list of theoretical diffraction peak intensities from %s" % (ciffile))
    
    # Storing theoretical peak intensities in a dictionnary
    intensityDic = {}
    for elt in cifpeaks:
        h = int(elt[0])
        k = int(elt[1])
        l = int(elt[2])
        # trick: Pnma to Pbnm
        # h = int(elt[2])
        # k = int(elt[0])
        # l = int(elt[1])
        intensityDic[(h,k,l)] = float(elt[4])

    # Preparing a list to store normalized grain sizes    
    grainsizes = []
    grainID = 0
    
    # Loop on grains
    for grain in grains:
        grainID += 1
        peaks = grain.getPeaks()
        grainintensity = [] # A temporary list containing the normalized intensities
        for peak in peaks:
            try:
                index = idlist.index(peak.getPeakID())
                intensity = peaksflt[index]['sum_intensity']
            except IndexError:
                print ("Failed to locate peak ID %d which was found in grain %s" % (peak.getPeakID(), grain.getName()))
                return
            hkl = peak.getHKL()
            try:
                cifintensity = intensityDic[(hkl[0],hkl[1],hkl[2])] 
                relat_intensity = float(intensity) / float(cifintensity) # Normalize by theoretical intensity
                grainintensity.append(relat_intensity) # Append to list
            except KeyError:
                print ("Failed to find theoretical peak intensity for %d %d %d" % (hkl[0],hkl[1],hkl[2]))
        # Make sure that you kick out outliers with surprizingly high intensities by comparing average and median    
        average = sum(grainintensity) / len(grainintensity) # av int for each grain
        kickout = kickoutfactor*numpy.median(grainintensity) # kickoutfactor times median of each grain
        if average >= kickout:
            print ("\nIntensities of grain %s are a bit shakey\n--> Grain %s is removed from the list." % (grainID,grainID))
        else:
            grainsizes.append(average) # Collect the averages in a list

    # Print the result to the console or save into a file
    string = ""
    for item in grainsizes:
        string += str(item) + "\n"
    if output is None:
        print (string)
        # pass
    else:
        f= open(output,"w+")
        f.write(string)
        f.close()
    print ("\nDetermined relative grain sizes for %d grains.\n" % len(grainsizes))
    return grainsizes



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
    parser.add_argument('logfile', help="Path and name for GrainSpotter output file")
    parser.add_argument('fltfile', help="Path and name for fltfile")
    # parser.add_argument('parfile', help="Path and name for ImageD11 paramter file") # par file is not needed
    parser.add_argument('ciffile', help="Path and name for CIF file")
    parser.add_argument('-w', '--wavelength', required=True, help="Wavelength (angstrom, required)", type=float)
    
    # Optionnal arguments
    parser.add_argument('-m', '--ttheta_min', help="If set, minimum 2 theta for calculation of peaks from cif file (degrees). Guessed from GS output file otherwise.", type=float, default=None) # Can be guessed from the grains
    parser.add_argument('-M', '--ttheta_max', help="If set, maxium 2 theta for calculation of peaks from cif file (degrees). Guessed from GS output file otherwise.", type=float, default=None) # Can be guessed from the grains
    parser.add_argument('-o', '--output', required=False, help="If set, saves result to this file. Default is %(default)s (no filter)", default=None, type=str)
    parser.add_argument('-HV', '--histogram_vol', required=False, help="If set, plots a histogram of the grain volumes on the screen. Default is %(default)s", default=False, type=bool)
    parser.add_argument('-HR', '--histogram_rad', required=False, help="If set, plots a histogram of the grain radii on the screen. Default is %(default)s", default=False, type=bool)
    parser.add_argument('-b', '--histogram_bins', required=False, help="Sets the number of histrogram bins. Only works if histogram is set True. Default is %(default)s", default=60, type=int)
    parser.add_argument('-r', '--reject', required=False, help="Reject grains with AverageI > reject*MedianI. Default is %(default)s", default=20., type=float)
    
    # Parse arguments
    args = vars(parser.parse_args())
    ciffile = args['ciffile']
    ttheta_min = args['ttheta_min']
    ttheta_max = args['ttheta_max']
    wavelength = args['wavelength']
    logfile = args['logfile']
    fltfile = args['fltfile']
    #parfile = args['parfile']
    output = args['output']
    histogram_vol = args['histogram_vol']
    histogram_rad = args['histogram_rad']
    histogram_bins = args['histogram_bins']
    reject = args['reject']
    
    try:
        with open(ciffile) as f:
            test = f.readlines()
        print ("It´s here")
    except IOError:
        print ("File not accessible")
        
    grainsizes = grainSizeEstimate(logfile, fltfile, ciffile, wavelength, ttheta_min = ttheta_min, ttheta_max = ttheta_max, output=output, kickoutfactor=reject)
    
    # Make a histogram
    #n_bins = 60 # Determines the number of columns in the histogram
    if histogram_vol == True:
        print ("Plotting histogram ...\n")
        plt.hist(grainsizes, bins = histogram_bins)
        plt.xlabel("Grain volumes ($\mu$m^3)")
        plt.ylabel("Number of grains")
        plt.title("n = %s" % len(grainsizes), fontsize = 20)
        plt.show()
    if histogram_rad == True:
        print ("Plotting histogram ...\n")
        radii = []
        for item in grainsizes:
            radius = (3*item/4/numpy.pi())^(1/3)
            radii.append(radius)
        plt.hist(radii, bins = histogram_bins)
        plt.xlabel("Grain radii ($\mu$m)")
        plt.ylabel("Number of grains")
        plt.title("n = %s" % len(grainsizes), fontsize = 20)
        plt.show()


# Calling method 1 (used when generating a binary in setup.py)
def run():
    main(sys.argv[1:])

# Calling method 2 (if run from the command line)
if __name__ == "__main__":
    main(sys.argv[1:])
