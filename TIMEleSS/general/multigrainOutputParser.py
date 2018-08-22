#!/usr/bin/env python
# -*- coding: utf-8 -*

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

"""
Parses various files used in 3D-XRD, dealing with grains
- grainspotter output files
- gff files
"""

# Mathematical stuff (for data array)
import scipy
import scipy.linalg

# OS and file names
import os

# Specific TIMEleSS code
from TIMEleSS.general import grain3DXRD
from TIMEleSS.general import indexedPeak3DXRD

############################################################################################# 

"""
General parser for grain files, based on file extensions.

If extension is .log, will parse for GrainSpotter log files
If extension is .gff, will parse for gff file

Returns 
	A list of grains

Parameters
	filename: name and path to the gff or the GrainSpotter log file
"""

def parseGrains(filename):
	fff, file_extension = os.path.splitext(filename)
	if (file_extension == ".gff"):
		return parse_gff(filename)
	elif (file_extension == ".log"):
		return parse_GrainSpotter_log(filename)
	return []


############################################################################################# 


"""
Parser for GrainSpotter log files

Returns 
	A list of grains

Parameters
	logfile: name and path to GrainSpotter log file
"""
def parse_GrainSpotter_log(logfile):
	# Read LOG file
	f = open(logfile, 'r')
	# reads number of grains found
	text = open(logfile).read()
	ngrains = int(text.count('Grain')-1)
	# reads all the lines and saves it to the array "content"
	logcontent = [line.strip() for line in f.readlines()]
	f.close()
	# empty array containing later the line number to all grains in the file
	headgrains = []
	# empty array containing later the peak ID and hkl to all peaks found for a grain in the file
	peakInfo = []
	logpeakid = []
	# Parsing data for each grain and putting them in a f list
	grainList = []

	# looks for the word "grain" in the file by scaning each line (starting from line 20 see above) and its corresponding line number
	i=0
	lineNumb=-1
	NumbGrain = 0 
	for i, line in enumerate(logcontent):
		i+=1
		lineNumb += 1
		items = line.find("Grain")
		if (items > -1):
			# adds the line Number to the emty array where "grain" was found
			NumbGrain += 1
			headgrains.append(lineNumb)
	headgrains.remove(headgrains[0])
	for lineindex in headgrains:
		line = logcontent[lineindex]
		# Getting number of peaks
		a=line.split()
		numbpeaks = int(a[2])  
		# for the Grainnumber I have to remove the comma from the value to use it as an interger
		GrainNumW=(a[1]).strip( ',' )
		GrainNum=int(GrainNumW[0])
		grain = grain3DXRD.Grain()
		grain.setFileName(logfile)
		grain.setNPeaks(numbpeaks)
		grain.setFileIndex(int(GrainNum))
		# Extracting U matrix
		U = scipy.empty([3,3])
		line1 = logcontent[lineindex+3].split()
		line2 = logcontent[lineindex+4].split()
		line3 = logcontent[lineindex+5].split()
		for i in range (0,3):
			U[0,i] = float(line1[i])
			U[1,i] = float(line2[i])
			U[2,i] = float(line3[i])
		# Extracting UBI matrix
		UBI = scipy.empty([3,3])
		line1 = logcontent[lineindex+7].split()
		line2 = logcontent[lineindex+8].split()
		line3 = logcontent[lineindex+9].split()
		for i in range (0,3):
			UBI[0,i] = float(line1[i])
			UBI[1,i] = float(line2[i])
			UBI[2,i] = float(line3[i])
		B =scipy.linalg.inv(scipy.dot(UBI,U))
		# Setting information
		grain.setUBBi(U,B,UBI)
		# extracting the Euler angles phi1 phi phi2
		euler = logcontent[lineindex+13].split()
		grain.setEulerAngles(float(euler[0]),float(euler[1]),float(euler[2]))
		# Reading and storing some of the peak information
		peakList = []
		for i in range (0,numbpeaks):
			thispeak = indexedPeak3DXRD.indexedPeak()
			peakinfo = logcontent[lineindex+17+i].split()
			thispeak.setNum(int(peakinfo[0]))
			thispeak.setGVEID(int(peakinfo[1]))
			thispeak.setPeakID(int(peakinfo[2]))
			thispeak.setHKL(int(peakinfo[3]), int(peakinfo[4]), int(peakinfo[5]))
			thispeak.setOmegaMeasured(float(peakinfo[14]))
			thispeak.setEtaMeasured(float(peakinfo[17]))
			thispeak.setTThetaMeasured(float(peakinfo[20]))
			peakList.append(thispeak)
		grain.setPeaks(peakList)
		grainList.append(grain)
	
	return grainList


#############################################################################################     


"""
Parser for Gff files

Returns 
	A list of grains

Parameters
	logfile: name and path to the GFF file
"""
def parse_gff(gfffile):


	# Read gff file
	g = open(gfffile, 'r')
	# reads all the lines and saves it to the array "content"
	gffcontent = [line.strip() for line in g.readlines()]
	g.close()
	# Parsing data for each grain and putting them in a f list
	grainList = []
	#######	
# grain_id phase_id grainsize grainvolume x y z phi1 PHI phi2 U11 U12 U13 U21 U22 U23 U31 U32 U33 UBI11 UBI12 UBI13 UBI21 UBI22 UBI23 UBI31 UBI32 UBI33 eps11 eps12 eps13 eps22 eps23 eps33 

	gffcontent.remove(gffcontent[0])
	j= 0
	for lineindex in enumerate(gffcontent):
		j += 1
		line = gffcontent[int(j-1)].split()
		grain = grain3DXRD.Grain()
		grain.setFileName(gfffile)
		grain.setFileIndex(int(j))
		# Extracting U matrix
		U = scipy.empty([3,3])
		for i in range (0,3):
			U[0,i] = float(line[10+i])
			U[1,i] = float(line[13+i])
			U[2,i] = float(line[16+i])
		# Extracting UBI matrix
		UBI = scipy.empty([3,3])
		for i in range (0,3):
			UBI[0,i] = float(line[19+i])
			UBI[1,i] = float(line[22+i])
			UBI[2,i] = float(line[23+i])
		# Extracting B
		B =scipy.linalg.inv( scipy.dot(UBI,U))
		# Setting information
		grain.setUBBi(U,B,UBI)
		# extracting the Euler angles phi1 phi phi2
		grain.setEulerAngles(float(line[7]), float(line[8]),  float(line[9]))
		# Reading and storing peak information
		grainList.append(grain)

	return grainList

#############################################################################################

"""
Parser for FLT (peaks from diffraction data)

Returns 
	Two lists
	- peaks
	- idlist
	
	idlist is the list of "spot3d_id" for each peak
	peaks is a collection of peaks, each of them is a dictionnary will all information from the flt file

Parameters
	fname: name and path to the FLT file
"""
def parseFLT(fname):
	strings = "sc  fc  omega  Number_of_pixels  avg_intensity  s_raw  f_raw  sigs  sigf  covsf  sigo  covso  covfo  sum_intensity  sum_intensity^2  IMax_int  IMax_s  IMax_f  IMax_o  Min_s  Max_s  Min_f  Max_f  Min_o  Max_o  dety  detz  onfirst  onlast  spot3d_id  xl  yl  zl  tth  eta  gx  gy  gz"
	stringlist = strings.split()
	peaks = []
	idlist = []
	# Read file
	f = open(fname, 'r')
	for line in f.readlines():
		li=line.strip()
		peak = {}
		if not li.startswith("#"):
			litxt = li.split()
			for i in range(0,len(stringlist)):
				peak[stringlist[i]] = litxt[i]
			thisid = int(peak["spot3d_id"])
			idlist.append(thisid)
			peaks.append(peak)
	f.close()
	print ("Parsed list of peaks from flt file %s, found %i peaks" % ( fname, len(peaks)))
	return [peaks,idlist]
