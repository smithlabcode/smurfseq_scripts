#!/usr/bin/env python
'''     
*   Copyright (C) 2018
*
*   Authors: Timour Baslan. et. al.
*
*   This program is free software: you can redistribute it and/or modify
*   it under the terms of the GNU General Public License as published by
*   the Free Software Foundation, either version 3 of the License, or
*   (at your option) any later version.
*
*   This program is distributed in the hope that it will be useful,
*   but WITHOUT ANY WARRANTY; without even the implied warranty of
*   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
*   GNU General Public License for more details.
*
*   You should have received a copy of the GNU General Public License
*   along with this program.  If not, see <http://www.gnu.org/licenses/>.
'''
'''
* GPL added by rishvanth
'''

import sys


def main():

	if (len(sys.argv) < 6):
		print "./get_bin_counts.py [sam_file] [chrom_info] [bin_boundaries] [out] [out_stats]" 
		sys.exit()

	infilename = sys.argv[1]
	chrominfo_file = sys.argv[2]
	bins_file = sys.argv[3]
	outfilename = sys.argv[4]
	statfilename = sys.argv[5]
  
	chrominfo = fileToDictionary(chrominfo_file, 0)
	bins = fileToArray(bins_file, 0)
	INFILE = open(infilename, "r")
	OUTFILE = open(outfilename, "w")
	STATFILE = open(statfilename, "w")

	binCounts = []
	for i in range(len(bins)):
		binCounts.append(0)

	# print len(binCounts)
	# print len(bins)

	counter = 0
	totalReads = 0
	prevChrompos = ""
	for x in INFILE:
		if (x[0] == '@'): continue
		arow = x.rstrip().split("\t")
		thisChrom = arow[2]
		thisChrompos = arow[3]
		if thisChrom.find("_") > -1:
			#print thisChrom
			continue
		if thisChrom == "chrM":
			#print thisChrom
			continue
		if thisChrom == "":
			continue
		if chrominfo.has_key(thisChrom):
			pass
		else:
			continue

		totalReads += 1
			
		thisChrominfo = chrominfo[thisChrom]
		thisAbspos = long(thisChrompos) + long(thisChrominfo[2])
		
		counter += 1
		
		indexUp = len(bins) - 1
		indexDown = 0
		indexMid = int((indexUp - indexDown) / 2.0)

		while True:
			if thisAbspos >= long(bins[indexMid][2]):
				indexDown = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexMid
			else:
				indexUp = indexMid + 0
				indexMid = int((indexUp - indexDown) / 2.0) + indexDown

			if indexUp - indexDown < 2:
				break

		binCounts[indexDown] += 1
		prevChrompos = thisChrompos
		
	for i in range(len(binCounts)):
		thisRatio = float(binCounts[i]) / (float(counter) / float(len(bins)))
		OUTFILE.write("\t".join(bins[i][0:3]))
		OUTFILE.write("\t")
		OUTFILE.write(str(binCounts[i]))
		OUTFILE.write("\t")
		OUTFILE.write(str(thisRatio))
		OUTFILE.write("\n")

	binCounts.sort()
	
	STATFILE.write("TotalReads\tMedianBinCount\n")
	STATFILE.write(str(totalReads))
	STATFILE.write("\t")
	STATFILE.write(str(binCounts[len(bins)/2]))
	STATFILE.write("\n")

	INFILE.close()
	OUTFILE.close()
	STATFILE.close()


def fileToDictionary(inputFile, indexColumn):
	input = open(inputFile, "r")

	rd = dict()
#	input.readline()
	for x in input:
		arow = x.rstrip().split("\t")
		id = arow[indexColumn]
		if rd.has_key(id):
			#rd[id].append(arow)
			print "duplicate knowngene id = " + id
			print "arow =   " + str(arow)
			print "rd[id] = " + str(rd[id])
		else:
			rd[id] = arow
		
	input.close()
	return(rd)


def fileToArray(inputFile, skipFirst):
	input = open(inputFile, "r")

	ra = []

	for i in range(skipFirst):
		input.readline()

	for x in input:
		arow = x.rstrip().split("\t")
		ra.append(arow)
		
	input.close()
	return(ra)


if __name__ == "__main__":
	main()
