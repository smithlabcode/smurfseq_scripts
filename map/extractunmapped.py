#!/usr/bin/env python
'''     
*   Copyright (C) 2018 University of Southern California
*
*   Authors: Rishvanth Prabakar
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

import sys
import pysam
from operator import itemgetter

if(len(sys.argv) != 4):
  print "./extractunmapped.py [mapped sam file] [read file.fastq] [min read len]"
  sys.exit()

mapfile = pysam.AlignmentFile(sys.argv[1], "r", check_sq=False)
minlen = int(sys.argv[3])

# process the map file
samcigar = {}
for read in mapfile.fetch():
  if(read.cigartuples == None): continue
  # print read.query_name
  # print read.query_alignment_length
  start = 0
  # print read.flag
  # print read.cigartuples
  if ((read.flag >> 4) & 1):
    if(read.cigartuples[len(read.cigartuples)-1][0]==4 or read.cigartuples[len(read.cigartuples)-1][0]==5):
      start = read.cigartuples[len(read.cigartuples)-1][1]
    else:
      start = 0
  else: 
    if(read.cigartuples[0][0]==4 or read.cigartuples[0][0]==5):
      start = read.cigartuples[0][1]
    else:
      start = 0
 
  # print start, read.query_alignment_length 
 
  if(samcigar.has_key(read.query_name)):
    samcigar[read.query_name].append([start, read.query_alignment_length])
  else:
    samcigar[read.query_name] = [[start, read.query_alignment_length]]

mapfile.close()

# process the read file and print the unmapped regions
readfile = open(sys.argv[2], 'r')
while True:
  nameline = readfile.readline()
  seqline = readfile.readline()
  dummyline = readfile.readline()
  qualline = readfile.readline()
  if not nameline: break
  name = nameline.rstrip().split(" ")[0]
  mapped = samcigar.get(name[1:len(name)]) 
  if (mapped == None):
    print nameline[0:len(nameline)-1]
    print seqline[0:len(seqline)-1]
    print dummyline[0:len(dummyline)-1]
    print qualline[0:len(qualline)-1]
  else:
    mapped.sort(key=itemgetter(0))
    # print mapped
    readlen = len(seqline)-1
    # print readlen
    gapstart = 0
    gapend = mapped[0][0]
    j=0
    k=0
    for i in mapped:
      # print i
      if((gapend-gapstart) > minlen):
        # print "st:", gapstart, "end:", gapend
        # print name+"-"+str(k)
        print name
        print seqline[gapstart:gapend]
        print dummyline[0:len(dummyline)-1]
        print qualline[gapstart:gapend]   
        k += 1

      gapstart = gapend + mapped[j][1]
      if((j+1)==len(mapped)):
        gapend = readlen
      else:
        gapend = mapped[j+1][0]
      j += 1 
      # print gapend
    if((gapend-gapstart) > minlen):
      # print "st:", gapstart, "end:", gapend
      # print name+"-"+str(k)
      print name
      print seqline[gapstart:gapend]
      print dummyline[0:len(dummyline)-1]
      print qualline[gapstart:gapend]

readfile.close()
