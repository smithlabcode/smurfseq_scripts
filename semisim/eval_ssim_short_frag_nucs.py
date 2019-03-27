#!/usr/bin/env python

## Used to evaluate mapped semi-simulated results.
## Generates (1) precision and recall of mapped fragment
## location, (2) fraction of unmapped nucleoties after mapping,
## and (3) mean, median and SD of mapped fragment lengths
 

import sys
import pysam
import numpy
from operator import itemgetter

if(len(sys.argv) != 9):
  print "./evalmapper [read file.fa] [mapped sam file] [frag len]"
  sys.exit()

readfile = open(sys.argv[1], 'r')
mapfile = pysam.AlignmentFile(sys.argv[2], 'r', check_sq=False)
frag_len = int(sys.argv[3])

## These passed in just for printing. The program does not make
## use of these
mat = int(sys.argv[4])
mmat = int(sys.argv[5])
opn = int(sys.argv[6])
indel = int(sys.argv[7])
min_as = int(sys.argv[8])

# process the sam file
readmap = {}
fraglen = []
mapped_frags = {}


for read in mapfile.fetch():
  if (read.is_unmapped): continue

  seq = read.query_alignment_sequence
  if not seq: continue 

  start = 0
  if (read.is_reverse):
    if(read.cigartuples[len(read.cigartuples)-1][0]==4 or read.cigartuples[len(read.cigartuples)-1][0]==5):
      start = read.cigartuples[len(read.cigartuples)-1][1] - 1 # sam file uses 1-offset
    else:
      start = 0
  else: 
    if(read.cigartuples[0][0]==4 or read.cigartuples[0][0]==5):
      start = read.cigartuples[0][1] - 1
    else:
      start = 0
 
  if(mapped_frags.has_key(read.query_name)):
    mapped_frags[read.query_name].append([start, read.query_alignment_length])
  else:
    mapped_frags[read.query_name] = [[start, read.query_alignment_length]]


  
  if (read.mapping_quality == 0): continue
  qname = read.query_name
  qchr = read.reference_name
  qpos = read.reference_start
  # print qname, qchr, qpos, seq
  fraglen.append(read.query_alignment_length)
  if (readmap.has_key(qname)):
    readmap[qname].append([qchr, qpos, len(seq)])
  else: 
    readmap[qname] = [[qchr, qpos, len(seq)]]
  

mapfile.close()


# process the read file
numfrag = 0
numfp = 0
numfn = 0
numtp = 0

while True:
  nameline = readfile.readline()
  dummyline = readfile.readline()
  if not nameline: break
  name = nameline.rstrip().split(" ")[0]
  pos = nameline.rstrip().split(" ")[1].split(",")
  mappos = readmap.get(name[1:len(name)], [])
  numfrag += len(pos)
  # print name
  # print len(pos)
  # print len(mappos)
  nummatch = 0
  for i, rval in enumerate(pos):
    rch = rval.split(":")[0]
    rpos = int(rval.split(":")[1])
    err = frag_len / 2
    # print rch, rpos
    matchfound = False
    for j, mval in enumerate(mappos):
      mch = mval[0]
      mpos = int(mval[1])
      mlen = int(mval[2])
      # print "m: ", mch, mpos
      if(rch == mch):
        overlap = 0
        if (mpos < rpos):
          overlap = mpos + mlen - rpos
        else:
          if (mpos+mlen < rpos+frag_len):
            overlap = mlen
          else:
            overlap = rpos + frag_len - mpos

        # if((mpos <= rpos+err) and (mpos+mlen >= rpos+err)):
        # if (abs(mpos - rpos) <= err):
        if (overlap >= frag_len/2):
          mval[1] = '0'
          numtp += 1
          nummatch += 1
          matchfound = True
          # print "match found"
          break
    if(matchfound == False):
      # print "match not found"
      numfn += 1
  numfp += (len(mappos) - nummatch)

# print "Total fragments:", numfrag
# print "True positives:", numtp
# print "False positives:", numfp
# print "False negatives:", numfn
# print "Precision:", float(numtp)/float(numtp+numfp)
# print "Recall:", float(numtp)/float(numtp+numfn)
# print "Total maps:", numtp+numfp 


readfile.seek(0)

unmap_bases = 0
total_bases = 0

while True:
  nameline = readfile.readline().rstrip()
  seqline = readfile.readline().rstrip()
  if not nameline: break
  total_bases += len(seqline)
  name = nameline.split(" ")[0]
  mapped = mapped_frags.get(name[1:len(name)]) 
  if (mapped == None):
    # print (name, len(seqline), 0, len(seqline))   
    unmap_bases += len(seqline)
  else:
    mapped.sort(key=itemgetter(0))
    readlen = len(seqline)
    gapstart = 0
    gapend = mapped[0][0]
    j=0
    k=0
    for i in mapped:
      if((gapend-gapstart) > 0):
        # print (name, len(seqline), gapstart, (gapend-gapstart))   
        unmap_bases += (gapend-gapstart)
        k += 1

      gapstart = gapend + mapped[j][1]
      if((j+1)==len(mapped)):
        gapend = readlen
      else:
        gapend = mapped[j+1][0]
      j += 1 
      # print gapend
    if((gapend-gapstart) > 0):
      # print (name, len(seqline), gapstart, (gapend-gapstart))   
      unmap_bases += (gapend-gapstart)



readfile.close()


pres = round(float(numtp)/float(numtp+numfp), 4)
rec = round(float(numtp)/float(numtp+numfn), 4)
mean_len = round(numpy.mean(fraglen), 4)

std_len = round(numpy.std(fraglen), 4)
median_len = round(numpy.median(fraglen), 4)

num_maps = numtp+numfp
frac_unmap = round((float(unmap_bases)/float(total_bases)) * 100, 4)

print min_as, mat, mmat, opn, indel, pres, rec, num_maps, frac_unmap, mean_len, std_len, median_len 
