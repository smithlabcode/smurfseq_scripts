#!/usr/bin/env python

## Generates semi-sim frags from mapped WGS reads.
## These reads need to be mapped before generating 
## semi-sim reads

from __future__ import print_function
from random import randint
import sys
import pysam

if (len(sys.argv) != 3):
  print("./gen_semisim_frags.py [sam file] [frag len]")
  sys.exit()

frag_len = int(sys.argv[2])

def rev_comp(in_seq):
  out_seq = ''
  for i in in_seq:
    if i == 'A':
      out_seq += 'T'
    elif i == 'T':
      out_seq += 'A'
    elif i == 'G':
      out_seq += 'C'
    elif i == 'C':
      out_seq += 'G'
    else:
      out_seq += 'N'
  return out_seq[::-1] 


read_count = 0

map_file = pysam.AlignmentFile(sys.argv[1], 'r', check_sq=False)
for read in map_file.fetch():
  # Make sure the alignmnet is unique 
  if (read.has_tag('XA') == False and read.has_tag('SA') == False and 
      read.is_unmapped == False and read.mapping_quality >= 1):  
    # print(read.query_name)
    seq = read.query_alignment_sequence
    if (len(seq) > frag_len):
  
      frag_st = randint(0, frag_len)
      # frag_st = 0

      frag_ch = read.reference_name 

      cig_readlen = 0
      num_ins = 0
      num_del = 0
      tup = read.cigartuples
      # print(tup)
      i = 0
      while (frag_st + frag_len < len(seq)):
        while (cig_readlen <= frag_st):    
          if (tup[i][0] == 0 or tup[i][0] == 1):
            cig_readlen += tup[i][1]
        
          if (tup[i][0] == 1):
            num_ins += tup[i][1]
          elif (tup[i][0] == 2):
            num_del += tup[i][1]     
          i += 1

        frag_pos = read.reference_start + frag_st - num_ins + num_del    
        frag_seq = seq[frag_st:frag_st+frag_len]
        if (read.is_reverse):
          frag_seq = rev_comp(frag_seq)

        # print(frag_ch, frag_pos, frag_seq)
        # print(">Read" + str(read_count) + "-" + frag_ch + ":" + str(frag_pos))
        # print(frag_seq)
        print(frag_ch, str(frag_pos), str(frag_pos+frag_len), frag_seq, sep="\t")
        frag_st += frag_len
        read_count += 1
