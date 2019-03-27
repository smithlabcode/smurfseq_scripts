#!/usr/bin/env python

## Takes a sam file of reads generated from gen_semisim_frags.py
## This program generated semi-simulated reads from the mapped
## fragments. 

from __future__ import print_function
from random import randint
import sys
import pysam

if (len(sys.argv) != 6):
  print("./foo [sam file] [frag count] [frag len] [read count] [out file]")
  sys.exit()

frag_count = int(sys.argv[2])
frag_len = int(sys.argv[3])
read_len = frag_count * frag_len
read_count = int(sys.argv[4])

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



seq_list = list()
pos_list = list()

map_file = pysam.AlignmentFile(sys.argv[1], 'r', check_sq=False)
for read in map_file.fetch():
  # Make sure the alignmnet is unique 
  if (read.has_tag('XA') == False and read.has_tag('SA') == False and 
      read.is_unmapped == False and read.mapping_quality > 1):
    # print(read.query_name)

    read_name = read.query_name
    print(read.query_name)
    read_map_ch = read_name.split("-")[1].split(":")[0]
    read_map_pos = int(read_name.split("-")[1].split(":")[1])
    if (read_map_ch == read.reference_name):
      if (abs(read_map_pos - read.reference_start) <= 10):
  
        seq = read.query_alignment_sequence
        if len(seq) == frag_len:
          frag_ch = read.reference_name 
          frag_pos = read.reference_start
          frag_seq = seq
          if (read.is_reverse):
            frag_seq = rev_comp(frag_seq)

          # print(frag_ch, frag_pos, read_name, frag_seq)
          seq_list.append(frag_seq)
          pos_list.append((frag_ch, str(frag_pos)))


print(len(seq_list))
print(len(pos_list))


out_file = open(sys.argv[5], 'w')
count = 0
while (count < read_count) and (len(seq_list) != 0):
  length = 0
  out_seq = ''
  out_pos = ''
  while length < read_len and (len(seq_list) != 0):
    j=randint(0,len(seq_list)-1)
    out_seq += seq_list[j]
    out_pos += (pos_list[j][0]+":"+pos_list[j][1]+",")
    length = length + len(seq_list[j])
    del seq_list[j]
    del pos_list[j]
  
  out_file.write(">Read%d" %(count) + " " + out_pos[0:len(out_pos)-1] + "\n")
  out_file.write(out_seq + "\n")

  count += 1
  

out_file.close()
