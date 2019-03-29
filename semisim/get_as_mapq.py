#!/usr/bin/env python

from __future__ import print_function
import sys
import pysam

if (len(sys.argv) != 4):
  print("./mappedfraglen [mapped sam file] [min as] [min mapq]")
  sys.exit()



map_file = pysam.AlignmentFile(sys.argv[1], 'r', check_sq=False)
min_as = int(sys.argv[2])
min_mapq = int(sys.argv[3])

of_name = sys.argv[1][0:len(sys.argv[1])-4] + ".as" + str(min_as) + \
  ".hq" + str(min_mapq) + ".sam" 
as_maps = pysam.AlignmentFile(of_name, "w", template=map_file)

for read in map_file.fetch():
  if (read.mapping_quality >= min_mapq):
    for i in read.get_tags():
      if (i[0] == "AS"):
        if (i[1] >= min_as):
          as_maps.write(read)    
          break


map_file.close()
as_maps.close()
