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

if(len(sys.argv) != 2):
  print "./get_unique_maps.py [sam file]"
  sys.exit()

mapfile = open(sys.argv[1], 'r')

# process the sam file
while True:
  mapline = mapfile.readline()
  if not mapline: break
  if (mapline[0] == '@'): 
    print mapline.rstrip()
    continue 
  line = mapline.rstrip().split("\t")
  # print mapline

  if (int(line[3]) == 0):
    continue

  found = 0
  if (int(line[4]) < 2):
    for i in range(0,len(line)):
      if (line[i][0:2] == 'XA'):
        found = 1

  if (found == 0):
    print mapline.rstrip()
 

mapfile.close()
