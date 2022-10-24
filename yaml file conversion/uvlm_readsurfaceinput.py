import numpy as np
import re

def uvlm_readsurfaceinput():
	fid1 = open("surfaceinput.txt","r")
    fid = -1	
    if (fid1 != -1):
        fid = fid1
    if (fid >= 0): 
        for i in range (1,4):
            tline = fid.readline().strip()
        numsurf = int(fid.readline().strip())
        for k in range (1,numsurf + 1):
            for i in range (1,4):
                tline = fid.readline.strip()
            tline = fid.readline.strip()
            # Implementation of sscanf funtion in python !!!
            a = list(map(int,tline.split("")))
            nn = a[1], ne = a[2], nnx = a[3], nny = a[4]
            for i in range (1,4):
                tline = fid.readline.strip()				
            for i in range (1,nn + 1):
                tline = fid.readline.strip()					
                tline = re.sub('d','e',tline)
                #Should we consider value of k or index of k !!!
                surface[k][(coord[i,1:4])] = list(map(float,tline.split("")))
            for i in range (1,4):
                tline = fid.readline.strip()				
            for i in range (1,ne + 1):
                tline = fid.readline.strip()
                surface[k][(inz[i,1:4])] = list(map(int,tline.split("")))
	return
"""   
 def file_reading(file):
    newline_break = ""
    for readline in reading:
        line_strip = line.rstrip('\n')
        newline_break += line_strip

def function(file):
    lines = []
    for line in f:
        lines.append(line)
    return lines 
"""