import re
import numpy as np
from matplotlib import pyplot as plt
def uvlm_readmodel():
	model = []
	fid1 = open("simulationinput_aero.txt","r")
	fid = -1
    range_4to6 = np.array([4,5,6])
		if (fid1 != -1):
			fid = fid1
		if (fid >= 0): 
			for i in range (1,4):
				tline = fid.readline()		
			tline = fid.readline().strip()
			str_n = list(map(str,tline.split("")))
			string = str_n
            count = len(str_n)
			if (count>1):
				str_n = list(str(str_n))
				string = dict(str_n[1][1])
				string = string[1]			
			model["strSimName"] = string
			for i in range (1,4):
				tline = fid.readline()			
			tline = fid.readline().strip()
            tline = re.sub('d','e',tline)
			model["simulationsettings"] = list(map(float,tline.split("")))
			for i in range (1,4):
				tline = fid.readline()                
			tline = fid.readline().strip()
            tline = re.sub('d','e',tline[9:-1])
			model["windData"] = list(map(float,tline.split(" ")))
			model["density"]  = model["windData"][1]
			model["winddir"]  = model["windData"][range_4to6]
			model["vinf"]     = model["windData"][2]*model["winddir"]
			plt.close(fid)
		
	fid3 = open("simulationinput_fsi.txt","r")
	if (fid3 >= 0):
		fid = fid3
		for i in range (1,4):
			tline = fid.readline()		
		tline = fid.readline().strip()
		str_n = list(map(str,tline.split("")))
        string = str_n
        count = len(str_n)
		if (count > 1):
			str_n = list(str(str_n))
            string = dict(str_n[1][1])
            string = string[1]		
		model["strSimName"] = string
		for i in range (1,4):
			tline = fid.readline()
		tline = fid.readline().strip()
        tline = re.sub('d','e',tline)
		model["simulationsettings"] = list(map(float,tline.split("")))
		for i in range (1,4):
			tline = fid.readline()
		tline = fid.readline().strip()
        tline = re.sub('d','e',tline[9:-1])
		model["windData"] = list(map(float,tline.split("")))
		model["density"]  = model["windData"][1]
		model["winddir"]  = model["windData"][range_4to6]
		model["vinf"]     = model["windData"][2]*model["winddir"]       
		plt.close(fid3)

	fid1 = open(model["strSimName"] + '_uvlm_models.dres',"r")
	fid2 = open(model["strSimName"] + '_uvlm_models.txt',"r")
	fid = -1
    range_1to3 = np.array([1,2,3])
		if (fid1 != -1):
			fid = fid1
		elseif (fid2 != -1):
			fid = fid2
	if (fid >= 0): 
		tline = readline(fid)
		model["nsurfaces"] = sscanf(tline,'%i')
		i3 = 0
		for i in range (1,model["nsurfaces"]+1):
			tline = fid.readline().strip()
			arrline = list(map(int,tline.split("")))
			model["surfaces"][i]["nnodes"] = arrline[1]
			model["surfaces"][i]["nnodes"] = arrline[4]
			for j in range (1,model["surfaces"][i]["nelements"] + 1):
				tline = fid.readline().strip()
				model["surfaces"][i]["connectivity"][j,:] = list(map(int,tline.split("")))
			for j in range (1,model["surfaces"][i]["nnodes"]):				
                model["surfaces"][i]["nodes"][j]["indices_q"] = [3*(j-1)+(range_1to3)*(j-1)+3] + i3           
			i3 = i3 + 3*model["surfaces"][i]["nnodes"]		
	plt.close(fid)
	return
    
    
