import shutil
import os
from os.path import exists
import numpy as np


def fun_writeDeSiOAeroInput(simu,mesh,wake):

	currDir = os.getcwd()
	# creating uvlm input files and directories
	caseDir = simu["currDir"] + "\\"  +  simu["strfilename"] + '\\DeSiO-Aero\\'
	os.mkdir(caseDir)

	# checking, if DeSiO.bat exists
	flag = "currDir" + '\\' + "DeSiO.bat"
	if exists(flag):
	  shutil.copy(flag,caseDir)
	else:
	  warningMessage = print(f"Warning: file does not exist:\n {flag}")
	

	# writing input files
	os.chdir(caseDir)
	with open("surfaceinput.txt","w") as fid:
        print("!! \n", file = fid)
        print("!! \n", file = fid)
        print("!! number of surfaces (1):\n", file = fid)
        print("%i \n"%(np.shape(mesh,ndmin=2)), file = fid)
        for i_s in range (1,np.shape(mesh,ndmin=2)+1):
            print("!! \n", file = fid)
            print("!! \n", file = fid)
            print("!! !! surface: number of nodes (1), number of rings (2), number of nodes along first dimension (3), and number of nodes along second direction (4). \n", file = fid)
            print(f"{mesh[i_s][nn] {mesh[i_s][mx]*mesh[i_s][my]} {mesh[i_s][mx]+1} {mesh[i_s][my]+1}\n"], file = fid)
            print("!! \n", file = fid)
            print("!! \n", file = fid)
            print("!! surface node coordinates: nodal phi (1, 2, 3). \n", file = fid)
            for i in range (1,np.shape(mesh[i_s]["nodes"],ndmin=1)):
                print("%20.10fd0 %20.10fd0 %20.10fd0 \n"%(mesh[i_s]["nodes"][i:1], mesh[i_s]["nodes"][i:2], mesh[i_s]["nodes"][i:3]), file = fid)
            
            print("!! \n", file = fid)
            print("!! \n", file = fid)
            print("!! surface: connectivities (1, 2, 3, 4). \n", file = fid)
            for i in range (1,np.shape(mesh[i_s]["connectivity"],ndmin=1)):
                print("%i %i %i %i \n"%(mesh[i_s]["connectivity"][i,:], file = fid)        
		
	
	fid.close()

	if not wake:
		n_wakes = 0
		nprop   = 0
	else:
		n_wakes = np.shape(wake["wakes"],ndmin=2)
		nprop   = np.shape(wake["property"],ndmin=1)

	  
	with open("wakeinput.txt","w") as fid:
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! number of wakes (1), number of wake properties (2)\n", file = fid)
        print("%i\t%i\n"%(n_wakes, nprop), file = fid)
        for i in range (1,n_wakes+1):
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! wake %i: number of surface (1), number of segments (2), property number (3)\n"%(i), file = fid)
            print("%i\t%i\t%i\n"%(wake["wakes"][i]["nsurf"], wake["wakes"][i]["nsegments"], wake["wakes"][i]["nproperty"]), file = fid)
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! wake %i: nodes (1, 2), ring (3)\n"%(i), file = fid)
            for j in range (1,np.shape(wake["wakes"][i]["inf"],ndmin=1)):
                print("%i\t%i\t%i\n"%(wake["wakes"][i]["inf"][j,:]), file = fid)
            
        
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! wake property %i: time cut (1), max. rows (2)\n"%(i), file = fid)
        for i in range (1,nprop +1):
            print("%i\t%i\t0.0d0\t\n"%(wake["property"][i,:]), file = fid)
        
	fid.close()

	with open("simulationinput_aero.txt","w") as fid:
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! filename\n", file = fid)
        print("solution\n", file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! simulation settings: totalt (1), deltat (2), cutoff (3)\n", file = fid)
        print("%20.15fd0\t%20.15fd0\t%20.15fd0\n"%(simu["time"],simu["deltat"],simu["cutoff"]), file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! wind data: sort(1), fluid density(2), intensity(3), duration(4), direction vector (5-7)\n", file = fid)
        print("constant\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\n"%(simu["density"], simu["i_vinf"], simu["time"], simu["d_vinf"]), file = fid)
	fid.close()
	
    os.chdir(simu["currDir"])
	print('creating DeSiO-Aero input files')
	
	return