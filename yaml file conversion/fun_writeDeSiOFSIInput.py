import os
import shutil


def fun_writeDeSiOFSIInput(simu,simu_struc,simu_aero):
    
    currDir = os.getcwd()

    # creating uvlm input files and directories
    caseDir = simu["currDir"] + "\\"  +  simu["strfilename"] + '\\DeSiO\\'
	os.mkdir(caseDir)

    flag = "currDir" + '\\' + "DeSiO.bat"
	if exists(flag):
	  shutil.copy(flag,caseDir)
	else:
	  warningMessage = print(f"Warning: file does not exist:\n {flag}")
.
    os.chdir(caseDir)
	with open("simulationinput_fsi.txt","w") as fid:
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! filename (1), prev filename aero (2), prev filename structure (3)\n:", file = fid)
        print("%s\n"%(simu["strfilename"]), file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! total time (1) deltat aero (2) cutoff (3) time for load factor (4) type of load factor function (5)\n", file = fid)
        
        lf_duration = simu_aero["time"]
        if "lf_duration" in simu:
            lf_duration = simu["lf_duration"]
        
        lf_type = 'constant'
        if "lf_type" in simu:
            lf_type = simu["lf_type"]
        
        print("%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%s\n"%(simu_aero["time"],simu_aero["deltat"],simu_aero["cutoff"],lf_duration,lf_type), file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! wind data: sort(1), fluid density(2), intensity(3), duration(4), direction vector (5-7)\n", file = fid)
        print("%s\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\n"%('constant',simu_aero["density"],simu_aero["i_vinf"],simu_aero["time"],simu_aero["d_vinf"]), file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! deltat (1), tolerance (2), iteration limit (3), gravity flag (4), row2: flag for writing matrices - on/off = 1/0\n", file = fid)
        
        struc_simType = 'dynamic'
        if "simType" in simu_struc:
            struc_simType = simu_struc["simType"]
        
        print("%s\n"%('struc_simType'), file = fid)
        print("%20.10fd0\t%20.10fd0\t%i\t%i\n"%(simu_struc["deltat"],simu_struc["tol"],simu_struc["niter"],simu_struc["grav"]), file = fid)
        
        flag_linearization = 0
        if "flag_linearization" in simu:
            flag_linearization = simu["flag_linearization"]
        
        print("%i\t%i\n"%(0,flag_linearization), file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! gravity vector (1, 2, 3)\n", file = fid)
        print("%20.10fd0\t%20.10fd0\t%20.10fd0\n"%(simu_struc["grav_vec"]), file = fid)
    fid.close()

    with open("fsi_input.txt","w") as fid:
        nfsi = 0
        if "data" in simu:
            nfsi = np.shape(simu["data"],ndmin=1)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! number of fsi (1)\n", file = fid)
        print("%i\n"%(nfsi), file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! input for fluid-structure interaction: from shell/beaminput (1), surface from surfaceinput (2) local search radius (3)\n", file = fid)
        for i in range (1:nfsi + 1):
            print("%s\t%i\t%i\t%20.10fd0\n"%(simu["data"][i,1],simu["data"][i,2],simu["data"][i,3],simu["data"][i,4]), file = fid)
    fid.close()

    os.chdir(simu["currDir"])
	print('creating DeSiO-FSI input files')
    
    return