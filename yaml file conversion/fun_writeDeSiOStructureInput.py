import os
import shutil
import numpy as np
 
def fun_writeDeSiOStructureInput(simu,mesh):

    currDir = os.getcwd()

    # creating structure input files and directories
    caseDir = simu["currDir"] + "\\"  +  simu["strfilename"] + '\\DeSiO-Aero\\'
    os.mkdi r(caseDir)

    # checking, if DeSiO.bat exists and 
	flag = "currDir" + '\\' + "DeSiO.bat"
	if exists(flag):
	  shutil.copy(flag,caseDir)
	else:
	  warningMessage = print(f"Warning: file does not exist:\n {flag}")

    # writing input files
    os.chdir(caseDir)
    with open('beaminput.txt','w') as fid:
        print("!! \n", file = fid)
        print("!! \n", file = fid)
        print("!!number of beams (1), number of cross-section properties (2), flag for property input type (3):\n", file = fid) 
        print("%i\t%i\t%i\n"%(mesh[2],simu["matbeam"][2], simu["flag_matbeam"]), file = fid)
        arr_node = []
        nnodes = 0
        for i_s in range (1,np.shape(mesh[2], ndmin = 2)+1):
            strbeamname = ''            
            if "strname" in mesh["i_s"]:
                strbeamname = mesh["i_s"]["strname"]       
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! beam %s: number of nodes (1), number of elements (2):\n"%('strbeamname'), file = fid)
            print("%i\t%i\n"%(mesh(i_s)["nn mesh(i_s)"]["mx"]), file = fid)
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! beam %s: nodal phi (1, 2, 3), nodal d1 (4, 5, 6), nodal d2 (7, 8, 9), nodal d3 (10, 11, 12):\n"%('strbeamname'), file - fid)
            for i in range (1,np.shape(mesh[i_s]["nodes"], ndmin = 1)+1):
                print("%20.15fd0\t"%([mesh[i_s]["nodes"][i,1:12)]), file = fid)
                print("!!\n", file = fid)
            
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! beam %s: connectivities (1, 2), cross-section property (3):\n"%('strbeamname'), file = fid)
            for i in range (1,np.shape(mesh[i_s]["connectivity"], ndmin = 1)):
                print("%i\t%i\t%i\n"%(mesh["i_s"]["connectivity"][i,:]), file = fid)        
            
            # global nodes
            nnodes = nnodes + mesh["i_s"][nn]
            arr_node[None+1:nnodes] = [len[arr_node]+1:nnodes]
        
        if (np.shape(simu["matbeam"]) !=0): 
            for i_mat in range (1, np.shape(simu["matbeam"],ndmin=2)):
                strbeamname = ''
                if "strname" in simu["matbeam"][il_mat]:
                    strbeamname = simu["matbeam"][i_mat]["strname"]                        
                if (simu["flag_matbeam"] == 1): # general DeSiO input-format (Voigt notation)
                    print("!! cross-section property %i (%s) (Voigt notation):\n"%(i_mat,strbeamname), file = fid) 
                    print("!! row1: cbeam\n", file = fid)
                    print("!! row2: cmass\n", file = fid)
                else: # isotropic material DeSiO input-format 
                    print("!! cross-section property %i (%s):\n"%(i_mat,strbeamname), file = fid)
                    print("!! row1: EA(1), GA1(2), GA2(3), EI1(4), EI2(5), GI3(6), ES1(7), ES2(8), GS1(9), GS2(10), EI12(11)\n", file = fid)
                    print("!! row2: rhoA3(1), rhoI1(2), rhoI2(3), rhoS1(4), rhoS2(5), rhoI12(6)\n", file = fid)
                
                
                string = re.sub('e','d',print(f"{simu['matbeam']['i_mat']['cmat']: .5f}"))
                print('%s'%('string'), file = fid) 
                print('\n', file = fid)
                
                
                string = re.sub('e','d',print(f"{simu['matbeam']['i_mat']['mmat']: .5f}"))
                print('%s'%('string'), file = fid) 
                print('\n', file = fid)
                
                string = re.sub('e','d',print(f"{simu['matbeam']['i_mat']['diss']: .5f}"))
                print('%s'%('string'), file = fid) 
                print('\n', file = fid)
            
        
        fid.close()

    if 'pointmass12' in simu:
        with open('pointmass12input.txt','w') as fid:
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print('!! number of point masses\n', file = fid)
            print('%i\n'%(simu["pointmass12"][2]), file = fid)
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print('!! (1) node, (2) point mass:\n', file = fid)
            for ipm12 in range (1,np.shape(simu.pointmass12,2)
                print("%i\t%10.5fd0\n"%(simu["pointmass12"][ipm12]["node"],simu["pointmass12"][ipm12]["mass"]), file = fid)
            
        fid.close()
    

    if "constraints" in simu:
        with open('constraint12input.txt','w') as fid:
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! number of constraints for nodes with 12 coordinates:\n", file = fid)
            print("%i\n"%(np.shape(simu["constraints"],ndmin2)), file = fid)
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! constraints for nodes with 12 coordinates: sort (1), nodes (2, 3), phi1 (4, 5, 6), phi2 (7, 8, 9), dir (10, 11, 12):\n", file = fid)
            for ic in range (1,np.shape(simu["constraints"],ndmin=2)+1):
                phi1 = [0,0,0]
                phi2 = [0,0,0]
                direc  = [0,0,0]
                print("%s\t"%(simu["constraints"][ic]["types"]), file = fid)
                print("%s\t"%(simu["constraints"][ic]["nodes"]), file = fid)
                if "phi1" in simu["constraints"][ic]:
                    if not (simu["constraints"][ic]["phi1"]):
                        phi1 = simu["constraints"][ic]["phi1"]                
                if "phi2" in simu["constraints"][ic]:
                    if not (simu["constraints"][ic]["phi2"]):
                        phi1 = simu["constraints"][ic]["phi2"] 
                if "dir" in simu["constraints"][ic]:
                    if not (simu["constraints"][ic]["direc"]):
                        phi1 = simu["constraints"][ic]["direc"] 
                
                
    #             print(fid,'%10.8fd0\t',phi1)
    #             print(fid,'%10.8fd0\t',phi2)
    #             print(fid,'%10.8fd0\t',direc)
                string = re.sub('e','d',print(f"{phi1: .8f}"))
                print("%s"%('string'), file = fid)
                
                string = re.sub('e','d',print(f"{phi2: .8f}"))
                print("%s"%('string'), file = fid)
                
                string = re.sub('e','d',print(f"{direc: .8f}"))
                print("%s"%('string'), file = fid)
                
                print("\n", file = fid)
            
        fid.close()
    

    if "rb" in simu:
        with open('rigidbodyinput.txt','w') as fid:
            print("!! rigid body input\n", file = fid)
            print("!!\n", file = fid)
            print("!! number of rigid bodies (1), number of body properties (2):\n", file = fid)
            print("%i\t%i\n"%(np.shape(simu["rb"],ndmin=2), np.shape(simu["rb"],ndmin=2)), file = fid)
            for irb in range (1,np.shape(simu["rb"],ndmin=2)+1):
                types = ''
                print("!!\n", file = fid)
                print("!!\n", file = fid)
                if "types" in simu["rb"][irb]:
                    types = simu["rb"][irb]["types"]
                print("!! rigid body %i, %s phi(1:3), d1(4:6) d2(7:9), d3(10:12)\n"%(irb,'types'), file = fid)
                print("%10.8fd0\t"%(simu["rb"][irb]["centre_mass"]), file = fid)
                print("%10.8fd0\t"%(simu["rb"][irb]["D1"]), file = fid)
                print("%10.8fd0\t"%(simu["rb"][irb]["D2"]), file = fid)
                print("%10.8fd0\t"%(simu["rb"][irb]["D3"]), file = fid)
                print("\n", file = fid)
                print("!!\n", file = fid)
                print("!!\n", file = fid)
                print("!! rigid body %i, %s property\n"%(irb,'types'),file = fid)
                print("%i\n"%(irb), file = fid)
            
            for irb in range (1,np.shape(simu["rb"],ndmin=2)+1):
                types = ''
                print("!!\n", file = fid)
                print("!!\n", file = fid)            
                if "types" in simu["rb"][irb]:
                    types = simu["rb"][irb]["types"]
                print("!! rigid body %i property, %s hub mass(1), J11(2), J22(3), J33(4), J23(5), J13(6), rhoxhi3(7), rhoxhi2(8), rhoxhi1(9), J12(10)\n"%(irb,'types'), file = fid)
                print("%10.8fd0\t"%(simu["rb"][irb]["mass_matrix"]), file = fid)
                print("\n", file = fid)
        
        fid.close()    
    

    with open('simulationinput_structure.txt','w') as fid:
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! filename\n", file = fid)
        print("solution\n", file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! number of simulations:\n", file = fid)
        print("%i\n"%(1), file = fid)
        print("!! Settings for nonlinear dynamic solver with constant time stepping\n", file = fid)
        print("!! row1: simulation settings: totalt (1), deltat (2), tolerance (3), iteration limit (4), gravity flag (5)\n", file = fid)
        print("!! row2: flag for writing matrices - on/off = 1/0\n", file = fid)
        print("modal\n", file = fid)
        print("10\t1.0d-8\t0.0d0\t1.0d8\t0\n", file = fid)
        print("0\n", file = fid)
        print("!!\n", file = fid)
        print("!!\n", file = fid)
        print("!! gravity vector (1, 2, 3):\n", file = fid)
        print("%10.5fd0"%(simu["grav_vec"]), file = fid)
        print("\n", file = fid)
    fid.close()

    if (np.shape(simu["loads"],ndmin=1)) != 0:
        with open('load12input.txt','w') as fid:
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! number of loads for nodes with 12 coordinates.\n", file=fid)
            print("%i\n"%(np.shape(simu["loads"], ndmin =1)), file = fid)
            print("!!\n", file = fid)
            print("!!\n", file = fid)
            print("!! loads for nodes with 12 coordinates: sort (1), intensity (2), duration (3), node (4), spatial (5, 6, 7, 8, 9, 10), material (11, 12, 13, 14, 15, 16).\n", file = fid)
            for i in range(1,np.shape(simu["loads"],ndmin=1)+1):
                print("%s"%(simu["loads"][i]["sort"]), file = fid)
                print("%20.15fd0\t%20.15fd0\t%i\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\n"%([simu["loads"][i]["intensity"],simu["loads"][1]["duration"],simu["loads"][1]["node"],simu["loads"][1]["spatial"],simu["loads"][1]["material"]), file = fid)
        
        fid.close()
    

    
    os.chdir(simu["currDir"])
    print('creating DeSiO-Structure input files')
    
    return