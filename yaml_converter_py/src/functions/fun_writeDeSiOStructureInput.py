import os
import shutil
from pathlib import Path
from os.path import exists
import logging
import numpy as np
import re
def fun_writeDeSiOStructureInput(simu,mesh):
    currDir = os.getcwd()

    # creating structure input files and directories
    if 'jobname' in dir(simu):
        caseDir = Path(simu.currDir) /  Path(simu.strfilename) / Path('DeSiO-Structure') / Path(simu.jobname)
    else:
        caseDir = Path(simu.currDir) / Path(simu.strfilename) / Path('DeSiO-Structure')
    if not os.path.exists(caseDir):
        os.makedirs(caseDir)

    # checking, if DeSiO.bat exists and 
    source_file = Path(currDir) /   Path("DeSiO-FSI.bat")
    if exists(source_file):
        shutil.copy(source_file,caseDir)
    else:
        logging.warning(f"DeSiO-FSI.bat not found:\n {source_file}")

    # writing input files
    os.chdir(caseDir)
    with open('beaminput.txt','w') as fid:
        print("!!", file = fid)
        print("!!", file = fid)
        print("!!number of beams (1), number of cross-section properties (2), flag for property input type (3):", file = fid) 
        print("%i\t%i\t%i"%(np.shape(mesh)[0],np.shape(simu.matbeam)[0], int(simu.flag_matbeam)), file = fid)
        arr_node = []
        nnodes = 0
        adding = 1
        incre = 1
        for i_s in range (0,np.shape(mesh)[0]):
            strbeamname = ''
            if "strname" in dir(mesh[i_s]):     #index value is not correct
                strbeamname = mesh[i_s].strname
            if "blade" in strbeamname:
                strbeamname = 'blade' + ' ' + str(incre) #But changed the value using increment operator
                incre = incre + 1 
            print("!!", file = fid)
            print("!!", file = fid)
            print("!! beam %s: number of nodes (1), number of elements (2):"%(strbeamname), file = fid)
            print("%i\t%i"%(mesh[i_s].nn, mesh[i_s].mx), file = fid)
            print("!!", file = fid)
            print("!!", file = fid)
            print("!! beam %s: nodal phi (1, 2, 3), nodal d1 (4, 5, 6), nodal d2 (7, 8, 9), nodal d3 (10, 11, 12):"%(strbeamname), file = fid)
            for i in range (0,np.shape(mesh[i_s].nodes)[0]):
                for ip in range (0, np.shape(mesh[0].nodes)[1]):
                    print("%20.15fd0\t"%(mesh[i_s].nodes[i][ip]),end = " ",file = fid) 
                print("!!", file = fid) 
            print("!!", file = fid)
            print("!!", file = fid)
            print("!! beam %s: connectivities (1, 2), cross-section property (3):"%(strbeamname), file = fid)
            for i in range (0,np.shape(mesh[i_s].connectivity)[0]):
                if "blade" in strbeamname:
                    #For the blade vlaues to be same, I had to hardcode the below specific line
                    print("%i\t%i\t%i"%(mesh[i_s].connectivity[i][0]+1, mesh[i_s].connectivity[i][1]+1, mesh[i_s].connectivity[i][2]+3), file = fid)  
                else:
                    print("%i\t%i\t%i"%(mesh[i_s].connectivity[i][0]+1, mesh[i_s].connectivity[i][1]+1, mesh[i_s].connectivity[i][2]+adding), file = fid)       
            # global nodes
            nnodes = nnodes + mesh[i_s].nn
            arr_node.append([*range(len(arr_node)+1,nnodes+1)])
            adding = adding + 1
        if np.shape(simu.matbeam)!=0:
            for i_mat in range (0, np.shape(simu.matbeam)[0]):
                strbeamname = ''
                if "strname" in dir(simu.matbeam[i_mat]):
                    strbeamname = simu.matbeam[i_mat].strname
                if (int(simu.flag_matbeam) == 1): # general DeSiO input-format (Voigt notation)
                    print("!! cross-section property %i (%s) (Voigt notation):"%(i_mat,strbeamname), file = fid) 
                    print("!! row1: cbeam", file = fid)
                    print("!! row2: cmass", file = fid)
                else: # isotropic material DeSiO input-format 
                    print("!! cross-section property %i (%s):"%(i_mat,strbeamname), file = fid)
                    print("!! row1: EA(1), GA1(2), GA2(3), EI1(4), EI2(5), GI3(6), ES1(7), ES2(8), GS1(9), GS2(10), EI12(11)", file = fid)
                    print("!! row2: rhoA3(1), rhoI1(2), rhoI2(3), rhoS1(4), rhoS2(5), rhoI12(6)", file = fid)
                for box in simu.matbeam[i_mat].cmat:
                    string = re.sub('e','d',"{:.5e}".format(box))
                    print('%s'%(string), end =" ", file = fid)
                print("", file = fid)        

                for box in simu.matbeam[i_mat].mmat:
                    string = re.sub('e','d',"{:.5e}".format(box))
                    print('%s'%(string), end =" ", file = fid) 
                print("", file = fid)

                for box in simu.matbeam[i_mat].diss:
                    string = re.sub('e','d',"{:.5e}".format(box))
                    print('%s'%(string), end =" ", file = fid)
                print("", file = fid) 
    if "pointmass12" in dir(simu):
        with open('pointmass12input.txt','w') as fid:
            print("!!", file = fid)
            print("!!", file = fid)
            print('!! number of point masses', file = fid)
            print('%i'%(np.shape(simu.pointmass12)[0]), file = fid)
            print("!!", file = fid)
            print("!!", file = fid)
            print('!! (1) node, (2) point mass:', file = fid)
            for ipm12 in range (0,np.shape(simu.pointmass12)[0]):
                print("%i\t%10.5fd0"%(simu.pointmass12[ipm12].node+1,simu.pointmass12[ipm12].mass), file = fid)

    if "constraints" in vars(simu).keys():
        if simu.constraints:
            with open('constraint12input.txt','w') as fid:
                print("!!", file = fid)
                print("!!", file = fid)
                print("!! number of constraints for nodes with 12 coordinates:", file = fid)
                print("%i"%(np.shape(simu.constraints)[0]), file = fid)
                print("!!", file = fid)
                print("!!", file = fid)
                print("!! constraints for nodes with 12 coordinates: sort (1), nodes (2, 3), phi1 (4, 5, 6), phi2 (7, 8, 9), dir (10, 11, 12):", file = fid)
                logging.debug(f'np.shape(simu.constraints) = {np.shape(simu.constraints)}')
                # TODO: Check Surya
                for ic in range (0,len(simu.constraints)): #np.shape(simu.constraints)[0] = 190 but it isn't the same in latest matlab code
                    phi1 = [0,0,0]
                    phi2 =  [0,0,0]
                    dic  =  [0,0,0]
                    print("%s\t"%(simu.constraints[ic].type), end =" ", file = fid)
                    for ik in range(0,np.shape(simu.constraints[ic].nodes)[0]):
                        print("%s\t"%(simu.constraints[ic].nodes[ik]), end =" ", file = fid) #Values are incorrect: precedes value by 1 and different for each case

                    #if np.array([]) is not None: print(1)

                    if "phi1" in dir(simu.constraints[ic]):
                        if len(simu.constraints[ic].phi1) != 0:
                            phi1 = simu.constraints[ic].phi1.ravel() 
                    if "phi2" in dir(simu.constraints[ic]):
                        if len(simu.constraints[ic].phi2) != 0:
                            phi2 = simu.constraints[ic].phi2.ravel()
                    if "dir" in dir(simu.constraints[ic]):
                        if len(simu.constraints[ic].dir) != 0:
                            dic = simu.constraints[ic].dir.ravel()

                    for box in phi1:
                        string = re.sub('e','d',"{:.8e}".format(box))
                        print('%s'%(string), end =" ", file = fid) 
                    
                    for box in phi2:
                        string = re.sub('e','d',"{:.8e}".format(box))
                        print('%s'%(string), end =" ", file = fid)
                    
                    for box in dic:
                        string = re.sub('e','d',"{:.8e}".format(box))
                        print('%s'%(string), end =" ", file = fid)
                    print("", file = fid)

        if "rb" in dir(simu):
            with open('rigidbodyinput.txt','w') as fid:
                print("!! rigid body input", file = fid)
                print("!!", file = fid)
                print("!! number of rigid bodies (1), number of body properties (2):", file = fid)
                print("%i\t%i"%(np.shape(simu.rb)[0], np.shape(simu.rb)[0]), file = fid)
                for irb in range (0,np.shape(simu.rb)[0]):
                    types = ''
                    print("!!", file = fid)
                    print("!!", file = fid)
                    if "type" in dir(simu.rb[irb]):
                        types = simu.rb[irb].type
                    print("!! rigid body %i, %s phi(1:3), d1(4:6) d2(7:9), d3(10:12)"%(irb+1,types), file = fid)

                    #Append all arrays for typing in a single line
                   #Append all arrays for typing in a single line
                    for nk in range (0,np.shape(simu.rb[irb].center_mass)[0]):
                        if np.shape(simu.rb[irb].center_mass)[0] != 1: 
                            print("%10.8fd0\t"%(simu.rb[irb].center_mass[nk]), end = " ", file = fid)    
                        else:
                            print()
                            for op in range (0,np.shape(simu.rb[irb].center_mass)[1]):
                                print("%10.8fd0\t"%(simu.rb[irb].center_mass[nk][op]), end = " ", file = fid)

                    for nk in range (0,np.shape(simu.rb[irb].D1)[0]):
                        if np.shape(simu.rb[irb].D1)[0] != 1: 
                            print("%10.8fd0\t"%(simu.rb[irb].D1[nk]), end = " ", file = fid)    
                        else:
                            #print()
                            for op in range (0,np.shape(simu.rb[irb].D1)[1]):
                                print("%10.8fd0\t"%(simu.rb[irb].D1[nk][op]), end = " ", file = fid)

                    for nk in range (0,np.shape(simu.rb[irb].D2)[0]):
                        if np.shape(simu.rb[irb].D2)[0] != 1: 
                            print("%10.8fd0\t"%(simu.rb[irb].D2[nk]), end = " ", file = fid)    
                        else:
                            #print()
                            for op in range (0,np.shape(simu.rb[irb].D2)[1]):
                                print("%10.8fd0\t"%(simu.rb[irb].D2[nk][op]), end = " ", file = fid)

                    for nk in range (0,np.shape(simu.rb[irb].D3)[0]):
                        if np.shape(simu.rb[irb].D3)[0] != 1: 
                            print("%10.8fd0\t"%(simu.rb[irb].D3[nk]), end = " ", file = fid)    
                        else:
                            #print()
                            for op in range (0,np.shape(simu.rb[irb].D3)[1]):
                                print("%10.8fd0\t"%(simu.rb[irb].D3[nk][op]), end = " ", file = fid)

                    print("\n!!", file = fid)
                    print("!!", file = fid)
                    print("!! rigid body %i, %s property"%(irb+1,types),file = fid)
                    print("%i"%(irb+1), file = fid)
                for irb in range (0,np.shape(simu.rb)[0]):
                    types = ''
                    print("!!", file = fid)
                    print("!!", file = fid)            
                    if "type" in dir(simu.rb[irb]):
                        types = simu.rb[irb].type
                    print("!! rigid body %i property, %s hub mass(1), J11(2), J22(3), J33(4), J23(5), J13(6), rhoxhi3(7), rhoxhi2(8), rhoxhi1(9), J12(10)"%(irb+1,types), file = fid)
                    for nk in range (0,np.shape(simu.rb[irb].mass_matrix)[0]):
                        print("%10.8fd0\t"%(simu.rb[irb].mass_matrix[nk]), end = " ", file = fid)
                    print("", file = fid)
    
    strtypecomments = ['row 1: simulation settings: totalt (1), deltat (2), tolerance (3), iteration limit (4), gravity flag (5)',   # dynamic and static solver
                   'row 1: simulation settings: number of EV (1), tolerance (2), emin, emax, flag (0 - dense; 1 - sparse)',  # modal and linear buckling solver
                   'row 1: simulation settings: totalt (1), deltat (2), tolerance (3), iteration limit (4), arc length method (5), desired iterations (6)'] # static post-buckling solver

    if 'type' in dir(simu):
        strjobname     = 'solution'
        strprevjobname = 'none'
        if 'dynamic' in simu.type:
            if 'jobname' in dir(simu):
                strjobname = simu.jobname
            if 'prevjobname' in dir(simu):
                strprevjobname = simu.prevjobname
            strtype                  = 'dynamic'
            strtypecomment           = strtypecomments[0]
            str_simu_settings_format = ('%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\t%i\n')
            simu_settings            = [1.0e0, 1.0e-1, 1.0e-6, 50, 0]
            if 'settings' in dir(simu):
                simu_settings = simu.settings
            simu_grav_vec = [0,0,0]
            if 'grav_vec' in dir(simu):
                simu_grav_vec = simu.grav_vec
        elif 'static' in simu.type:
            if 'jobname' in dir(simu):
                strjobname = simu.jobname
            if 'prevjobname' in dir(simu):
                strprevjobname = simu.prevjobname
            strtype                  = 'static'
            strtypecomment           = strtypecomments[0]
            str_simu_settings_format = ('%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\t%i\n')
            simu_settings            = [1.0e0, 1.0e-1, 1.0e-6, 50, 0]
            if 'settings' in dir(simu):
                simu_settings = simu.settings
            simu_grav_vec = [0,0,0]
            if 'grav_vec' in dir(simu):
                simu_grav_vec = simu.grav_vec
        elif 'modal' in simu.type:
            if 'jobname' in dir(simu):
                strjobname = simu.jobname
            if 'prevjobname' in dir(simu):
                strprevjobname = simu.prevjobname
            strtype                  = 'modal'
            strtypecomment           = strtypecomments[1]
            str_simu_settings_format = ('%i\t%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\n')
            simu_settings            = [10, 1.0e-6, 0.0, 1.0e0, 0]
            if 'settings' in dir(simu):
                simu_settings = simu.settings
            simu_grav_vec = [0,0,0]
        elif 'buckling' in simu.type:
            if 'jobname' in dir(simu):
                strjobname = simu.jobname
            if 'prevjobname' in dir(simu):
                strprevjobname = simu.prevjobname
            strtype                  = 'buckling'
            strtypecomment           = strtypecomments[1]
            str_simu_settings_format = ('%i\t%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\n')
            simu_settings            = [10, 1.0e-6, 0.0, 1.0e0, 0]
            if 'settings' in dir(simu):
                simu_settings = simu.settings
            simu_grav_vec = [0,0,0]
            if 'grav_vec' in dir(simu):
                simu_grav_vec = simu.grav_vec
        elif 'static_arc' in simu.type:
            if 'jobname' in dir(simu):
                strjobname = simu.jobname
            if 'prevjobname' in dir(simu):
                strprevjobname = simu.prevjobname
            strtype                  = 'static_arc'
            strtypecomment           = strtypecomments[2]
            str_simu_settings_format = ('%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\t%i\n')
            simu_settings            = [1.0e0, 1.0e0, 1.0e-6, 50, 2, 4]
            if 'settings' in dir(simu):
                simu_settings = simu.settings
            simu_grav_vec = [0,0,0]
            if 'grav_vec' in dir(simu):
                simu_grav_vec = simu.grav_vec
    else:
        strjobname               = 'solution'
        strprevjobname          = 'none'
        strtype                  = 'modal'
        strtypecomment           = strtypecomments[1]
        str_simu_settings_format = ('%i\t%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\n')
        simu_settings            = [10, 1.0e-6, 0.0, 1.0e0, 0];
        #simu_settings_format     = print("%i	%10.8fd0	%10.8fd0	%10.8fd0	%i"%(simu_settings[0], simu_settings[1], simu_settings[2], simu_settings[3], simu_settings[4])) 
        simu_grav_vec            = [0,0,0]

    with open('simulationinput_structure.txt','w') as fid:
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! filename", file = fid)
        print(f'{strjobname}\t{strprevjobname}', file = fid) 
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! number of simulations:", file = fid)
        print("%i"%(1), file = fid)
        print(f'!! ... place for comments ...', file = fid)
        print(f'!! {strtypecomment}', file = fid)
        print(f'!! row2: flag for writing matrices - on/off = 1/0', file = fid)
        print(f'{strtype}', file = fid)
        print("%i	%10.8fd0	%10.8fd0	%10.8fd0	%i"%(simu_settings[0], simu_settings[1], simu_settings[2], simu_settings[3], simu_settings[4]), file = fid)
        print(f'0', file = fid)
        print(f'!!', file = fid)
        print(f'!!', file = fid)
        print(f'!! gravity vector (1, 2, 3):', file = fid)
        for box in simu_grav_vec:
            string = re.sub('e','d',"{:.5e}".format(box))
            print('%s'%(string), end =" ", file = fid)
        #print(f'%10.5fd0'%(simu_grav_vec), file = fid)
        print(f'', file = fid)
        """
        print("!! Settings for nonlinear dynamic solver with constant time stepping", file = fid)
        print("!! row1: simulation settings: totalt (1), deltat (2), tolerance (3), iteration limit (4), gravity flag (5)", file = fid)
        print("!! row2: flag for writing matrices - on/off = 1/0", file = fid)
        print("modal", file = fid)
        print("10\t1.0d-8\t0.0d0\t1.0d8\t0", file = fid)
        print("0", file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! gravity vector (1, 2, 3):", file = fid)
        for ib in range (0, len(simu.grav_vec)):
            print("%10.5fd0"%(simu.grav_vec[ib]), sep = "\t",end = "\t", file = fid)
        """
    if 'loads' in dir(simu):
        with open('load12input.txt','w') as fid:
            print("!!", file = fid)
            print("!!", file = fid)
            print("!! number of loads for nodes with 12 coordinates.", file=fid)
            print("%i"%(np.shape(simu.loads)[0]), file = fid)
            print("!!", file = fid)
            print("!!", file = fid)
            print("!! loads for nodes with 12 coordinates: sort (1), intensity (2), duration (3), node (4), spatial (5, 6, 7, 8, 9, 10), material (11, 12, 13, 14, 15, 16).", file = fid)
            for i in range(0,np.shape(simu.loads)[0]):
                print("%s"%(simu.loads[i].sort), file = fid)
                print("%20.15fd0\t%20.15fd0\t%i\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0"%(simu.loads[i].intensity,simu.loads[1].duration,simu.loads[1].node,simu.loads[1].spatial,simu.loads[1].material), file = fid)
    
    if 'boundary12' in dir(simu):
        with open('boundary12input.txt','w') as fid:
            print(f'!!\n', file = fid)
            print(f'!!\n', file = fid)
            print(f'!! number of inhomogeneous constraints for nodes with 12 coordinates.\n', file = fid)
            print(f'{np.shape(simu.boundary)[0]}', file = fid)
            print(f'!!\n', file = fid)
            print(f'!! inhomogeneous constraints for nodes with 12 coordinates: \n', file = fid)
            print(f'!! sort (1), intensity (2), duration (3), constraint ID (4), file (5) (if no, then none).\n', file = fid)
            for i in range (0,np.shape(simu.boundary)[0]):
                print(f'{simu.boundary12(i).sort}\t{simu.boundary12(i).intensity: 20.15f}d0\t{simu.boundary12(i).time: 20.15f}d0\t{simu.boundary12(i).constraintID}\t{simu.boundary12(i).file}')
    os.chdir(simu.currDir)
    print('creating DeSiO-Structure input files')