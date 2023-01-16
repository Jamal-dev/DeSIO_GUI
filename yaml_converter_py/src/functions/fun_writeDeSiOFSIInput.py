import os
import shutil
from pathlib import Path
from os.path import exists
import logging
import numpy as np
import re
def fun_writeDeSiOFSIInput(simu,simu_struc,simu_aero):
    currDir = os.getcwd()

    # creating uvlm input files and directories
    caseDir = Path(simu.currDir) /  Path(simu.strfilename) / Path('DeSiO')
    if not os.path.exists(caseDir):
        os.makedirs(caseDir)

    source_file = currDir / Path("DeSiO-FSI.bat")
    if exists(source_file):
        shutil.copy(source_file,caseDir)
    else:
        logging.warning(f"DeSiO-FSI.bat not found:\n {source_file}")
    

    os.chdir(caseDir)
    with open("simulationinput_fsi.txt","w") as fid:
        print("!!", file = fid)
        print("!!", file = fid)

        structprevjobname = 'none'
        aeroprevjobname   = 'none'
        if 'structprevjobname' in dir(simu):
            structprevjobname = simu.structprevjobname
        if 'aeroprevjobname' in dir(simu):
            aeroprevjobname = simu.aeroprevjobname


        print("!! filename (1), prev filename aero (2), prev filename structure (3):", file = fid)
        print("%s\t%s\t%s"%(simu.strfilename,aeroprevjobname,structprevjobname), file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! total time (1) deltat aero (2) cutoff (3) time for load factor (4) type of load factor function (5)", file = fid)
        
        lf_duration = 0.0
        if "lf_duration" in dir(simu):
            lf_duration = simu.lf_duration
        
        lf_type = 'constant'
        if "lf_type" in dir(simu):
            lf_type = simu.lf_type
        
        #for box in simu_aero.time:
        simu_aero_t = re.sub('e','d',"{:.5e}".format(simu_aero.time))
            #print('%s'%(string), end =" ", file = fid) 

        #for box in simu_aero.deltat:
        simu_aero_d = re.sub('e','d',"{:.5e}".format(simu_aero.deltat))
            #print('%s'%(string), end =" ", file = fid) 

        #for box in simu_aero.cutoff:
        simu_aero_cut = re.sub('e','d',"{:.5e}".format(simu_aero.cutoff))
            #print('%s'%(string), end =" ", file = fid) 

        #for box in lf_duration:
        lf_d = re.sub('e','d',"{:.5e}".format(lf_duration))
            #print('%s'%(string), end =" ", file = fid) 


        print("%s\t%s\t%s\t%s\t%s"%(simu_aero_t,simu_aero_d,simu_aero_cut,lf_d,lf_type), file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! wind data: sort(1), fluid density(2), intensity(3), duration(4), direction vector (5-7)", file = fid)

        density = 1.0; i_vinf  = 1.0; time = 1.0; d_vinf = [1,0,0]; str_sort = 'constant';

        if 'density' in dir(simu_aero):
            density = simu_aero.density
        if 'i_vinf' in dir(simu_aero):
            i_vinf = simu_aero.i_vinf
        if 'time' in dir(simu_aero):
            time = simu_aero.time
        if 'd_vinf' in dir(simu_aero):
            d_vinf = simu_aero.d_vinf
        if 'sort' in dir(simu_aero):
            str_sort = simu_aero.sort

        strrep_density = re.sub('e','d',"{:.8e}".format(density))
        strrep_i_vinf = re.sub('e','d',"{:.8e}".format(i_vinf))
        strrep_time = re.sub('e','d',"{:.8e}".format(time))
        
        print("%s\t%s\t%s\t%s"%(str_sort,strrep_density,strrep_i_vinf,strrep_time), end = " ", file = fid)
        for box in d_vinf:
            strrep_d_vinf = re.sub('e','d',"{:.8e}".format(box))
            print("\t%s"%(strrep_d_vinf), end = " ", file = fid)
        print("", file = fid)

        print("!!", file = fid)
        print("!!", file = fid)
        print("!! deltat (1), tolerance (2), iteration limit (3), gravity flag (4), row2: flag for writing matrices - on/off = 1/0", file = fid)
        
        struc_simType = 'dynamic'
        if "simType" in dir(simu_struc):
            struc_simType = simu_struc.simType
        
        print("%s"%(struc_simType), file = fid)

        simu_struc_d = re.sub('e','d',"{:.5e}".format(simu_struc.deltat))
        simu_struc_tol = re.sub('e','d',"{:.5e}".format(simu_struc.tol))


        print("%s\t%s\t%i\t%i"%(simu_struc_d,simu_struc_tol,simu_struc.niter,simu_struc.grav), file = fid)
        
        flag_linearization = 0
        if "flag_linearization" in dir(simu):
            flag_linearization = simu.flag_linearization
        
        print("%i\t%i"%(0,flag_linearization), file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! gravity vector (1, 2, 3)", file = fid)
        #print("%20.10fd0\t%20.10fd0\t%20.10fd0"%(simu_struc.grav_vec[0],simu_struc.grav_vec[1],simu_struc.grav_vec[2]), file = fid)
        for box in simu_struc.grav_vec:
            string = re.sub('e','d',"{:.5e}".format(box))
            print('%s'%(string), end =" ", file = fid) 
        print("",file = fid)
        #inz_Ext = []
        if 'file' in simu_aero.sort:
            if 'wind_field_file' in dir(simu_aero):
                source_field_file = Path(currDir) / Path(simu_aero.wind_field_file)
                if exists(source_field_file):
                    shutil.copy(source_field_file,caseDir)
                else:
                    logging.warning(f"Warning: Wind data file does not exist:\n{source_field_file}")
                print(f'!!', file = fid)
                print(f'!! if sort = file then consider these addtional lines for inflow settings', file = fid)
                print(f'!! filename for wind field: filename with file extension (1)', file = fid)
                
                wind_field_file_ext = simu_aero.wind_field_file.split('.')[-1]
                if 'wnd' in wind_field_file_ext:
                    wind_field_type = 1
                elif 'bts' in wind_field_file_ext:
                    wind_field_type = 2
                else:
                    wind_field_type = 1
                    logging.warning(f"Warning: wrong file type for wind_field_file! Needed (.wnd or .bts)")
                
                print(f'{simu_aero.wind_field_file}', file = fid)
                print(f'!!', file = fid)
                print(f'!! inflow settings', file = fid)
                print(f'!! grids center (1-3)', file = fid)
                grid_center = np.array([0,0,0])
                if hasattr(simu_aero,'grid_center') and not type(simu_aero.grid_center) is list:
                    logging.debug(f"Reading grid center from simu_aero\ngrid_center: {simu_aero.grid_center}")
                    grid_center = np.asarray(simu_aero.grid_center)
                logging.debug(f'grid_center: {grid_center}')
                if grid_center.ndim>1:
                    grid_center = grid_center.ravel()

                for box in grid_center:
                    strrep_grid_center = re.sub('e','d',"{:.8e}".format(box))
                    print(f"\t{strrep_grid_center}", end = " ", file = fid)
            else:
                logging.warning(f"Warning: no field: wind_field_file found!")


    with open("fsi_input.txt","w") as fid:
        """
        Here the value of simu.data[1][2] and simu.data[2][2] are the same.
        For rectification, I think we should change the root values of simu.data[i][2]
        """
        nfsi = 0
        if "data" in dir(simu):
            nfsi = np.shape(simu.data)[0]
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! number of fsi (1)", file = fid)
        print("%i"%(nfsi), file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! input for fluid-structure interaction: from shell/beaminput (1), surface from surfaceinput (2) local search radius (3)", file = fid)
        for i in range (0,nfsi):
            if isinstance(simu.data[0][3], str):
                print("%s\t%i\t%i\t%s"%(simu.data[i][0],simu.data[i][1],simu.data[i][2],simu.data[i][3]), file = fid)
            else:
                print("%s\t%i\t%i\t%20.10fd0"%(simu.data[i][0],simu.data[i][1],simu.data[i][2],simu.data[i][3]), file = fid)



    if 'data_input' in dir(simu):
        if simu.data_input != 0:
            for i in range (0,nfsi):
                with open("simu.data[i][3]") as fid:
                    print(f'!! Input search radius for weights using radial-based function.\n', file = fid)
                    print(f'!! fsi {i}: {simu.data[i][0]} {simu.data[i][1]} to surface {simu.data[i][2]}\n', file = fid)
                    print(f'!! number of nodes for fsi search radius\n', file = fid)
                    print(f'{np.shape(simu.data_input(i).node_fsi_radius)[0]}', file = fid)
                    print(f'!!\n', file = fid)
                    print(f'!!\n', file = fid)
                    print(f'!! fsi search radius (1)\n', file = fid)
                    strrep_fsi_radius = re.sub('e','d',"{:.8e}".format(simu.data_input(i).node_fsi_radius))
                    print("\t%s"%(strrep_fsi_radius), end = " ", file = fid)
                    print("", file = fid)

    os.chdir(simu.currDir)
    print('created DeSiO-FSI input files')