import os
from pathlib import Path
from os.path import exists
import shutil
import logging
import numpy as np
import re
def fun_writeDeSiOAeroInput(simu,mesh,wake):
    """ Write the input file for DeSiOAero
    """
    currDir = os.getcwd()
    # creating uvlm input files and directories
    caseDir = Path(simu.currDir) /  Path(simu.strfilename) / Path('DeSiO-Aero')
    if not exists(caseDir):
        os.makedirs(caseDir)

    # checking, if DeSiO.bat exists
    source_file = currDir / Path("DeSiO.bat")
    if exists(source_file):
        shutil.copy(source_file,caseDir)
    else:
        logging.warning(f"DeSiO.bat not found in current directory:\n {source_file}")
    os.chdir(caseDir)
    with open("surfaceinput.txt","w") as fid:
        print("!! ", file = fid)
        print("!! ", file = fid)
        print("!! number of surfaces (1):", file = fid)
        #To match the number of surfaces as 6, I subtracted the length by 2. 
        #Actual output: 8
        print("%i"%(np.shape(mesh)[0]), file = fid)
        for i_s in range (0, len(mesh)): #2, len(mesh) is correct
            print("!! ", file = fid)
            print("!! ", file = fid)
            print("!! !! surface: number of nodes (1), number of rings (2), number of nodes along first dimension (3), and number of nodes along second direction (4). ", file = fid)
            print(f"{mesh[i_s].nn} {mesh[i_s].mx*mesh[i_s].my} {mesh[i_s].mx+1} {mesh[i_s].my+1}", file = fid)
            print("!! ", file = fid)
            print("!! ", file = fid)
            print("!! surface node coordinates: nodal phi (1, 2, 3). ", file = fid)
            for i in range (0,np.shape(mesh[i_s].nodes)[0]):
                print("%20.10fd0 %20.10fd0 %20.10fd0 "%(mesh[i_s].nodes[i][0], mesh[i_s].nodes[i][1], mesh[i_s].nodes[i][2]), file = fid)
            
            print("!! ", file = fid)
            print("!! ", file = fid)
            print("!! surface: connectivities (1, 2, 3, 4). ", file = fid)
            for i in range (1,np.shape(mesh[i_s].connectivity)[0]):
                print("%i %i %i %i "%(mesh[i_s].connectivity[i][0],mesh[i_s].connectivity[i][1],mesh[i_s].connectivity[i][2], mesh[i_s].connectivity[i][3]), file = fid)
    if not wake:
        n_wakes = 0
        nprop   = 0
    else:
        n_wakes = np.shape(wake.wakes)[0]
        nprop   = np.shape(wake.property)[0]    
  
    with open("wakeinput.txt","w") as fid:
        """
        The values within 'wake' are changing alterantively, and it seems that the origin value (value of attributes within wake) is something to be changed
        rather than tampering within function of writing the input file.   
        """
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! number of wakes (1), number of wake properties (2)", file = fid)
        print("%i\t%i"%(n_wakes, nprop), file = fid)
        for i in range (0,n_wakes):
            print("!!", file = fid)
            print("!!", file = fid)
            print("!! wake %i: number of surface (1), number of segments (2), property number (3)"%(i), file = fid)
            print("%i\t%i\t%i"%(wake.wakes[i].nsurf, wake.wakes[i].nsegments, wake.wakes[i].nproperty), file = fid)
            print("!!", file = fid)
            print("!!", file = fid)
            print("!! wake %i: nodes (1, 2), ring (3)"%(i), file = fid)
            for j in range (0,np.shape(wake.wakes[i].inf)[0]):
                print("%i\t%i\t%i"%(wake.wakes[i].inf[j][0]+1,wake.wakes[i].inf[j][1]+1,wake.wakes[i].inf[j][2]+1), file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! wake property %i: time cut (1), max. rows (2)"%(i), file = fid)
        for i in range (0,nprop):
            print("%i\t%i\t0.0d0\t"%(wake.property[i][0],wake.property[i][1]), file = fid)
    # writing simulationinput_aero.txt
    with open("simulationinput_aero.txt","w") as fid:
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! filename", file = fid)
        print("solution", file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! simulation settings: totalt (1), deltat (2), cutoff (3)", file = fid)
        aero_setting = [simu.time,simu.deltat,simu.cutoff]

        for box in aero_setting:
            strrep_aero_setting = re.sub('e','d',"{:.8e}".format(box))
            print("%s"%(strrep_aero_setting), end = " ", file = fid)
        print("", file = fid)

        #print("%20.15fd0\t%20.15fd0\t%20.15fd0"%(simu.time,simu.deltat,simu.cutoff), file = fid)
        print("!!", file = fid)
        print("!!", file = fid)
        print("!! wind data: sort(1), fluid density(2), intensity(3), duration(4), direction vector (5-7)", file = fid)

        density = 1.0; i_vinf  = 1.0; time = 1.0; d_vinf = [1,0,0]; str_sort = 'constant';
        if 'density' in dir(simu):
            density = simu.density
        if 'i_vinf' in dir(simu):
            i_vinf = simu.i_vinf
        if 'time' in dir(simu):
            time = simu.time
        if 'd_vinf' in dir(simu):
            d_vinf = simu.d_vinf
        if 'sort' in dir(simu):
            str_sort = simu.sort

        strrep_density = re.sub('e','d',"{:.8e}".format(density))
        strrep_i_vinf = re.sub('e','d',"{:.8e}".format(i_vinf))
        strrep_time = re.sub('e','d',"{:.8e}".format(time))
        
        print("%s\t%s\t%s\t%s"%(str_sort,strrep_density,strrep_i_vinf,strrep_time), end = " ", file = fid)
        for box in d_vinf:
            strrep_d_vinf = re.sub('e','d',"{:.8e}".format(box))
            print("\t%s"%(strrep_d_vinf), end = " ", file = fid)
        print("", file = fid)

        if 'file' in simu.sort:
            if 'wind_field_file' in dir(simu):
                source_field_file = Path(currDir) / Path(simu.wind_field_file)
                if exists(source_field_file):
                    shutil.copy(source_field_file,caseDir)
                else:
                    logging.warning(f"Warning: Wind data file does not exist:\n{source_field_file}")
            print(f'!!', file = fid)
            print(f'!! if sort = file then consider these addtional lines for inflow settings', file = fid)
            print(f'!! filename for wind field: filename with file extension (1)', file = fid)

        # Extenstion of file
            wind_field_file_ext = simu.wind_field_file.split('.')[-1]
            if 'wnd' in wind_field_file_ext:
                wind_field_type = 1
            elif 'bts' in wind_field_file_ext:
                wind_field_type = 2
            else:
                wind_field_type = 1
                logging.warning(f"Warning: wrong file type for wind_field_file! Needed (.wnd or .bts)")
            
            print(f'{simu.wind_field_file}', file = fid)
            print(f'!!', file = fid)
            print(f'!! inflow settings', file = fid)
            print(f'!! grids center (1-3)', file = fid)
            grid_center = [0,0,0]
            if 'grid_center' in dir(simu):
                grid_center = simu.grid_center
            if np.shape(grid_center)[1] != 0:
                grid_center = grid_center.ravel()

            for box in grid_center:
                strrep_grid_center = re.sub('e','d',"{:.8e}".format(box))
                print(f"\t{strrep_grid_center}", end = " ", file = fid)
        else:
            logging.warning(f"Warning: no field: wind_field_file found!")
    os.chdir(simu.currDir)
    print('creating DeSiO-Aero input files')
