import os
import shutil
from pathlib import Path
from os.path import exists
import logging
import numpy as np
import re
def fun_create_batchrunfile(path_DeSiO, path_mtee, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic,flag_init_self_weight,flag_init_rotor_velocity):
    os.chdir(path_fsi)
    with open('DeSiO.bat','w') as fid:
        print(f"set path_DeSiO={path_DeSiO}", file = fid)
        print(f"set path_mtee={path_mtee}", file = fid)
        print("%s"%('REM'), file = fid)
        print(f"set path_fsi={path_fsi}", file = fid)
        print(f"set path_struc={path_struc}", file = fid)
        if flag_init_self_weight == 1:
            print("%s"%('REM Start static self-weight analysis with DeSiO'), file = fid)
            print("cd %path_struc%" + "/" + "%s"%(simu_struct_static.jobname), file = fid)
            print(f'"{path_DeSiO}\DeSiO.exe" | "{path_mtee}\mtee.exe" "%path_struc%\{simu_struct_static.jobname}\output.txt"', file = fid)
            if flag_init_rotor_velocity == 1:
                print("%s"%('REM Copying results files of self weight to dynamic init_rotor_vel simulation'), file = fid)
                print(f"copy%path_struc%\{simu_struct_static.jobname}\{simu_struct_static.jobname}_q.dres' ' '... %path_struc%\{simu_struct_dynamic.jobname}\{simu_struct_static.jobname}_qs.dres", file = fid)
                print(f"copy%path_struc%\{simu_struct_static.jobname}\{simu_struct_static.jobname}_q.dres' ' '... %path_struc%\{simu_struct_dynamic.jobname}\{simu_struct_static.jobname}_v.dres", file = fid)
                print(f"copy%path_struc%\{simu_struct_static.jobname}\{simu_struct_static.jobname}_q.dres' ' '... %path_struc%\{simu_struct_dynamic.jobname}\{simu_struct_static.jobname}_lambda.dres", file = fid)
            else:
                print("%s"%('REM Copying results files of self weight to fsi simulation'), file = fid)
                print(f"copy%path_struc%\{simu_struct_static.jobname}\{simu_struct_static.jobname}_q.dres' ' '... %path_fsi%\{simu_struct_static.jobname}_qs.dres", file = fid)
                print(f"copy%path_struc%\{simu_struct_static.jobname}\{simu_struct_static.jobname}_q.dres' ' '... %path_fsi%\{simu_struct_static.jobname}_v.dres", file = fid)
                print(f"copy%path_struc%\{simu_struct_static.jobname}\{simu_struct_static.jobname}_q.dres' ' '... %path_fsi%\{simu_struct_static.jobname}_lambda.dres", file = fid)

        if flag_init_rotor_velocity == 1:
            print("%s"%('REM Start dynamic init_rotor_vel analysis with DeSiO'), file = fid)
            print(f"cd%path_struc%\{simu_struct_dynamic.jobname}")
            print(f'"{path_DeSiO}\DeSiO.exe" | "{path_mtee}\mtee.exe" "%path_struc%\{simu_struct_dynamic.jobname}\output.txt"', file = fid)
            print("%s"%('REM Copying results files of dynamic struc simulation for init_rotor_vel to fsi simulation'), file = fid)
            print(f"copy%path_struc%\{simu_struct_dynamic.jobname}\{simu_struct_dynamic.jobname}_q.dres' ' '... %path_fsi%\{simu_struct_dynamic.jobname}_qs.dres", file = fid)
            print(f"copy%path_struc%\{simu_struct_dynamic.jobname}\{simu_struct_dynamic.jobname}_q.dres' ' '... %path_fsi%\{simu_struct_dynamic.jobname}_v.dres", file = fid)
            print(f"copy%path_struc%\{simu_struct_dynamic.jobname}\{simu_struct_dynamic.jobname}_q.dres' ' '... %path_fsi%\{simu_struct_dynamic.jobname}_lambda.dres", file = fid)

        print("%s"%('REM Start fsi simulation with DeSiO'), file = fid)
        print(f"cd %path_fsi%", file = fid)
        print("%s"%('REM Boost thread priority'), file = fid)
        print(f"SET desio_exe={path_DeSiO}\DeSiO.exe", file = fid)
        print(f'start "" /REALTIME /B /W  %desio_exe% | "{path_mtee}\mtee.exe" "%path_fsi%\output.txt"', file = fid)
        #print(f'%desio_exe%'  ' | ' '"' {path_mtee} '\mtee.exe' '"' ' "%path_fsi%\output.txt"',file = fid)
        print("%s"%('pause'), file = fid)


#==============================================================================

def fun_create_bashrunfile(path_DeSiO, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic, flag_init_self_weight, flag_init_rotor_velocity):
# write batch file for fsi calculation    
    os.chdir(Path(simu_fsi.currDir)/Path(simu_fsi.strfilename))


    caseDir = Path(simu_fsi.currDir) /  Path(simu_fsi.strfilename) / Path('DeSiO-Aero')
    if not exists(caseDir):
        os.makedirs(caseDir)

    flag1 = Path(simu_fsi.currDir) / Path('my_slurm_set.txt')
    if not exists(flag1):
        # create file with default settings for Luis Cluster at LUH
        with open('my_slurm_set.txt','w') as fid:        
            print("%s"%('#!/bin/bash -l'), file = fid)
            print("%s"%('#SBATCH --ntasks=10'), file = fid)
            print("%s"%('#SBATCH --nodes=1'), file = fid)
            print("%s"%('#SBATCH --ntasks-per-node=1'), file = fid)
            print("%s"%('#SBATCH --mem=20G'), file = fid)
            print("%s"%('#SBATCH --time=12:00:00'), file = fid)
            print("%s"%('#SBATCH --partition=amo'), file = fid)
            print("%s"%('#SBATCH --output output.out'), file = fid)
            print("%s"%('#SBATCH --error error.err'), file = fid)
    else:
        src_path = flag1
        dst_path = Path(simu_fsi.currDir) / Path(simu_fsi.strfilename) / Path(my_slurm_set.txt)
        shutil.copyfile(src_path, dst_path)
        #copyfile([simu_fsi.currDir '\my_slurm_set.txt'],[simu_fsi.currDir '\' simu_fsi.strfilename '\my_slurm_set.txt'])

    
    with open('DeSiO.sh','w') as fid:
        # set path variables
        print("%s"%('#!/bin/bash -l'), file = fid)
        print("%s"%('# set path directories'), file = fid)
        # set directories
        print("%s"%('mypath="$(pwd)"'), file = fid)
        print(f'path_desio="{path_DeSiO}"', file = fid)
        print(f'path_fsi="$mypath/DeSiO" "# to fsi simulation"', file = fid)
        print(f'path_struc="$mypath/DeSiO-Structure" "# to structural simulation"', file = fid)
        # set jobnames
        print(f'#', file = fid)
        print(f'# set jobnames', file = fid)
        print(f'jobname_fsi="{simu_fsi.strfilename}"', file = fid)
        if flag_init_self_weight == 1:
            print(f'jobname_structure_self_weight="{simu_struct_static.jobname}"', file = fid)
        
        if flag_init_rotor_velocity == 1:
            print(f'jobname_structure_init_vel="{simu_struct_dynamic.jobname}"', file = fid)
        
        # set job settings in fsi, e.g. for parallelization. create runDeSiO.sh
        print(f'#', file = fid)
        print(f'# copy my_slurm settings to files', file = fid)
        print(f'cp "$mypath/my_slurm_set.txt" "$path_fsi/runDeSiO.sh"', file = fid)
        print(f'echo "#SBATCH --job-name=$jobname_fsi" >> "$path_fsi/runDeSiO.sh"', file = fid)
        print(f'echo "module load intel/2021a" >> "$path_fsi/runDeSiO.sh"', file = fid)
        print(f'echo "$path_desio" >> "$path_fsi/runDeSiO.sh"', file = fid)
        #print(f'cat "$mypath/version_settings.txt" >> "$path_fsi/runDeSiO.sh"', file = fid)
        print(f'#', file = fid)
        if flag_init_self_weight == 1:
            # set job settings in structure self-weight, e.g. for parallelization. create runDeSiO.sh
            print(f'cp "$mypath/my_slurm_set.txt" "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"', file = fid)
            print(f'echo "#SBATCH --job-name=$jobname_structure_self_weight" >> "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"', file = fid)
            print(f'echo "module load intel/2021a" >> "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"', file = fid)
            print(f'echo "$path_desio" >> "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"', file = fid)
            print(f'#', file = fid)
        
        if flag_init_rotor_velocity == 1:
            # set job settings in structure initial velocity, e.g. for parallelization. create runDeSiO.sh
            print(f'cp "$mypath/my_slurm_set.txt" "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"', file = fid)
            print(f'echo "#SBATCH --job-name=$jobname_structure_init_vel" >> "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"', file = fid)
            print(f'echo "module load intel/2021a" >> "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"', file = fid)
            print(f'echo "$path_desio" >> "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"', file = fid)
            print(f'#', file = fid)
        
        if flag_init_self_weight == 1:
            # changing dir for static self-weight calculation
            print(f'# switching to static self-weight calculations', file = fid)
            print(f'cd "$path_struc/$jobname_structure_self_weight"', file = fid)
            print(f'rm "check.log"', file = fid)
            print(f'#sh runDeSiO.sh', file = fid)
            print(f'sbatch runDeSiO.sh', file = fid)
            print(f'until [ -f "check.log" ]', file = fid)
            print(f'do', file = fid)
            print(f'   sleep 5', file = fid)
            print(f'done', file = fid)
            if flag_init_rotor_velocity == 1:
                print(f'# copying result files to initial rotor velocity', file = fid)
                print('cp "$%s_q.dres" "$path_struc/$jobname_structure_init_vel/$%s_qs.dres"'%('jobname_structure_self_weight','jobname_structure_self_weight'), file = fid)
                print('cp "$%s_v.dres" "$path_struc/$jobname_structure_init_vel/$%s_vs.dres"'%('jobname_structure_self_weight','jobname_structure_self_weight'), file = fid)
                print('cp "$%s_lambda.dres" "$path_struc/$jobname_structure_init_vel/$%s_lambdas.dres"'%('jobname_structure_self_weight','jobname_structure_self_weight'), file = fid)
                print(f'#', file = fid)
            else:
                print(f'# copying result files to initial rotor velocity', file = fid)
                print('cp "$%s_q.dres" "$path_fsi/$%s_qs.dres"'%('jobname_structure_self_weight','jobname_structure_self_weight'), file = fid)
                print('cp "$%s_v.dres" "$path_fsi/$%s_vs.dres"'%('jobname_structure_self_weight','jobname_structure_self_weight'), file = fid)
                print('cp "$%s_lambda.dres" "$path_fsi/$%s_lambdas.dres"'%('jobname_structure_self_weight','jobname_structure_self_weight'), file = fid)
                print(f'#', file = fid)
            
        if flag_init_rotor_velocity == 1:
            # changing dir for dynamic initial velocity calculation
            print(f'# switching to dynamic initial rotor velocity calculations', file = fid)
            print(f'cd "$path_struc/$jobname_structure_init_vel"', file = fid)
            print(f'rm "check.log"', file = fid)
            print(f'#sh runDeSiO.sh', file = fid)
            print(f'sbatch runDeSiO.sh', file = fid)
            print(f'until [ -f "check.log" ]', file = fid)
            print(f'do', file = fid)
            print(f'sleep 5', file = fid)
            print(f'done', file = fid)
            print(f'# copying result files to initial rotor velocity', file = fid)
            print('cp "$%s_q.dres" "$path_fsi/$%s_qs.dres"'%('jobname_structure_init_vel','jobname_structure_init_vel'), file = fid)
            print('cp "$%s_v.dres" "$path_fsi/$%s_vs.dres"'%('jobname_structure_init_vel','jobname_structure_init_vel'), file = fid)
            print('cp "$%s_lambda.dres" "$path_fsi/$%s_lambdas.dres"'%('jobname_structure_init_vel','jobname_structure_init_vel'), file = fid)
            print(f'#', file = fid)
        
        # changing dir for fsi calculation
        print(f'#switching to fsi calculation', file = fid)
        print(f'cd "$path_fsi"', file = fid)
        print(f'sbatch runDeSiO.sh', file = fid)
        print(f'#sh runDeSiO.sh', file = fid)
