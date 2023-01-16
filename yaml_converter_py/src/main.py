from __future__ import annotations


import numpy as np
from copy import deepcopy

from dataclasses import field
import sys
import json
import logging
from utils import yaml_loader, Helpers, skew
import yaml


from numpy import eye, zeros, ones
from numpy.linalg import norm
from math import sin, cos, pi
import os
from pathlib import Path
import glob

current_file_dir = Path(os.path.dirname(os.path.abspath(__file__)))
user_inputs = json.load(open(current_file_dir/Path('user_inputs.json'), 'r'))


file_handler = logging.FileHandler(filename=  Path('app.log'))
stdout_handler = logging.StreamHandler(stream=sys.stdout)
handlers = [file_handler, stdout_handler]

logging.basicConfig(format='%(asctime)s %(levelname)-8s %(message)s',
                    encoding='utf-8', 
                    level=logging.DEBUG,
                    datefmt='%Y-%m-%d %H:%M:%S',
                    handlers=handlers)
logging.info("read_yaml conversion started!")



# reading yml file
data = yaml_loader(user_inputs['yaml_file_path'])
  

# reading user inputs
scaling_blade = user_inputs['scaling_blade']
scaling_tower = user_inputs['scaling_tower']
cos_msl       = np.array(user_inputs['cos_msl'])                                                 # MSL position in global coordinate system

n_yaw         = np.array([user_inputs['n_yaw']])                                                 # yaw rotation axis
n_tilt        = np.array([user_inputs['n_tilt']])                                                # tilt rotation axis

flag_tower       = user_inputs['flag_tower']                                                  # - tower on(1)/off(0)
flag_tower_aero  = user_inputs['flag_tower_aero']                                                  # - tower aero grid on(1)/off(0)
flag_tower_struc = user_inputs['flag_tower_struc']                                                  # - tower aero grid on(1)/off(0)

flag_blade       = user_inputs['flag_blade']                                  # - blades on(1)/off(0)
flag_blade_aero  = user_inputs['flag_blade_aero']                                                  # - blades aero grid on(1)/off(0)
flag_blade_struc = user_inputs['flag_blade_struc']                                                  # - blades structure mesh on(1)/off(0)

flag_foundation       = user_inputs['flag_foundation']                                             # - foundation on(1)/off(0)
flag_foundation_aero  = user_inputs['flag_foundation_aero']                                             # - foundation aero grid on(1)/off(0)
flag_foundation_struc = user_inputs['flag_foundation_struc']                                             # - foundation structure mesh on(1)/off(0)

flag_hub     = user_inputs['flag_hub']                                                      # - hub on(1)/off(0)
flag_nacelle = user_inputs['flag_nacelle']

# setting jobname
jobname = user_inputs['jobname']
if "jobname" in data.keys():
    strfilename = data["jobname"]
else:
    strfilename = jobname
strfilename = strfilename + \
                '_pitch' +  \
                str(data["environment"]["pitch_angle"]) + \
                '_vel' +    \
                str(data["environment"]["vinf"]) + \
                '_nspan' + \
                str(data["components"]["blade"]["DeSiO"]["uvlm"]["M_aero"]) + \
                '_test'
path_DeSiO = str(Path(user_inputs['desio_file_path']))
path_mtee  = str(Path(user_inputs['mtee_file_path']))
if 'path_DeSiO' in data["simulationparamter"].keys():
    path_DeSiO = data["simulationparamter"]["path_DeSiO"]
else:
    logging.warning(f"No path to DeSiO executable given in the yaml file. DesioPath = {path_DeSiO}")


def write_simu_settings():
    dict_data = dict(
                strfilename = strfilename,
                time = data["simulationparamter"]["time"],
                deltat = data["simulationparamter"]["deltat_struc"],
                cutoff = data["simulationparamter"]["cutoff"],
                density = data["environment"]["air_density"],
                i_vinf = data["environment"]["vinf"],
                d_vinf = data["environment"]["winddir"],
                tol = data["simulationparamter"]["tolerance"],
                niter = data["simulationparamter"]["niter"],
                grav = data["simulationparamter"]["gravity"],
                currDir = os.getcwd(),
                grav_vector = data["simulationparamter"]["grav_vector"]
                )
    print(dict_data)
    with open(Path(os.getcwd()) / 'simu_settings.json','w') as outfile:
        json.dump(dict_data,outfile)
write_simu_settings()    
from functions import * 

'''
  Environment Conditions
'''
yaw_angle = 0.0;
pitch_angle = 0.0;
if "environment" in data.keys():
    if "yaw_angle" in data["environment"].keys():
        yaw_angle = data["environment"]["yaw_angle"]
    if "pitch_angle" in data["environment"].keys():
        pitch_angle = data["environment"]["pitch_angle"]
fsi_radius_rbf = 0;
if 'fsi_radius_rbf' in data["simulationparamter"].keys():
    fsi_radius_rbf = data["simulationparamter"]["fsi_radius_rbf"];


'''
    Blades
'''

n_blades = 0;
if 'number_of_blades' in data["assembly"].keys():
    n_blades = data["assembly"]["number_of_blades"];
if "blade" in data["components"].keys() and flag_blade == 1:
#     blade_Aero  = fun_get_uvlm_geometry('blade',data["components"]["blade"]["DeSiO"]["uvlm"],data["airfoils"],scaling_blade);
    blade_Struc = fun_get_beam_model('blade',data["components"]["blade"]["DeSiO"]["beam"],data["materials"],scaling_blade)
    blade_Aero  = fun_get_uvlm_geometry('blade',
                                        data["components"]["blade"]["DeSiO"]["uvlm"],
                                        data["airfoils"],scaling_blade)
    wake_diffusion = 0.0
    if "wake_diffusion" in data["simulationparamter"].keys():
        wake_diffusion = data["simulationparamter"]["wake_diffusion"]

'''
    Tower
'''
tower_top = 0;
tower_bot = 0;
if 'tower' in data["components"].keys() and flag_tower == 1:
	tower_Aero   = fun_get_uvlm_geometry('pipe',data["components"]["tower"]["DeSiO"]["uvlm"],
                                         data["airfoils"],
                                         scaling_tower);    # create 2d mesh for aero grid
	tower_Struc  = fun_get_beam_model('pipe',data["components"]["tower"]["DeSiO"]["beam"],
                                      data["materials"],scaling_tower);      # create structural mesh
	tower_top    = np.max(data["components"]["tower"]["DeSiO"]["beam"]["reference_axis"]["z"]["values"])*scaling_tower;
	tower_bot    = np.min(data["components"]["tower"]["DeSiO"]["beam"]["reference_axis"]["z"]["values"])*scaling_tower;

'''
    Monopile
'''
tp_mass = 0;
if 'monopile' in data["components"].keys():
    monopile_Aero  = fun_get_uvlm_geometry('pipe',
                                            data["components"]["monopile"]["DeSiO"]["uvlm"],
                                            data["airfoils"],
                                            scaling_tower)
    monopile_Struc = fun_get_beam_model('pipe',
                                                data["components"]["monopile"]["DeSiO"]["beam"],
                                                data["materials"],
                                                scaling_tower)
    tp_mass        = data["components"]["monopile"]["transition_piece_mass"]

'''
    Naccele
'''
i1 = np.array([[1,0,0]]) 
i2 = np.array([[0,1,0]]) 
i3 = np.array([[0,0,1]])

uptilt_angle       = 0.0;
overhang           = 0.0;
Twr2Shft           = 0.0;
e_shaft            = np.array([0,0,0]);
nacelle_centerm_tt = np.array([0,0,0])
yaw_mass           = 0.0;
yaw_center         = np.array([[0],[0],[0]])

if 'nacelle' in data["components"].keys() and flag_nacelle == 1:
    uptilt_angle       = data["components"]["nacelle"]["DeSiO"]["drivetrain"]["uptilt"]; # in rad
    overhang           = data["components"]["nacelle"]["DeSiO"]["drivetrain"]["overhang"];
    Twr2Shft           = data["components"]["nacelle"]["DeSiO"]["drivetrain"]["Twr2Shft"];
    yaw_mass           = data["components"]["nacelle"]["DeSiO"]["drivetrain"]["yaw_mass"];
    e_shaft            = -((cos(uptilt_angle)*eye(3) + sin(uptilt_angle)*skew(n_tilt)+(1-cos(uptilt_angle))*n_tilt.T@n_tilt)@i1.T).T
    nacelle_centerm_tt = data["components"]["nacelle"]["DeSiO"]["elastic_properties_mb"]["center_mass"];

shaft_rb1_tt = nacelle_centerm_tt

'''
    Hub
'''
hub_centerm_tt = np.array([0,0,0])
hub_2Apex_tt   = 0.0;
hub_diameter   = 0.0;
if 'hub' in data["components"].keys() and flag_hub == 1:
    hub_diameter   = data["components"]["hub"]["DeSiO"]["diameter"];
    hub_2Apex_tt   = data["components"]["hub"]["DeSiO"]["Hub2Apex"];
    hub_centerm_tt = np.asarray(data["components"]["hub"]["DeSiO"]["elastic_properties_mb"]["center_mass"])
    logging.debug(f"Read from yaml file: hub_centerm_tt:  {hub_centerm_tt}, hub_centerm_tt.shape: {hub_centerm_tt.shape}")

shaft_rb2_tt = hub_centerm_tt;
    
if 'nacelle' in data["components"].keys() and flag_nacelle == 1:
    if 'hub' in data["components"].keys() and flag_hub == 1:
        shaft_rb1_tt   = np.array([0.0, 0.0, Twr2Shft])
        e_shaft_copy  = e_shaft.copy().flatten()
        logging.debug(f'shaft_rb1_tt.shape = {shaft_rb1_tt.shape}, e_shaft_copy.shape = {e_shaft_copy.shape}')
        assert shaft_rb1_tt.ndim == e_shaft_copy.ndim, "shaft_rb1_tt.ndim != e_shaft_copy.ndim"
        shaft_rb2_tt   = shaft_rb1_tt + e_shaft_copy*overhang
        hub_centerm_tt = shaft_rb1_tt + e_shaft_copy*(overhang+hub_2Apex_tt)

cos_tt = cos_msl + np.array([0,0,tower_top])                                   # global coordinates of tower top
cos_tb = cos_msl + np.array([0,0,tower_bot])                                    # global coordinates of tower bottom
logging.debug(f'cost_tt.shape = {cos_tt.shape}, cos_tb.shape = {cos_tb.shape}')
logging.debug(f'hub_centerm_tt.shape = {hub_centerm_tt.shape}')
assert cos_tt.ndim == hub_centerm_tt.ndim, "cos_tt.ndim != hub_centerm_tt.ndim"
cos_hc = cos_tt + hub_centerm_tt;                                      # global coordinates of hub center of mass
cos_na = cos_tt + nacelle_centerm_tt; 


from classes.data.data_structs import *
from classes.data.simulation_settings import *



nrb = -1; nhub = 0; nnac = 0;
att_remove_list = ["cutoff", "density", 
                    "i_vinf", "d_vinf",
                    "sort", "wind_field_file",
                    "grid_center"]
simu_struct = deepcopy(SimulationSettings(rb = [],
                                 pointmass12 = [],
                                 constraints = [],
                                 deltat=data["simulationparamter"]["deltat_struc"],
                                 matbeam = [])
                                 )
simu_aero =  deepcopy(SimulationSettings(rb = [],
                                 pointmass12 = [],
                                 constraints = [],
                                 deltat=data["simulationparamter"]["deltat_aero"],
                                 matbeam = [])
                                 )

for at in att_remove_list:
    delattr(simu_struct, at)

att_remove_list = set(simu_struct.__dict__.keys()) - set(["currDir", "strfilename", "time", "deltat"])
for at in att_remove_list:
    delattr(simu_aero, at)

if 'nacelle' in data["components"].keys() and flag_nacelle :
    simu_struct.rb.append(HubNacceleSettings())
    nrb  = nrb + 1;
    nnac = 0;
    simu_struct.rb[nnac].type = 'nacelle';
    simu_struct.rb[nnac].center_mass = cos_na;
    simu_struct.rb[nnac].D1 = np.array([1,0,0]) 
    simu_struct.rb[nnac].D2 = np.array([0,1,0])
    simu_struct.rb[nnac].D3 = np.array([0,0,1])
    simu_struct.rb[nnac].mass_matrix = zeros((1,10))
    if 'mass_matrix' in data["components"]["nacelle"]["DeSiO"]["elastic_properties_mb"].keys():
        simu_struct.rb[nnac].mass_matrix = data["components"]["nacelle"]["DeSiO"]["elastic_properties_mb"]["mass_matrix"]
    
# =================================================================================================================    
if 'hub' in data["components"].keys() and flag_hub :
    nrb  = nrb + 1;
    if 'nacelle' in data["components"].keys() and flag_nacelle:
        nhub = nnac + 1;
    else:
        nhub = 0;
    simu_struct.rb.append(HubNacceleSettings())

    alpha_hub = 45;
    R_hub     = cos(alpha_hub*pi/180)*eye(3) + sin(alpha_hub*pi/180)*skew(n_yaw)+(1-cos(alpha_hub*pi/180))*n_yaw.T @ n_yaw

    simu_struct.rb[nhub].type = 'hub';
    simu_struct.rb[nhub].center_mass = cos_hc;
    simu_struct.rb[nhub].D1 = (R_hub @ i1.T).T
    simu_struct.rb[nhub].D2 = (R_hub @ i2.T).T
    simu_struct.rb[nhub].D3 = (R_hub @ i3.T).T
    simu_struct.rb[nhub].mass_matrix = zeros((1,10));
    if 'mass_matrix' in data["components"]["hub"]["DeSiO"]["elastic_properties_mb"]:
        simu_struct.rb[nhub].mass_matrix = data["components"]["hub"]["DeSiO"]["elastic_properties_mb"]["mass_matrix"]

simu_fsi = FSISettings()
# simulation setting flags Matlab code line 564 to 594
# flags for pre-calculations
flag_init_self_weight = 0;
if 'flag_init_self_weight' in data["simulationparamter"].keys(): 
    flag_init_self_weight = data["simulationparamter"]["flag_init_self_weight"]

simu_struct.grav = flag_init_self_weight;

flag_init_rotor_velocity = 0;
if 'flag_init_rotor_velocity' in data["simulationparamter"].keys():
    flag_init_rotor_velocity = data["simulationparamter"]["flag_init_rotor_velocity"]


# simulation settings regarding DeSiO executeable
operating_sys = 'windows';
if 'os'  in data["simulationparamter"].keys():
    operating_sys = data["simulationparamter"]["os"]



ifort_version = ''
if 'ifort_version' in data["simulationparamter"].keys():
    ifort_version = data["simulationparamter"]["ifort_version"]

simu_fsi.radius_rbf  = 0;
if 'fsi_radius_rbf' in data["simulationparamter"].keys():
    simu_fsi.radius_rbf  = data["simulationparamter"]["fsi_radius_rbf"]


# external flow field data
simu_aero.sort  = 'constant';
if 'wind_field_file' in data["environment"].keys():
    wind_field_file = data["environment"]["wind_field_file"]
    simu_aero.sort  = 'file';
    simu_aero.wind_field_file = wind_field_file
    if flag_hub:
        simu_aero.grid_center = simu_struct.rb[nhub].center_mass
    else:
        simu_aero.grid_center = []


'''
    Foundation Structure
'''
imat0    = 0; # couner for material
nn12     = 0; # global beam node12 counter
mesh_struc = []
mesh_aero = []
n_struc  = 0
n_surf   = 0
monopile = [Monopile()]
if 'monopile' in data["components"].keys() and flag_foundation :
    # creating structural mesh in DeSiO-Format
    if flag_foundation_struc == 1: 
        monopile[0].grid.struc.init_struc_params()
        X_R = np.hstack([ones((monopile_Struc.M+1,1))*cos_msl[0],                          ones((monopile_Struc.M+1,1))*cos_msl[1],                          ones((monopile_Struc.M+1,1))*cos_msl[2]])
        assert X_R.shape == (monopile_Struc.M+1,3), "X_R shape is wrong"
        monopile[0].grid.struc.X_RE         = monopile_Struc.arr_coordinates[:,:3] + X_R
        assert monopile[0].grid.struc.X_RE.shape[1] == 3, "X_RE shape is wrong"
        monopile[0].grid.struc.D1           = monopile_Struc.arr_coordinates[:,3:6]
        monopile[0].grid.struc.D2           = monopile_Struc.arr_coordinates[:,6:9]
        monopile[0].grid.struc.D3           = monopile_Struc.arr_coordinates[:,9:12]
        monopile[0].grid.struc.M            = monopile_Struc.M;
        monopile[0].grid.struc.connectivity = monopile_Struc.connectivity;

        monopile[0].grid.struc.matBeam.n_gl = nrb + nn12 + np.array(range(monopile_Struc.M + 1)) +1
        monopile[0].grid.struc.matBeam.inz  = np.array(range(monopile_Struc.M))
        monopile[0].grid.struc.matBeam.Cmat = monopile_Struc.arr_stiff_matrix;
        monopile[0].grid.struc.matBeam.mmat = monopile_Struc.arr_mass_matrix;
        monopile[0].grid.struc.matBeam.diss = monopile_Struc.dissipation;

        mesh_struc                             = fun_set_struc_mesh(mesh_struc,monopile)
        nmesh = len(mesh_struc)
        nmesh -= 1
        # last structure which was added named monopile
        mesh_struc[nmesh].strname = 'monopile';
        imat0                                  = mesh_struc[nmesh].imat;
        n_struc                                = n_struc + 1;
        nn12          = monopile[0].grid.struc.matBeam.n_gl[-1] - nrb;

    # creating aero grid in DeSiO-Format
    if flag_foundation_aero : 
        monopile[0].grid.struc.init_aero_params()
        temp = ones(((monopile_Aero.M+1)*(2*monopile_Aero.N+1),1))
        X_R                             = np.hstack([temp * cos_msl[0],                                            temp * cos_msl[1],                                            temp * cos_msl[2]])
        monopile[0].grid.aero.X            = monopile_Aero.X_W + X_R;
        monopile[0].grid.aero.M            = monopile_Aero.M; 
        monopile[0].grid.aero.N = 2*monopile_Aero.N;
        monopile[0].grid.aero.connectivity = monopile_Aero.connectivity_w;
        monopile[0].grid.aero.node_fsi_radius = monopile_Aero.arr_node_fsi_radius

        mesh_aero = fun_set_aero_mesh(mesh_aero,monopile); 
        n_surf    = n_surf + 1
        
        fun_plot_3Dmesh('2D',mesh_aero);


    # set input for FSI
    if flag_foundation_aero  and flag_foundation_struc :
        if not np.any(monopile[0].grid.aero.node_fsi_radius):
            simu_fsi.data.append(['beam',n_struc,n_surf, fsi_radius_rbf])
        else:
            fsi_radius_filename = 'file_fsi_radius_' + str(n_surf) + '_input.txt'
            simu_fsi.data.append( ['beam',n_struc,n_surf, fsi_radius_filename])
            simu_fsi.data_input_node_fsi_radius.append( monopile[0].grid.aero.node_fsi_radius)

    
'''
    Tower
'''
n_struc = len(mesh_struc)
n_surf  = len(mesh_aero)
tower = [deepcopy(Tower())];
if 'tower' in data["components"].keys() and flag_tower == 1:
    # creating structural mesh in DeSiO-Format
    if flag_tower_struc == 1: 
        tower[0].grid.struc.init_struc_params()
        temp = ones((tower_Struc.M+1,1))
        X_R         = np.hstack([temp * cos_msl[0],                                 temp * cos_msl[1],                                 temp * cos_msl[2]])

        tower[0].grid.struc.X_RE         = tower_Struc.arr_coordinates[:,:3] + X_R;
        tower[0].grid.struc.D1           = tower_Struc.arr_coordinates[:,3:6]
        tower[0].grid.struc.D2           = tower_Struc.arr_coordinates[:,6:9]
        tower[0].grid.struc.D3           = tower_Struc.arr_coordinates[:,9:12]
        tower[0].grid.struc.M            = tower_Struc.M;
        tower[0].grid.struc.connectivity = tower_Struc.connectivity;

        tower[0].grid.struc.matBeam.n_gl = nrb + nn12 +                                         np.array(range(tower_Struc.M + 1)) + 1
        tower[0].grid.struc.matBeam.inz  = np.array(range(tower_Struc.M ))
        tower[0].grid.struc.matBeam.Cmat = tower_Struc.arr_stiff_matrix;
        tower[0].grid.struc.matBeam.mmat = tower_Struc.arr_mass_matrix;
        tower[0].grid.struc.matBeam.diss = tower_Struc.dissipation;

        mesh_struc                             = fun_set_struc_mesh(mesh_struc,tower);
        nmesh = len(mesh_struc)
        nmesh -= 1
        mesh_struc[nmesh].strname = 'tower';
        imat0                                  = mesh_struc[nmesh].imat;
        n_struc                                = n_struc + 1;
        nn12          = tower[0].grid.struc.matBeam.n_gl[-1]-nrb;


    # creating aero grid in DeSiO-Format
    if flag_tower_aero == 1: 
        tower[0].grid.struc.init_aero_params()
        temp = ones(((tower_Aero.M+1)*(2*tower_Aero.N+1),1))
        X_R                             = np.hstack([temp * cos_msl[0],                                            temp * cos_msl[1],                                            temp * cos_msl[2]])
        tower[0].grid.aero.X            = tower_Aero.X_W + X_R;
        tower[0].grid.aero.M            = tower_Aero.M; 
        tower[0].grid.aero.N = 2*tower_Aero.N;
        tower[0].grid.aero.connectivity = tower_Aero.connectivity_w;
        tower[0].grid.aero.node_fsi_radius = tower_Aero.arr_node_fsi_radius

        mesh_aero = fun_set_aero_mesh(mesh_aero,tower); 
#         print('mx=',vars(mesh_aero[0])['mx'],',my=',vars(mesh_aero[0])['my'],'nn=', vars(mesh_aero[0])['nn'])
#         print('mx=',vars(mesh_aero[1])['mx'],',my=',vars(mesh_aero[1])['my'],'nn=', vars(mesh_aero[1])['nn'])
        n_surf    = n_surf + 1

        fun_plot_3Dmesh('2D',mesh_aero);

    # set input for FSI
    if flag_tower_aero == 1 and flag_tower_struc == 1:
        if not np.any(tower[0].grid.aero.node_fsi_radius):
            simu_fsi.data.append(['beam',n_struc,n_surf, fsi_radius_rbf])
        else:
            fsi_radius_filename = 'file_fsi_radius_' + str(n_surf) + '_input.txt'
            simu_fsi.data.append( ['beam',n_struc,n_surf, fsi_radius_filename])
            simu_fsi.data_input_node_fsi_radius.append( tower[0].grid.aero.node_fsi_radius)
    
utils_local = Helpers()

'''
    Blades
'''
nbl = 0; imatb = imat0;
blades = [deepcopy(Blade()) for i in range(n_blades)]
blade_obj = []
wake = []
if 'blade' in data["components"] and flag_blade == 1:
    ir1    = -i2 
    ir2    =  i1
    ir3    =  i3

    # Blade root and rest of the blade
    nn = 0; 
    ne = 0; 
    inz_blade_root = [];
    if blade_Aero.nbr != 0:
        # getting indices for extracting blade root from arrays
        arr_inz_airfoil = blade_Aero.arr_inz_airfoil;
        inz = np.nonzero(arr_inz_airfoil<=blade_Aero.nbr -1)[0]
        inz_blade_root = np.array(range(inz[-1]+2)) 

        M_br = len(inz_blade_root)-1;
        N_br = 2*blade_Aero.N;

        nn_br = (M_br+1)*(blade_Aero.N+1);
        ne_br = M_br*blade_Aero.N;
        for i  in range(n_blades):
            n_struc = len(mesh_struc);
            n_surf  = len(mesh_aero);

            # position angle of blade i in rotor
            alpha_blades = (i)*2*np.pi/n_blades*180/np.pi;

            # positioning blade in rotor
            blades[i].grid = deepcopy(fun_blades2rotor(blade_Aero, 
                                                blade_Struc, 
                                                ir3, 
                                                -pitch_angle, 
                                                ir2, 
                                                alpha_blades, 
                                                n_tilt, 
                                                uptilt_angle*180/pi, 
                                                [0,0,hub_diameter/2])); 


            
            # creating structural mesh in DeSiO-Format
            logging.debug(f'cos_hc = {cos_hc}, cos_hc.shape = {cos_hc.shape}')
            if flag_blade_struc == 1:
                nbl = blade_Struc.M+1;
                temp = ones((blade_Struc.M+1,1))
                X_R = np.hstack([temp * cos_hc[0],
                                temp * cos_hc[1],
                                temp * cos_hc[2]])

                blades[i].grid.struc.X_RE = blades[i].grid.struc.X_RE + X_R;
                blades[i].grid.struc.M = blade_Struc.M;
                blades[i].grid.struc.connectivity = blade_Struc.connectivity;
#                 print(f'nrb = {nrb}')
#                 print(f'nn12 = {nn12}')
#                 print(f'nbl = {nbl}')
#                 print(f'np.array(range(nbl)) = {np.array(range(nbl))}')
                blades[i].grid.struc.matBeam.n_gl = (nrb+1) + nn12 +                                                     (i)*nbl + np.array(range(nbl))
                # logging.debug(f'blades[{0}].grid.struc.matBeam.n_gl = {blades[0].grid.struc.matBeam.n_gl}')
                logging.debug(f'blades[{i}].grid.struc.matBeam.n_gl = {blades[i].grid.struc.matBeam.n_gl}')
                blades[i].grid.struc.matBeam.inz  = np.array(range(blade_Struc.M))
                blades[i].grid.struc.matBeam.Cmat = blade_Struc.arr_stiff_matrix
                blades[i].grid.struc.matBeam.mmat = blade_Struc.arr_mass_matrix
                blades[i].grid.struc.matBeam.diss = blade_Struc.dissipation

                mesh_struc = fun_set_struc_mesh(mesh_struc,[blades[i]]);
                nmesh = len(mesh_struc)
                nmesh -= 1
                mesh_struc[nmesh].strname = 'blade ' + str(i)
                mesh_struc[nmesh].imat = imatb;

                n_struc = n_struc + 1;
                nbl = n_blades*nbl;

            
            # creating aero grid in DeSiO-Format
            if flag_blade_aero == 1:

                temp = ones(((blade_Aero.M+1)*(blade_Aero.N+1),1))
                X_R = utils_local.hstack(temp, cos_hc)
                
                X_C = blades[i].grid.aero.X_C  + X_R

                temp = ones(((blade_Aero.M+1)*(2*blade_Aero.N+1),1))
                X_R = utils_local.hstack(temp, cos_hc)
                X_W = blades[i].grid.aero.X_W  + X_R;

                # Blade root and rest of the blade
                if len(blade_obj)<2:
                    logging.info(f'i = {i},len(blade_obj) = {len(blade_obj)}')
                    # there are only two blades
                    blade_obj.append(deepcopy(Blade2()))
                if blade_Aero.nbr != 0:
                    #print(f'blade_obj[0].grid.aero.N = {blade_obj[0].grid.aero.N}')
                    # extracting blade root from arrays
                    blade_obj[0].grid.aero.N = N_br;
                    
                    blade_obj[0].grid.aero.M = M_br;
                    blade_obj[0].grid.aero.X = X_W[:(N_br+1)*(M_br+1),:]
                    blade_obj[0].grid.aero.connectivity =                                 blade_Aero.connectivity_w[:M_br*N_br,:]
                    # TODO: check if this is correct
                    blade_obj[0].grid.aero.node_fsi_radius = blade_Aero.arr_node_fsi_radius_w[:(N_br+1)*(M_br+1)]

                    
                    if len(blade_obj)<2:
                        logging.info(f'i = {i},len(blade_obj) = {len(blade_obj)}')
                        # there are only two blades
                        blade_obj.append(deepcopy(Blade2()))
                    # extracting rest of blade from arrays
                    
                    blade_obj[1].grid.aero.N = blade_Aero.N;
                    blade_obj[1].grid.aero.M = blade_Aero.M-M_br;
#                     print('X_C=',X_C)
                    blade_obj[1].grid.aero.X = X_C[(nn_br-blade_Aero.N-1):,:]
                    blade_obj[1].grid.aero.connectivity = blade_Aero.connectivity_c[ne_br:,:]-(nn_br-blade_Aero.N)+1;
                    blade_obj[1].grid.aero.node_fsi_radius = blade_Aero.arr_node_fsi_radius_w[:(N_br+1)*(M_br+1)]
                else:
                    blade_obj[0].grid.aero.N = blade_Aero.N;
                    blade_obj[0].grid.aero.M = blade_Aero.M;
                    blade_obj[0].grid.aero.X = X_C;
                    blade_obj[0].grid.aero.connectivity = blade_Aero.connectivity_c;
                    blade_obj[0].grid.aero.node_fsi_radius = blade_Aero.arr_node_fsi_radius


                
                # creating mesh in DeSiO-Format
                # print('mx=',vars(mesh_aero[0])['mx'],',my=',vars(mesh_aero[0])['my'],'nn=', vars(mesh_aero[0])['nn'])
                # print('mx=',vars(mesh_aero[1])['mx'],',my=',vars(mesh_aero[1])['my'],'nn=', vars(mesh_aero[1])['nn'])
                
                mesh_aero = fun_set_aero_mesh(mesh_aero,blade_obj);
                print(f'len(mesh_aero) = {len(mesh_aero)}')
                for obj in mesh_aero:
                    if hasattr(obj, 'imat'):
                        # print(f'obj.__dict__ = {obj.__dict__}')
                        if 'imat' in obj.__dict__.keys():
                            delattr(obj,'imat')
                            # del obj.__dict__['imat']
                    if hasattr(obj, 'strname'):
                        if 'strname' in obj.__dict__.keys():
                            delattr(obj,'strname')
                            # del obj.__dict__['strname']
                
                # start troubleshooting from line 396 before that it's correct

                # logging.debug('Inputs to fun_blade_wake')
                # logging.debug(f'n_surf + len(blade_obj)= {n_surf + len(blade_obj)}')
                # logging.debug(f'data["simulationparamter"]["nwakerows"] = {data["simulationparamter"]["nwakerows"]}')
                # logging.debug(f'wake = {wake}')
                # logging.debug(f'mesh_aero[n_surf + len(blade_obj)-1] = {mesh_aero[n_surf + len(blade_obj)-1]}')
                # logging.debug('Inputs to fun_blade_wake end')
                wake = fun_blade_wake(wake,
                                    n_surf + len(blade_obj),
                                    mesh_aero[n_surf + len(blade_obj)-1],
                                    data["simulationparamter"]["nwakerows"],
                                    data["simulationparamter"]["nwakerows"],
                                    wake_diffusion);
                

            
            # setting input for FSI
            if flag_blade_aero == 1 and flag_blade_struc == 1:
                logging.debug('Setting input for FSI')
                logging.debug(f'len(blade_obj) = {len(blade_obj)}')
                for j in range(len(blade_obj)):
                    if not np.any(blade_obj[j].grid.aero.node_fsi_radius):
                        simu_fsi.data.append(['beam',n_struc,n_surf+j+1, fsi_radius_rbf])
                    else:
                        fsi_radius_filename = 'file_fsi_radius_' + str(n_surf+j) + '_input.txt'
                        simu_fsi.data.append( ['beam',n_struc,n_surf+j, fsi_radius_filename])
                        simu_fsi.data_input_node_fsi_radius.append( blade_obj[j].grid.aero.node_fsi_radius)
                    

logging.debug(f'len(simu_fsi.data) = {len(simu_fsi.data)}')
logging.debug(f'simu_fsi.data = {simu_fsi.data}')


'''
    Additional Masses
'''
logging.info("Additional masses sections")
logging.info("*"*50)
nam = 0;
if monopile and flag_foundation_struc == 1:
	simu_struct.pointmass12.append(deepcopy(simu_struct.pointmass))
	simu_struct.pointmass12[nam].mass = tp_mass;
	simu_struct.pointmass12[nam].node = monopile[0].grid.struc.matBeam.n_gl[-1];
	nam = nam + 1;

# additional mass for yaw bearing at tower top
if 'nacelle' in data["components"].keys() and flag_nacelle == 1:
	if tower and flag_tower_struc == 1:
		simu_struct.pointmass12.append(deepcopy(simu_struct.pointmass))
		simu_struct.pointmass12[nam].mass = yaw_mass;
		simu_struct.pointmass12[nam].node = tower[0].grid.struc.matBeam.n_gl[-1];
		nam = nam + 1;
logging.debug(f'len(simu_struct.pointmass12) = {len(simu_struct.pointmass12)}')
logging.debug(f'simu_struct.pointmass12 = {simu_struct.pointmass12}')


'''
    Constraints
'''
logging.info("Additional masses sections end")
logging.info("*"*50)
logging.info("Constraints sections")
logging.info("*"*50)
simu_struct.constraints = [];
remo_internal_const = [];
nco = 0;
if tower and flag_tower_struc :
    if monopile and flag_foundation_struc :
        
        simu_struct.constraints = utils_local.try_append(simu_struct.constraints, 
                                                    nco,
                                                    deepcopy(simu_struct.constraint), 
                                                    )
        # simu_struct.constraints.append(deepcopy(simu_struct.constraint))
        simu_struct.constraints[nco].type  = 'rigidsupport'
        simu_struct.constraints[nco].nodes = np.array([monopile[0].grid.struc.matBeam.n_gl[0], 0])
        simu_struct.constraints[nco].phi1 = []
        simu_struct.constraints[nco].phi2 = []
        simu_struct.constraints[nco].dir = []
        remo_internal_const = utils_local.try_append(remo_internal_const,
                                                nco,
                                                monopile[0].grid.struc.matBeam.n_gl[0]
                                                )
        nco = nco + 1;
        
        
        
        simu_struct.constraints = utils_local.try_append(simu_struct.constraints, 
                                                    nco,
                                                    deepcopy(simu_struct.constraint), 
                                                    )
        
        simu_struct.constraints[nco].type  = 'rigidconnection';
        simu_struct.constraints[nco].nodes = np.array([monopile[0].grid.struc.matBeam.n_gl[-1], tower[0].grid.struc.matBeam.n_gl[0]])
        simu_struct.constraints[nco].phi1 = []
        simu_struct.constraints[nco].phi2 = []
        simu_struct.constraints[nco].dir = []
        remo_internal_const = utils_local.try_append(remo_internal_const,
                                                nco,
                                                monopile[0].grid.struc.matBeam.n_gl[-1]
                                                )
        nco = nco + 1;
        # remo_internal_const.append(monopile[0].grid.struc.matBeam.n_gl[-1] )
    else:
        
        simu_struct.constraints = utils_local.try_append(simu_struct.constraints, 
                                                    nco,
                                                    deepcopy(simu_struct.constraint), 
                                                    )
        simu_struct.constraints[nco].type  = 'rigidsupport';
        simu_struct.constraints[nco].nodes = np.array([tower[0].grid.struc.matBeam.n_gl[0], 0])
        simu_struct.constraints[nco].phi1 = []
        simu_struct.constraints[nco].phi2 = []
        simu_struct.constraints[nco].dir = []
        remo_internal_const = utils_local.try_append(remo_internal_const,
                                                nco,    
                                                tower[0].grid.struc.matBeam.n_gl[0]
                                                )
        nco = nco + 1;
        # remo_internal_const.append( tower[0].grid.struc.matBeam.n_gl[0])
        
    if 'nacelle' in data["components"] and flag_nacelle == 1:
        
        simu_struct.constraints = utils_local.try_append(simu_struct.constraints, 
                                                    nco,
                                                    deepcopy(simu_struct.constraint), 
                                                    )
        simu_struct.constraints[nco].type  = 'rigidconnection';
        simu_struct.constraints[nco].nodes = np.array([tower[0].grid.struc.matBeam.n_gl[-1], nnac])
        
        #             phi1 = nacelle_centerm_tt;
        simu_struct.constraints[nco].phi1 = []
        simu_struct.constraints[nco].phi2 = []
        simu_struct.constraints[nco].dir = []
        remo_internal_const = utils_local.try_append(remo_internal_const,
                                                nco,
                                                tower[0].grid.struc.matBeam.n_gl[-1]
                                                )
        nco = nco + 1
        # remo_internal_const.append( tower[0].grid.struc.matBeam.n_gl[-1])
logging.debug(f'len(simu_struct.constraints) = {len(simu_struct.constraints)}')
logging.debug(f'simu_struct.constraints = {simu_struct.constraints}')
logging.debug(f'remo_internal_const = {remo_internal_const}')
dir_lo = []
if "hub" in data["components"] and flag_hub == 1:
    if "nacelle" in data["components"] and flag_nacelle == 1:
        
        nhub = nnac + 1;
        simu_struct.constraints = utils_local.try_append(simu_struct.constraints, 
                                                    nco,
                                                    deepcopy(simu_struct.constraint), 
                                                    )
        simu_struct.constraints[nco].type  = 'revolutejoint';

        simu_struct.constraints[nco].nodes = np.array([nnac, nhub])
        # simu_struct.constraints[nco].dir  = np.array([e_shaft @ simu_struct.rb[nnac].D1.T, 
        #                                     e_shaft @ simu_struct.rb[nnac].D2.T, 
        #                                     e_shaft @ simu_struct.rb[nnac].D3.T])

        simu_struct.constraints[nco].dir  = []
        #          phi1 = [cos_tt + shaft_rb2_tt] - simu_struct.rb[nnac].center_mass ;
        if cos_tt.ndim == 1:
            cos_tt = np.expand_dims(cos_tt, axis=0)
        phi1 = cos_tt + shaft_rb1_tt - simu_struct.rb[nnac].center_mass ;
        simu_struct.constraints[nco].phi1 = np.array([phi1 @ simu_struct.rb[nnac].D1.T, 
                                                    phi1 @ simu_struct.rb[nnac].D2.T, 
                                                    phi1 @ simu_struct.rb[nnac].D3.T])

        logging.debug(f'phi1 = {phi1}')
        phi2 = cos_tt + shaft_rb2_tt - simu_struct.rb[nhub].center_mass;
        logging.debug(f'cos_tt = {cos_tt}')
        logging.debug(f'shaft_rb2_tt = {shaft_rb2_tt}')
        logging.debug(f'simu_struct.rb[nhub].center_mass = {simu_struct.rb[nhub].center_mass}')
        simu_struct.constraints[nco].phi2 = np.array([phi2 @ simu_struct.rb[nhub].D1.T, 
                                                    phi2 @ simu_struct.rb[nhub].D2.T, 
                                                    phi2 @ simu_struct.rb[nhub].D3.T])
        logging.debug(f'phi2 = {phi2}')
        logging.debug(f'simu_struct.constraints[{nco}] = {simu_struct.constraints[nco]}')
        # rotor rotation axis in global cos and local hub cos:
        #TODO: check if this is correct
        dir_gl = simu_struct.rb[nhub].center_mass -                     (simu_struct.rb[nnac].center_mass + phi1)
        dir_gl = dir_gl/norm(dir_gl);
        dir_lo = np.array([ dir_gl @ simu_struct.rb[nhub].D1.T, 
                    dir_gl @ simu_struct.rb[nhub].D2.T, 
                    dir_gl @ simu_struct.rb[nhub].D3.T])
        nco  = nco + 1;

if blades and flag_blade_struc == 1:
    if "hub" in data["components"].keys() and flag_hub == 1:
        nco0 = nco;
        j_temp =1
        logging.debug(f'len(blades) = {len(blades)}')
        for i in range(n_blades):
            
            nhub = nnac + 1;
            simu_struct.constraints = utils_local.try_append(simu_struct.constraints,
                                                        nco,
                                                        deepcopy(simu_struct.constraint),
                                                        )
            simu_struct.constraints[nco].type  = 'rigidconnection';
            simu_struct.constraints[nco].nodes = np.array([nhub, blades[i].grid.struc.matBeam.n_gl[0]])
            # TODO: check if this is correct
            simu_struct.constraints[nco].phi1  = np.array([0,0,0])
            logging.debug(f'simu_struct.constraints[{nco}].nodes = {simu_struct.constraints[nco].nodes}')
            
#                 phi1 = blades[i].grid.struc.X_RE[0,:] - simu_struct.rb[nhub].center_mass;
            phi1 = np.array([[0,0,0]]);
            simu_struct.constraints[nco].phi1 = np.array([phi1 @ simu_struct.rb[nhub].D1.T, 
                                                phi1 @ simu_struct.rb[nhub].D2.T, 
                                                phi1 @ simu_struct.rb[nhub].D3.T])
            simu_struct.constraints[nco].phi2 = []
            simu_struct.constraints[nco].dir = []

            remo_internal_const = utils_local.try_append(remo_internal_const,
                                                nco,
                                                blades[i].grid.struc.matBeam.n_gl[0]
                                                )
            logging.debug(f'simu_struct.constraints[{nco}] = {simu_struct.constraints[nco]}')
            nco  = nco0 + j_temp ;
            j_temp += 1
        nco += 1
        logging.debug(f'remo_internal_const = {remo_internal_const}')
        if (not tower) or (flag_tower_struc == 0):
            # TODO: check if this is correct
            simu_struct.constraints = utils_local.try_append(simu_struct.constraints,
                                                        nco,
                                                        deepcopy(simu_struct.constraint),
                                                        )
            simu_struct.constraints[nco].type  = 'sphericalsupport';
            simu_struct.constraints[nco].nodes = np.array([nhub, 0])
            simu_struct.constraints[nco].phi1 = []
            simu_struct.constraints[nco].phi2 = []
            simu_struct.constraints[nco].dir = []
            nco = nco + 1;
            simu_struct.constraints = utils_local.try_append(simu_struct.constraints,
                                                        nco,
                                                        deepcopy(simu_struct.constraint),
                                                        )
            simu_struct.constraints[nco].type  = 'rotation_global';
            simu_struct.constraints[nco].nodes = np.array([nhub, 0])
            simu_struct.constraints[nco].dir   = np.array([[0.00,1.00,0.00]])
            simu_struct.constraints[nco].phi1 = []
            simu_struct.constraints[nco].phi2 = []
            nco = nco + 1;
            simu_struct.constraints = utils_local.try_append(simu_struct.constraints,
                                                        nco,
                                                        deepcopy(simu_struct.constraint),
                                                        )
            simu_struct.constraints[nco].type  = 'rotation_global';
            simu_struct.constraints[nco].nodes = np.array([nhub, 0])
            simu_struct.constraints[nco].dir   = np.array([[0.00,0.00,1.00]])
            simu_struct.constraints[nco].phi1 = []
            simu_struct.constraints[nco].phi2 = []
            nco = nco + 1;
logging.debug(f'len(simu_struct.constraints) = {len(simu_struct.constraints)}')
logging.debug(f'simu_struct.constraints = {simu_struct.constraints}')


if 'nbl' in locals():
    nn = nrb + nn12 + nbl + 1
else:
    nn = nrb + nn12 + 1 
nodes = np.array(range(nn))
logging.debug(f'node: {nodes}')
ai = len(simu_struct.constraints)
for i in range(nn):
    if all(remo_internal_const-nodes[i]):
        simu_struct.constraints = utils_local.try_append(simu_struct.constraints,
                                                    ai,
                                                    deepcopy(simu_struct.constraint),
                                                    )
        simu_struct.constraints[ai].type  = 'internal';
        simu_struct.constraints[ai].nodes = np.array([nodes[i], 0])
        simu_struct.constraints[ai].phi1 = []
        simu_struct.constraints[ai].phi2 = []
        simu_struct.constraints[ai].dir = []
        logging.debug(f'simu_struct.constraints[{ai}] = {simu_struct.constraints[ai]}')
        
        ai = ai + 1;

imat = 0
if monopile and flag_foundation_struc :
    for i in range(monopile[0].grid.struc.matBeam.Cmat.shape[0]):
        simu_struct.matbeam = utils_local.try_append(simu_struct.matbeam,
                                                i,
                                                deepcopy(simu_struct.materialbeam),
                                                )
        simu_struct.matbeam[i].cmat    = monopile[0].grid.struc.matBeam.Cmat[i,:];
        simu_struct.matbeam[i].mmat    = monopile[0].grid.struc.matBeam.mmat[i,:];
        simu_struct.matbeam[i].diss    = monopile[0].grid.struc.matBeam.diss;
        simu_struct.matbeam[i].strname = 'monopile';
    imat = i;
    logging.debug(f'simu_struct.matbeam = {simu_struct.matbeam}')
    logging.debug(f'len(simu_struct.matbeam) = {len(simu_struct.matbeam)}')
    logging.debug(f'len(simu_struct.matbeam[i].cmat) = {len(simu_struct.matbeam[0].cmat)}')
    logging.debug(f'len(simu_struct.matbeam[i].mmat) = {len(simu_struct.matbeam[0].mmat)}')
    logging.debug(f'len(simu_struct.matbeam[i].diss) = {len(simu_struct.matbeam[0].diss)}')
    logging.debug(f'simu_struct.matbeam[i].strname = {simu_struct.matbeam[0].strname}')


if flag_foundation_struc:
    imat += 1
# TODO: check if this is correct
if tower and flag_tower_struc == 1:
    for i in range(tower[0].grid.struc.matBeam.Cmat.shape[0]):
        simu_struct.matbeam = utils_local.try_append(simu_struct.matbeam,
                                                i+imat,
                                                deepcopy(simu_struct.materialbeam),
                                                )
        simu_struct.matbeam[i+imat].cmat    = tower[0].grid.struc.matBeam.Cmat[i,:];
        simu_struct.matbeam[i+imat].mmat    = tower[0].grid.struc.matBeam.mmat[i,:];
        simu_struct.matbeam[i+imat].diss    = tower[0].grid.struc.matBeam.diss;
        simu_struct.matbeam[i+imat].strname = 'tower';
    imat = i + imat;
    logging.debug(f'simu_struct.matbeam = {simu_struct.matbeam}')
    logging.debug(f'len(simu_struct.matbeam) = {len(simu_struct.matbeam)}')
    logging.debug(f'len(simu_struct.matbeam[i].cmat) = {len(simu_struct.matbeam[0].cmat)}')
    logging.debug(f'len(simu_struct.matbeam[i].mmat) = {len(simu_struct.matbeam[0].mmat)}')
    logging.debug(f'len(simu_struct.matbeam[i].diss) = {len(simu_struct.matbeam[0].diss)}')
    logging.debug(f'simu_struct.matbeam[i].strname = {simu_struct.matbeam[0].strname}')


if flag_tower_struc:
    imat += 1
# TODO: check if this is correct
if blades and flag_blade_struc == 1:      
    for i in range(blades[0].grid.struc.matBeam.Cmat.shape[0]):
        simu_struct.matbeam = utils_local.try_append(simu_struct.matbeam,
                                                i+imat,
                                                deepcopy(simu_struct.materialbeam),
                                                )
        simu_struct.matbeam[i+imat].cmat = blades[0].grid.struc.matBeam.Cmat[i,:];
        simu_struct.matbeam[i+imat].mmat = blades[0].grid.struc.matBeam.mmat[i,:];
        simu_struct.matbeam[i+imat].diss = blades[0].grid.struc.matBeam.diss;
        simu_struct.matbeam[i+imat].strname = 'blade';
    logging.debug(f'simu_struct.matbeam = {simu_struct.matbeam}')
    logging.debug(f'len(simu_struct.matbeam) = {len(simu_struct.matbeam)}')
    logging.debug(f'len(simu_struct.matbeam[i].cmat) = {len(simu_struct.matbeam[0].cmat)}')
    logging.debug(f'len(simu_struct.matbeam[i].mmat) = {len(simu_struct.matbeam[0].mmat)}')



'''
    WRITING DESIO-INPUT FILES
'''
#==============================================================================
# creating new objects
simu_struct_static = deepcopy(SimulationSetting_speical())
simu_struct_static.copy(simu_struct)

simu_struct_dynamic = deepcopy(SimulationSetting_speical())
simu_struct_dynamic.add_dynamic_attributes()

simu_struct_modal = deepcopy(simu_struct_static)

if simu_aero and np.shape(mesh_aero)[0] !=0:
    fun_writeDeSiOAeroInput(simu_aero,mesh_aero,wake)
if simu_struct and np.shape(mesh_struc)[0] !=0:
    # initial DeSiO-Structure input
    path_struc = Path(simu_struct.currDir) / Path(simu_struct.strfilename) / Path('DeSiO-Structure')
    simu_struct.jobname  = 'initial_data';
    fun_writeDeSiOStructureInput(simu_struct,mesh_struc)
    # writing structural input file for static pre-analysis for self-weight
    simu_struct_static.type     = 'static';
    simu_struct_static.jobname  = 'self_weight';
    simu_struct_static.settings = np.array([1.0, 1.0e-1, 1.0e-6, 50, 1])
    # adding constraint to fix rotor-axis
    temp_idx = len(simu_struct_static.constraints)
    simu_struct_static.constraints = utils_local.try_append(simu_struct_static.constraints, 
                                                temp_idx,
                                                deepcopy(simu_struct_static.constraint), 
                                                )
    simu_struct_static.constraints[temp_idx].type  = 'rotation_local';
    simu_struct_static.constraints[temp_idx].nodes = np.array([nhub, 0]);
    simu_struct_static.constraints[temp_idx].dir   = dir_lo
    simu_struct_static.constraints[temp_idx].phi1   = []
    simu_struct_static.constraints[temp_idx].phi2   = []
    fun_writeDeSiOStructureInput(simu_struct_static,mesh_struc);

    # writing structural input file for dynamic pre-analysis for initial rotor speed

    if flag_init_rotor_velocity == 1:
        simu_struct_dynamic.type     = 'dynamic';
        simu_struct_dynamic.jobname  = 'init_rotor_speed';
        if flag_init_self_weight == 1:
            simu_struct_dynamic.settings = np.array([10.0 ,1.0e-1, 1.0e-6, 50, 1])
        else:
            simu_struct_dynamic.settings = np.array([10.0 ,1.0e-1, 1.0e-6, 50, 0])

        # asign simulation parameter
        temp_idx = len(simu_struct_dynamic.constraints)
        simu_struct_dynamic.constraints = utils_local.try_append(simu_struct_dynamic.constraints,
                                                    temp_idx,
                                                    deepcopy(simu_struct_dynamic.constraint),
                                                    )
        simu_struct_dynamic.constraints[temp_idx].type  = 'angularvelocity_local';
        simu_struct_dynamic.constraints[temp_idx].nodes   = np.array([nhub, 0])
        simu_struct_dynamic.constraints[temp_idx].dir     = dir_lo
        simu_struct_dynamic.constraints[temp_idx].phi1    = []
        simu_struct_dynamic.constraints[temp_idx].phi2    = []
        init_rotor_velocity = 0;
        if 'init_rotor_velocity' in data["environment"].keys():
            init_rotor_velocity = -data["environment"]["init_rotor_velocity"] * 2 * pi/60;

        simu_struct_dynamic.boundary12.sort         = 'linearC'
        simu_struct_dynamic.boundary12.intensity    = init_rotor_velocity;
        simu_struct_dynamic.boundary12.time         = 10.0;
        simu_struct_dynamic.boundary12.constraintID = len(simu_struct_dynamic.constraints);
        simu_struct_dynamic.boundary12.file         = 'none';
        # create files but statement to ask if self weight is on or not
        if flag_init_self_weight == 1:
            simu_struct_dynamic.prevjobname = simu_struct_static.jobname;

        fun_writeDeSiOStructureInput(simu_struct_dynamic,mesh_struc);


    # writing structural input file for modal analysis
    simu_struct_modal.type     = 'modal';
    simu_struct_modal.jobname  = 'modalanalysis';
    simu_struct_modal.settings = np.array([10, 1.0e-6, 0.0, 0.0, 0])
    fun_writeDeSiOStructureInput(simu_struct_modal,mesh_struc);
    
if simu_fsi:
    fun_writeDeSiOFSIInput(simu_fsi,simu_struct,simu_aero)
    if simu_aero and simu_struct:
        path_fsi = Path(simu_fsi.currDir + "/"  +  simu_fsi.strfilename + '/DeSiO')
        path_aero = Path(simu_aero.currDir + "/"  +  simu_aero.strfilename + '/DeSiO-Aero')
        path_struc = Path(simu_struct.currDir + "/"  +  simu_struct.strfilename + '/DeSiO-Structure')

        filePattern = path_aero  / Path('*.txt')
        print(filePattern)
        for file in glob.glob(str(filePattern)):
            file = Path(file)
            if not 'simulationinput' in file.name:
                src_path = path_aero / file.name
                dst_path = path_fsi / file.name
                shutil.copyfile(src_path, dst_path)
        filePattern = path_struc / Path('initial_data') / Path('*.txt')
        for file in glob.glob(str(filePattern)):
            file = Path(file)
            if not 'simulationinput' in file.name:
                src_path = path_struc / Path('initial_data')/ file.name
                dst_path = path_fsi / file.name
                shutil.copyfile(src_path, dst_path)

# create fsi input files
if flag_init_rotor_velocity == 1: # set jobname for fsi calculation
    simu_fsi.structprevjobname = simu_struct_dynamic.jobname;
    simu_fsi.aeroprevjobname  = 'none';
elif flag_init_rotor_velocity == 0 and flag_init_self_weight == 1:
    simu_fsi.structprevjobname = simu_struct_static.jobname;
    simu_fsi.aeroprevjobname  = 'none';

fun_writeDeSiOFSIInput(simu_fsi,simu_struct,simu_aero);

#  create batch or bash files for os windows/linux
if 'windows' in operating_sys:
    fun_create_batchrunfile(path_DeSiO, path_mtee, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic, flag_init_self_weight, flag_init_rotor_velocity);
else:
    fun_create_bashrunfile(path_DeSiO, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic, flag_init_self_weight, flag_init_rotor_velocity);
