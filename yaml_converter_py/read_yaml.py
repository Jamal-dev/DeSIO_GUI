import argparse
import numpy as np
import json
from pathlib import Path
import os
import subprocess

current_file_dir = os.path.dirname(os.path.abspath(__file__))
class Store_as_array(argparse._StoreAction):
    def __call__(self, parser, namespace, values, option_string=None):
        values = np.array(values)
        return super().__call__(parser, namespace, values, option_string)
def check_args():
    parser = argparse.ArgumentParser(prog='This is a program to read the yaml file and create desired inputs files for DeSIO')
    parser.add_argument("--yaml_file_path", type=str, help="Enter the path for yaml file",
                        nargs='?', default=str(Path(current_file_dir)/Path('IEA-15-240-RWT_new.yaml')) )
    parser.add_argument("--scaling_blade", type=float, help="Enter the scaling for the blade",
                        nargs='?', default=1 )
    parser.add_argument("--scaling_tower", type=float, help="Enter the scaling for the tower",
                        nargs='?', default=1 )
    parser.add_argument('--cos_msl', action=Store_as_array, type=float, 
                        nargs='?',help = "Enter the position of the MSL ex: --cos_ml 0 0 0",
                        default=[0.0, 0.0, 0.0] )
    parser.add_argument('--n_yaw', action=Store_as_array, type=float, 
                        nargs='+',help = "Enter the rotation axis of yaw ex: --n_yaw 0 0 1",
                        default=[0.0, 0.0, 1.0] )
    parser.add_argument('--n_tilt', action=Store_as_array, type=float, 
                        nargs='?',help = "Enter the rotation axis of tilt ex: --n_tilt 0 1 0",
                        default=[0.0, 1.0, 0.0] )
    parser.add_argument('--flag_tower', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate tower file',
                        default=True )
    parser.add_argument('--flag_blade', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate blade file',
                        default=True )
    parser.add_argument('--flag_tower_aero', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate tower aero file',
                        default=True )
    parser.add_argument('--flag_tower_struc', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate tower struc file',
                        default=True )
    parser.add_argument('--flag_blade_aero', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate blade aero file',
                        default=True )
    parser.add_argument('--flag_blade_struc', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate blade struc file',
                        default=True )
    parser.add_argument('flag_foundation', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate foundation file',
                        default=True )
    parser.add_argument('--flag_foundation_aero', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate foundation aero file',
                        default=True )
    parser.add_argument('--flag_foundation_struc', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate foundation struc file',
                        default=True )
    parser.add_argument('--flag_hub', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate hub file',
                        default=True )
    parser.add_argument('--flag_nacelle', action=argparse.BooleanOptionalAction,help = '--no-flag will not generate file and --flag will generate nacelle file',
                        default=True )
    parser.add_argument("--jobname", type=str, help="Enter the name of the job if not specified in yaml file",
                        nargs='?', default='mesh_tests_blade' )
    

    args = parser.parse_args()
    
    return args

if __name__ == '__main__':
    user_inputs = check_args()
    with open(Path('src/user_inputs.json'), 'w') as f:
        json.dump(vars(user_inputs), f)
    
    subprocess.run(['python', 'src/main.py'])
