from typing import List
from copy import deepcopy
from classes.data.data_structs import MeshStruct
def fun_set_aero_mesh(mesh:List,struct_var:List):
    nmesh = len(mesh)
    print(f"fun_set_aero_mesh: nmesh = {nmesh}")
    if nmesh > 0:
        imat = mesh[nmesh-1].imat
        mesh.extend([deepcopy(MeshStruct()) for i in range(len(struct_var))])
    else:
        mesh = [deepcopy(MeshStruct()) for i in range(len(struct_var))]
    for i in range(len(struct_var)):
        mesh[nmesh+i].nodes        = struct_var[i].grid.aero.X; 
        mesh[nmesh+i].connectivity = struct_var[i].grid.aero.connectivity; 
        mesh[nmesh+i].mx           = struct_var[i].grid.aero.N; 
        mesh[nmesh+i].my           = struct_var[i].grid.aero.M
        mesh[nmesh+i].nn           = (mesh[nmesh+i].mx+1)*(mesh[nmesh+i].my+1);
    return mesh