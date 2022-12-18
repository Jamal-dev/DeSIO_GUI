import numpy as np
from copy import deepcopy
from classes.data.data_structs import MeshStruct
def  fun_set_struc_mesh(mesh,struct_var):


    nmesh = len(mesh)
    print(f"fun_set_struc_mesh: nmesh = {nmesh}")
    imat  = 0
    if nmesh > 0:
        imat = mesh[nmesh-1].imat
        mesh.append(deepcopy(MeshStruct()))
    else:
        mesh = [deepcopy(MeshStruct())]
    
    i = 0
    mesh[nmesh].nodes               = np.hstack([struct_var[i].grid.struc.X_RE,                                             struct_var[i].grid.struc.D1,                                             struct_var[i].grid.struc.D2,                                             struct_var[i].grid.struc.D3])
    assert mesh[nmesh].nodes.shape[1] == 12, "fun_set_struc_mesh, nodes dimensions are not correct"
    mesh[nmesh].mx                  = struct_var[i].grid.struc.M; 
    mesh[nmesh].nn                  = (struct_var[i].grid.struc.M+1);
    mesh[nmesh].connectivity = np.zeros(shape = struct_var[i].grid.struc.connectivity.shape)
    mesh[nmesh].connectivity[:,:2] = struct_var[i].grid.struc.connectivity[:,:2]
    # element connectivity and material
    mesh[nmesh].connectivity[:,2]   = np.asarray(struct_var[i].grid.struc.matBeam.inz) + imat;
    mesh[nmesh].imat                = max(mesh[nmesh].connectivity[:,2]);
    return mesh