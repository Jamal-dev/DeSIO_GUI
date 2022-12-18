from math import cos, sin
import numpy as np
from functions import fun_extract_beam_data
def fun_get_beam_model(strName,model,materials,scale_opt=1):
    """
         input:
           strName - type of surface of (airfoil) cross-section
           model - component beam object specified in WindIO
           materials - material data specified in WindIO
           scale_opt - scaling factor for scaling in longitudinal direction (optional input)
         output:
           beam - beam object containing coordinates and connectivity for creating structural mesh
    """
    beam  = fun_extract_beam_data(strName,model,materials)
    # print(beam)
    scale = scale_opt
    # twist angle around pitch axis
    arr_twist = np.zeros((beam.M+1,1))
    if 'twist' in model.keys():
        arr_twist = -beam.arr_twist;


    # reference axis
    arr_xre   = np.vstack([beam.arr_xre_x,                             beam.arr_xre_y,                             beam.arr_xre_z*scale]).T

    # beam position and director coordinates of the reference axis
    beam.arr_coordinates = [];
    for i in range(model["M_struc"]+1):

        # save the old directors from the previous segment
        if i == 0:
            d1_old = np.asarray([[1], [0], [0]]);
            d2_old = np.asarray([[0], [1], [0]]);
            d3_old = np.asarray([[0], [0], [1]]);
            alpha_old = 0;
        else:

            alpha_old = arr_twist[i-1];
            d1_old = d1;
            d2_old = d2;
            d3_old = d3;

        assert d1_old.ndim == 2, "d1 should have 2 dim array"
        assert d2_old.ndim == 2, "d2 should have 2 dim array"
        assert d3_old.ndim == 2, "d3 should have 2 dim array"
        # compute the new 3-director
        n3 = np.zeros((3, 1));   # the current connection vector
        if i == model["M_struc"]:
            n3[:,0] = arr_xre[i,:] - arr_xre[i-1,:];  # exeption at last node: connection interpolated backwards
        else:
            n3[:,0] = arr_xre[i+1,:] - arr_xre[i,:];  # usual case: connection interpolated forward

        d3 = n3/np.linalg.norm(n3);   # 3-director as normed connection vector

        # compute the angle between old and new 3-director and the
        # corresponding rotation axis
        theta = np.arccos(np.dot(d3_old.flatten(), d3.flatten()));   # rotation angle
        # if i == model["M_struc"]:
        #     print(f"theta = {theta}")
        # compute the rotation axis for the computation with the new
        # 3-director if the rotation angle is not very small
        if theta <= 1e-6 or np.isnan(theta):   # special case (for the last segment): the old and new director are very similar --> then just use one of the old directors as 
            rot_matrix = np.eye(3);
        else:
            temp = np.cross(d3_old.flatten(), d3.flatten())
            rot_axis = temp/np.linalg.norm(temp);   # direction vector of the rotation axis
            u_x = rot_axis[0];  # 1-component of the rotation axis
            u_y = rot_axis[1];  # 2-component of the rotation axis
            u_z = rot_axis[2];  # 3-component of the rotation axis

        # compute the rotation matrix
            rot_matrix = np.asarray([
            [cos(theta) + u_x**2*(1-cos(theta)),          u_x*u_y*(1-cos(theta)) - u_z*sin(theta),	u_x*u_z*(1-cos(theta)) + u_y*sin(theta)],
            [u_y*u_x*(1-cos(theta)) + u_z * sin(theta),    cos(theta) + u_y**2*(1-cos(theta)),          u_y*u_z*(1-cos(theta)) - u_x*sin(theta)],
            [u_z*u_x*(1-cos(theta)) - u_y*sin(theta),    u_z*u_y*(1-cos(theta)) + u_x*sin(theta),    cos(theta) + u_z**2*(1-cos(theta))]
                    ]);


        # rotate the directors to get the new COS (without twist)
        d3_check = rot_matrix @ d3_old;
        d3_diff = d3 - d3_check;
        if np.linalg.norm(d3_diff)>1e-10:
            print('ERROR: rotation does not work as inted');
            print(d3_old);
            print(d3);
            print(d3_check);
            return

        d2_noTw = rot_matrix @ d2_old
        d1_noTw = rot_matrix @ d1_old

        # include the twist
        alpha_tw = arr_twist[i] - alpha_old;   # current twist angle

        # rotation axis for the twist rotation: the current 3-director
        rot_axis = d3;
        u_x = float(rot_axis[0])  # 1-component of the rotation axis
        u_y = float(rot_axis[1])  # 2-component of the rotation axis
        u_z = float(rot_axis[2])  # 3-component of the rotation axis

        # compute the rotation matrix for the twist rotation
        rot_matrix_Tw = np.asarray([
        [cos(alpha_tw) + u_x**2*(1-cos(alpha_tw)),          u_x*u_y*(1-cos(alpha_tw)) - u_z*sin(alpha_tw),	u_x*u_z*(1-cos(alpha_tw)) + u_y*sin(alpha_tw)],
        [u_y*u_x*(1-cos(alpha_tw)) + u_z*sin(alpha_tw),	cos(alpha_tw) + u_y**2*(1-cos(alpha_tw)),          u_y*u_z*(1-cos(alpha_tw)) - u_x*sin(alpha_tw)],
        [u_z*u_x*(1-cos(alpha_tw)) - u_y*sin(alpha_tw),    u_z*u_y*(1-cos(alpha_tw)) + u_x*sin(alpha_tw),    cos(alpha_tw) + u_z**2*(1-cos(alpha_tw))]
        ]);

        # compute the new directors
        d1 = rot_matrix_Tw @ d1_noTw
        d2 = rot_matrix_Tw @ d2_noTw

        # transpose the vectors for output
#         d1 = d1';
#         d2 = d2';
#         d3 = d3';

        # write the computed current director and reference point to the output-array
        temp = np.expand_dims(arr_xre[i,:],0)
        temp = np.hstack([temp,d1.T,d2.T,d3.T])
        # print(f"temp.shape = {temp.shape}")
        beam.arr_coordinates.append(temp[0,:])
        



    beam.arr_coordinates = np.asarray(beam.arr_coordinates)
    # beam connectivities
    for i in range(model["M_struc"]):
        beam.connectivity[i,:] = [i,i+1,i];
    return beam
