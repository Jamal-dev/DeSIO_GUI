import numpy as np
from functions import fun_extract_uvlm_data
from math import cos, sin
import logging
def fun_get_uvlm_geometry(strName,model,airfoils,scale_opt=1):
    """
         Input:
           airfoils - airfoil data specified in WindIO
           scale_opt - scaling factor for scaling in longitudinal direction (optional input)
         output:
             uvlm_ob - uvlm object containing coordinates and connectivity for creating aerodynamic grid
         =================================================================================================================

    """
    # extracting uvlm data from WindIO
    uvlm_ob = fun_extract_uvlm_data(strName,model,airfoils)
    
    # elements in x (span) and y (chord) direction
    mx =  len(uvlm_ob.arr_xhi_x)-1; 
    my =  uvlm_ob.N;
    
    # natural coordinates in chord-wise direction of cross-section
    xhi_y_c   = uvlm_ob.xhi_y_c;                                            # camber surface
    xhi_y_w   = np.append(uvlm_ob.xhi_y_w, uvlm_ob.xhi_y_w[-2::-1])       # whole surface of cross-section
    
    # some extracted geometry arrays interpolated to spatial discretization
    arr_chord      = uvlm_ob.arr_c;                                         # chord length c
    arr_pitch_ax   = uvlm_ob.arr_pitch_ax;                                  # position of pitch axis in c
    arr_xhi_airf_c = uvlm_ob.arr_xhi_airf_c;                                # camber surface of (ariofoil) cross-section in c
    arr_xhi_airf_w = uvlm_ob.arr_xhi_airf_w;                                # whole surface of (ariofoil) cross-section in c
    arr_twist      = -uvlm_ob.arr_twist;                                    # twist angle around pitch axis
    arr_xre        = np.vstack([uvlm_ob.arr_xre_x,uvlm_ob.arr_xre_y, uvlm_ob.arr_xre_z*scale_opt]).T       # reference axis
    arr_coordinates = np.zeros((mx+1,12))
    # position and director coordinates of the reference axis
    for i in range(mx+1):
        # save the old directors from the previous segment
        if i == 0:
            d1_old = np.asarray([[1], [0], [0]])
            d2_old = np.asarray([[0], [1], [0]])
            d3_old = np.asarray([[0], [0], [1]])
            alpha_old = 0;
        else:
            d1_old = d1   # transposed because they were transposed for the output
            d2_old = d2
            d3_old = d3
            alpha_old = arr_twist[i-1]


        # compute the new 3-director
        n3 = np.zeros((3, 1))   # the current connection vector
        if i == mx:
            n3[:,0] = arr_xre[i,:] - arr_xre[i-1,:];  # exeption at last node: connection interpolated backwards
        else:
            n3[:,0] = arr_xre[i+1,:] - arr_xre[i,:];  # usual case: connection interpolated forward

        d3 = n3/np.linalg.norm(n3);   # 3-director as normed connection vector

        # compute the angle between old and new 3-director and the
        # corresponding rotation axis
        theta = np.arccos(np.dot(d3_old[:,0], d3[:,0]));   # rotation angle
        #         if i == model.M_struc+1
        #             print('rotation angle:');
        #             print(theta);
        #         
        
        if theta <= 1e-6 or np.isnan(theta):   # special case (for the last segment): the old and new director are very similar --> then just use one of the old directors as 
            rot_matrix = np.eye(3)
        else:
            temp = np.cross(d3_old[:,0], d3[:,0])
            rot_axis = temp/np.linalg.norm(temp);   # direction vector of the rotation axis
            u_x = float(rot_axis[0])  # 1-component of the rotation axis
            u_y = float(rot_axis[1])  # 2-component of the rotation axis
            u_z = float(rot_axis[2])  # 3-component of the rotation axis
            
        #             print(rot_axis);
        #         if i == model.M_struc+1
        #             print('rotation axis:');
        #             print(rot_axis);
        #         
        # compute the rotation matrix
            rot_matrix = [
            [cos(theta) + u_x**2*(1-cos(theta)),          u_x*u_y*(1-cos(theta)) - u_z*sin(theta),	u_x*u_z*(1-cos(theta)) + u_y*sin(theta)],
            [u_y*u_x*(1-cos(theta)) + u_z*sin(theta),	cos(theta) + u_y**2*(1-cos(theta)),          u_y*u_z*(1-cos(theta)) - u_x*sin(theta)],
            [u_z*u_x*(1-cos(theta)) - u_y*sin(theta),    u_z*u_y*(1-cos(theta)) + u_x*sin(theta),    cos(theta) + u_z**2*(1-cos(theta))]]
            rot_matrix = np.asarray(rot_matrix)

        #         print(rot_matrix);

        #         if i == model.M_struc+1
        #             print('rotation matrix:');
        #             print(rot_matrix);
        #         

        # rotate the directors to get the new COS (without twist)
        d3_check = np.matmul(rot_matrix, d3_old);
        d3_diff = d3 - d3_check;
        if np.linalg.norm(d3_diff)>1e-10:
            print('ERROR: rotation does not work as inteded');
            print(d3_old);
            print(d3);
            print(d3_check);
            return

        d2_noTw = rot_matrix @ d2_old
        d1_noTw = rot_matrix @ d1_old

        # include the twist
        alpha_tw = arr_twist[i] - alpha_old;   # current twist angle

        # rotation axis for the twist rotation: the current 3-director
        rot_axis = d3
        u_x = float(rot_axis[0])  # 1-component of the rotation axis
        u_y = float(rot_axis[1])  # 2-component of the rotation axis
        u_z = float(rot_axis[2])  # 3-component of the rotation axis

        # compute the rotation matrix for the twist rotation
        rot_matrix_Tw = [
        [cos(alpha_tw) + u_x**2*(1-cos(alpha_tw)),          u_x*u_y*(1-cos(alpha_tw)) - u_z*sin(alpha_tw),	u_x*u_z*(1-cos(alpha_tw)) + u_y*sin(alpha_tw)],
        [u_y*u_x*(1-cos(alpha_tw)) + u_z*sin(alpha_tw),	cos(alpha_tw) + u_y**2*(1-cos(alpha_tw)),          u_y*u_z*(1-cos(alpha_tw)) - u_x*sin(alpha_tw)],
        [u_z*u_x*(1-cos(alpha_tw)) - u_y*sin(alpha_tw),    u_z*u_y*(1-cos(alpha_tw)) + u_x*sin(alpha_tw),    cos(alpha_tw) + u_z**2*(1-cos(alpha_tw))]]
        rot_matrix_Tw = np.asarray(rot_matrix_Tw)

        # print(f'i={i}, d1_noTw = {d1_noTw}')
        # print(f'i={i}, d1_noTw = {d2_noTw}')
        # compute the new directors
        d1 = rot_matrix_Tw @ d1_noTw
        d2 = rot_matrix_Tw @ d2_noTw

        # transpose the vectors for output
        

        # write the computed current director and reference point to the output-array
        arr_coordinates[i,:] = np.hstack([arr_xre[i,:],d1[:,0],d2[:,0],d3[:,0]])
    
    
    # calculating aerodynamic grid for DeSiO-Aero on the reference axis
    d10 = np.asarray([[1],[0],[0]] )
    d20 = np.asarray([[0],[1],[0]])
    d30 = np.asarray([[0],[0],[1]])
    for i_l in range(arr_coordinates.shape[0]):
        
        # reference position and director of/at center line
        xcenter = arr_coordinates[i_l:i_l+1,:3].T 
        d1 = arr_coordinates[i_l:i_l+1,3:6].T 
        d2 = arr_coordinates[i_l:i_l+1,6:9].T 
        d3 = arr_coordinates[i_l:i_l+1,9:12].T

        # rotation matrix for cross section orientation
        R_t = d3@d30.T + d2@d20.T + d1@d10.T
        
        # chord length and location of pitch axis
        c        = arr_chord[i_l] 
        xhi_p_ax = arr_pitch_ax[i_l]
        xhi_pn   = np.asarray([[0],[xhi_p_ax],[0]])

        # print(f'xcenter={xcenter, xcenter.shape}, xhi_p_ax={xhi_p_ax}, xhi_pn= {xhi_pn, xhi_pn.shape}, c={c}')
        # print(f'R_t={R_t}')
        # print(f'd1={d1,d1.shape}')
        # print(f'd2={d2, d2.shape}')
        # print(f'd3={d3, d3.shape}')
        
        # coordinates of whole cross-section
        dim_v = len(xhi_y_w)
        xcO = np.zeros((3,dim_v)); 
        xhi_pnO = np.zeros((3,dim_v))
        O              = np.ones((3,dim_v)) 
        xcO[:3,:]      = np.asarray([xcenter[0]*O[0,:],      xcenter[1]*O[1,:],  xcenter[2]*O[2,:]])
        xhi_pnO[:3,:]  = np.asarray([xhi_pn[0]*O[0,:],       xhi_pn[1]*O[1,:],   xhi_pn[2]*O[2,:]])
        xhi_airf_w     = np.asarray([arr_xhi_airf_w[i_l,:],  xhi_y_w,            0*O[2,:]])
        x_airf_w       = np.asarray([xhi_airf_w[0,:]*c,      xhi_airf_w[1,:]*c,  0*O[2,:]])
        
        # print(f'xcO={xcO, xcO.shape}')
        # print(f'xhi_pnO={xhi_pnO, xhi_pnO.shape}')
        # print(f'x_airf_w={x_airf_w, x_airf_w.shape}')
        # print(f'xhi_airf_w={xhi_airf_w, xhi_airf_w.shape}')


        # rotation of cross-section around pitch axis and recalcuating profile 
        # to middle axis (middle axis is on location of pitch axis)
        xp_w =  (R_t @ (x_airf_w - xhi_pnO*c)) + xcO 
        
        # print(f'xp_w={xp_w, xp_w.shape}')

        # coordinates of camber line
        dim_v = len(xhi_y_c)
        xcO = np.zeros((3,dim_v)); 
        xhi_pnO = np.zeros((3,dim_v))
        O              = np.ones((3,dim_v)) 
        xcO[:3,:]      = np.asarray([xcenter[0]*O[0,:],          xcenter[1]*O[1,:],         xcenter[2]*O[2,:]])
        xhi_pnO[:3,:]  = np.asarray([xhi_pn[0]*O[0,:],           xhi_pn[1]*O[1,:],          xhi_pn[2]*O[2,:]])
        xhi_airf_c     = np.asarray([arr_xhi_airf_c[i_l,:],      xhi_y_c,                   0*O[2,:]])
        x_airf_c       = np.asarray([xhi_airf_c[0,:]*c,          xhi_airf_c[1,:]*c,         0*O[2,:]])
        # rotation of cross-section around pitch axis and recalcuating profile 
        # to middle axis (middle axis is on location of pitch axis)
        xp_c =  (R_t @ (x_airf_c - xhi_pnO*c) )+ xcO 
        xp_0 = xcO
                
        # global node indizes for assembling and meshing
        inz_global_w = np.arange((i_l)*(2*my+1),  (i_l+1)*(2*my+1))
        inz_global_c = np.arange((i_l)*(my+1)  ,  (i_l+1)*(my+1))

        # writing into global coordinates
        # print(f'i = {i_l}, xp_w = {xp_w}')
        uvlm_ob.X_W[inz_global_w,:3] = xp_w.T
        uvlm_ob.X_C[inz_global_c,:3] = xp_c.T
        # this is only for assigning chord-length each node for fsi radius
        uvlm_ob.X_0[inz_global_c,:3]   = xp_0.T
        # logging.debug(f'uvlm_ob.X_0.shape = {uvlm_ob.X_0.shape}')
        temp = np.ones((2*my+1,1))
        # logging.debug(f'xp_0.shape = {xp_0.shape}')
        uvlm_ob.X_0_W[inz_global_w,:3] = np.hstack([xp_0[0,0] * temp,
                                            xp_0[1,0] * temp,
                                            xp_0[2,0] * temp])
        # logging.debug(f'uvlm_ob.X_0_W.shape = {uvlm_ob.X_0_W.shape}')
        
    
    logging.debug(f'uvlm_ob.X_0.shape = {uvlm_ob.X_0.shape}')
    logging.debug(f'uvlm_ob.X_0_W.shape = {uvlm_ob.X_0_W.shape}')
    logging.debug(f'xp_0.shape = {xp_0.shape}')
    # creating camber surface connectivities
    #       2 ___________ 1
    #        |           |
    #        |    BE     |   
    #        |___________|
    #        3           4
    surfaces = np.arange(0,mx*(my), dtype=int) 
    # print('surfaces.shape = ',surfaces.shape)
    connectivity_c = np.zeros((surfaces.size,4), dtype = int)    
    for j  in range(mx):
        n2 = surfaces[ j*my : j*my+ my ] + my + j+1
        n1 = n2+1
        n3 = n2-(my+1)
        n4 = n2-my
        temp = np.vstack((n1,n2,n3,n4)).T
        connectivity_c[j*my:j*my+my,:] = temp
    # creating surface connectivities of whole cross-section
    surfaces = np.arange(0,mx*(2*my), dtype=int)
    # print('surfaces.shape = ',surfaces.shape)
    connectivity_w = np.zeros((surfaces.size,4), dtype = int)

    for j  in range(mx):
        n2 = surfaces[ j * 2 * my : j * 2 * my + 2 * my ] + 2 * my + j + 1
        n1 = n2+1
        n3 = n2-(2*my+1)
        n4 = n2-2*my
        temp = np.vstack((n1,n2,n3,n4)).T
        connectivity_w[j*2*my:j*2*my+2*my,:] = temp
    # print('connectivity_w.shape=',connectivity_w.shape)
    uvlm_ob.connectivity_c = connectivity_c;
    uvlm_ob.connectivity_w = connectivity_w;
    uvlm_ob.N = my;
    uvlm_ob.M = mx;
    uvlm_ob.arr_node_fsi_radius_w = []
    for i in range((mx+1)*(2*my+1)):
        inz = np.nonzero(abs(uvlm_ob.X_0_W[i,2]) >= uvlm_ob.arr_xre_z)[0]
    
        uvlm_ob.arr_node_fsi_radius_w.append(uvlm_ob.arr_c[inz[-1]])
    uvlm_ob.arr_node_fsi_radius_w = np.asarray(uvlm_ob.arr_node_fsi_radius_w)
    
    logging.debug(f'uvlm_ob.arr_node_fsi_radius_w.shape = {uvlm_ob.arr_node_fsi_radius_w.shape}')
    logging.debug(f'uvlm_ob.arr_node_fsi_radius_w = {uvlm_ob.arr_node_fsi_radius_w}')
    uvlm_ob.arr_node_fsi_radius = [];
    for i in range((mx+1)*(my+1)):
        inz = np.nonzero(abs(uvlm_ob.X_0[i,2]) >= uvlm_ob.arr_xre_z)[0]
        uvlm_ob.arr_node_fsi_radius.append(uvlm_ob.arr_c[inz[-1]])
    uvlm_ob.arr_node_fsi_radius = np.asarray(uvlm_ob.arr_node_fsi_radius)
    logging.debug(f'uvlm_ob.arr_node_fsi_radius.shape = {uvlm_ob.arr_node_fsi_radius.shape}')
    logging.debug(f'uvlm_ob.arr_node_fsi_radius = {uvlm_ob.arr_node_fsi_radius}')

    return uvlm_ob