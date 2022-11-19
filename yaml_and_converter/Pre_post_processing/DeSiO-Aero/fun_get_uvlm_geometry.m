% =================================================================================================================
% Function to extract geometry from WindIO for aerodynamical grid in DeSiO-Format.
% 
% Author: Christian Hente
% Date: 05.05.2022
%
% input:
%   strName - type of surface of (airfoil) cross-section
%   model - component uvlm object specified in WindIO
%   airfoils - airfoil data specified in WindIO
%   scale_opt - scaling factor for scaling in longitudinal direction (optional input)
% output:
%     uvlm_ob - uvlm object containing coordinates and connectivity for creating aerodynamic grid
% =================================================================================================================
function uvlm_ob = fun_get_uvlm_geometry_test(strName,model,airfoils,scale_opt)
% =================================================================================================================
    switch nargin
        case 3
            scale = 1;
        case 4
            scale = scale_opt;
    end
    
    % extracting uvlm data from WindIO
    uvlm_ob = fun_extract_uvlm_data(strName,model,airfoils);
    
    % elements in x (span) and y (chord) direction
    mx =  length(uvlm_ob.arr_xhi_x)-1; 
    my =  uvlm_ob.N;
    
    % natural coordinates in chord-wise direction of cross-section
    xhi_y_c   = uvlm_ob.xhi_y_c;                                            % camber surface
    xhi_y_w   = [uvlm_ob.xhi_y_w(1:end),uvlm_ob.xhi_y_w(end-1:-1:1)];       % whole surface of cross-section
    
    % some extracted geometry arrays interpolated to spatial discretization
    arr_chord      = uvlm_ob.arr_c;                                         % chord length c
    arr_pitch_ax   = uvlm_ob.arr_pitch_ax;                                  % position of pitch axis in c
    arr_xhi_airf_c = uvlm_ob.arr_xhi_airf_c;                                % camber surface of (ariofoil) cross-section in c
    arr_xhi_airf_w = uvlm_ob.arr_xhi_airf_w;                                % whole surface of (ariofoil) cross-section in c
    arr_twist      = -uvlm_ob.arr_twist;                                    % twist angle around pitch axis
    arr_xre        = [uvlm_ob.arr_xre_x,uvlm_ob.arr_xre_y,...
                                            uvlm_ob.arr_xre_z*scale];       % reference axis
    
    % position and director coordinates of the reference axis
    for i = 1:mx+1
        % save the old directors from the previous segment
        if i == 1
            d1_old = [1; 0; 0];
            d2_old = [0; 1; 0];
            d3_old = [0; 0; 1];
            alpha_old = 0;
        else
            d1_old = d1';   % transposed because they were transposed for the output
            d2_old = d2';
            d3_old = d3';
            alpha_old = arr_twist(i-1);
        end
        
        % compute the new 3-director
        n3 = zeros(3, 1);   % the current connection vector
        if i == mx+1
            n3(:) = arr_xre(i,:) - arr_xre(i-1,:);  % exeption at last node: connection interpolated backwards
        else
            n3(:) = arr_xre(i+1,:) - arr_xre(i,:);  % usual case: connection interpolated forward
        end
        d3 = n3/norm(n3);   % 3-director as normed connection vector
        
        % compute the angle between old and new 3-director and the
        % corresponding rotation axis
        theta = acos(dot(d3_old, d3));   % rotation angle
%         if i == model.M_struc+1
%             disp('rotation angle:');
%             disp(theta);
%         end
        if theta <= 1e-6   % special case (for the last segment): the old and new director are very similar --> then just use one of the old directors as 
            rot_matrix = eye(3);
        else
            rot_axis = cross(d3_old, d3)/norm(cross(d3_old, d3));   % direction vector of the rotation axis
            u_x = rot_axis(1);  % 1-component of the rotation axis
            u_y = rot_axis(2);  % 2-component of the rotation axis
            u_z = rot_axis(3);  % 3-component of the rotation axis
            
%             disp(rot_axis);
%         if i == model.M_struc+1
%             disp('rotation axis:');
%             disp(rot_axis);
%         end
        % compute the rotation matrix
            rot_matrix = [
            cos(theta) + u_x^2*(1-cos(theta)),          u_x*u_y*(1-cos(theta)) - u_z*sin(theta),	u_x*u_z*(1-cos(theta)) + u_y*sin(theta);
            u_y*u_x*(1-cos(theta)) + u_z*sin(theta),	cos(theta) + u_y^2*(1-cos(theta)),          u_y*u_z*(1-cos(theta)) - u_x*sin(theta);
            u_z*u_x*(1-cos(theta)) - u_y*sin(theta),    u_z*u_y*(1-cos(theta)) + u_x*sin(theta),    cos(theta) + u_z^2*(1-cos(theta))];
        end
%         disp(rot_matrix);
    
%         if i == model.M_struc+1
%             disp('rotation matrix:');
%             disp(rot_matrix);
%         end
        
        % rotate the directors to get the new COS (without twist)
        d3_check = mtimes(rot_matrix, d3_old);
        d3_diff = d3 - d3_check;
        if norm(d3_diff)>1e-10
            disp('ERROR: rotation does not work as intended');
            disp(d3_old);
            disp(d3);
            disp(d3_check);
            return
        end
        d2_noTw = mtimes(rot_matrix, d2_old);
        d1_noTw = mtimes(rot_matrix, d1_old);
        
        % include the twist
        alpha_tw = arr_twist(i) - alpha_old;   % current twist angle
        
        % rotation axis for the twist rotation: the current 3-director
        rot_axis = d3;
        u_x = rot_axis(1);  % 1-component of the rotation axis
        u_y = rot_axis(2);  % 2-component of the rotation axis
        u_z = rot_axis(3);  % 3-component of the rotation axis
        
        % compute the rotation matrix for the twist rotation
        rot_matrix_Tw = [
        cos(alpha_tw) + u_x^2*(1-cos(alpha_tw)),          u_x*u_y*(1-cos(alpha_tw)) - u_z*sin(alpha_tw),	u_x*u_z*(1-cos(alpha_tw)) + u_y*sin(alpha_tw);
        u_y*u_x*(1-cos(alpha_tw)) + u_z*sin(alpha_tw),	cos(alpha_tw) + u_y^2*(1-cos(alpha_tw)),          u_y*u_z*(1-cos(alpha_tw)) - u_x*sin(alpha_tw);
        u_z*u_x*(1-cos(alpha_tw)) - u_y*sin(alpha_tw),    u_z*u_y*(1-cos(alpha_tw)) + u_x*sin(alpha_tw),    cos(alpha_tw) + u_z^2*(1-cos(alpha_tw))];
        
        % compute the new directors
        d1 = mtimes(rot_matrix_Tw, d1_noTw);
        d2 = mtimes(rot_matrix_Tw, d2_noTw);
        
        % transpose the vectors for output
        d1 = d1';
        d2 = d2';
        d3 = d3';
        
        % write the computed current director and reference point to the output-array
        arr_coordinates(i,1:12) = [arr_xre(i,1:3),d1,d2,d3];
    end
    
    % calculating aerodynamic grid for DeSiO-Aero on the reference axis
    d10 = [1;0;0]; d20 = [0;1;0]; d30 = [0;0;1];
    for i_l = 1:size(arr_coordinates,1)
        
        % reference position and director of/at center line
        xcenter = arr_coordinates(i_l,1:3)'; 
        d1 = arr_coordinates(i_l,4:6)'; 
        d2 = arr_coordinates(i_l,7:9)'; 
        d3 = arr_coordinates(i_l,10:12)';

        % rotation matrix for cross section orientation
        R_t = d3*d30' + d2*d20' + d1*d10';
        
        % chord length and location of pitch axis
        c        = arr_chord(i_l); 
        xhi_p_ax = arr_pitch_ax(i_l);
        xhi_pn   = [0;xhi_p_ax;0];
        
        % coordinates of whole cross-section
        xcO = []; xhi_pnO = [];
        O              = ones(3,length(xhi_y_w)); 
        xcO(1:3,:)     = [xcenter(1)*O(1,:);xcenter(2)*O(2,:);xcenter(3)*O(3,:)];
        xhi_pnO(1:3,:) = [xhi_pn(1)*O(1,:);xhi_pn(2)*O(2,:);xhi_pn(3)*O(3,:)];
        xhi_airf_w     = [arr_xhi_airf_w(i_l,:); xhi_y_w; 0*O(3,:)];
        x_airf_w       = [xhi_airf_w(1,:)*c; xhi_airf_w(2,:)*c; 0*O(3,:)];
        % rotation of cross-section around pitch axis and recalcuating profile 
        % to middle axis (middle axis is on location of pitch axis)
        xp_w   =  R_t*(x_airf_w - xhi_pnO*c) + xcO ;
        
        % coordinates of camber line
        xcO = []; xhi_pnO = [];
        O              = ones(3,length(xhi_y_c)); 
        xcO(1:3,:)     = [xcenter(1)*O(1,:);xcenter(2)*O(2,:);xcenter(3)*O(3,:)];
        xhi_pnO(1:3,:) = [xhi_pn(1)*O(1,:);xhi_pn(2)*O(2,:);xhi_pn(3)*O(3,:)];
        xhi_airf_c     = [arr_xhi_airf_c(i_l,:); xhi_y_c; 0*O(3,:)];
        x_airf_c       = [xhi_airf_c(1,:)*c; xhi_airf_c(2,:)*c; 0*O(3,:)];
        % rotation of cross-section around pitch axis and recalcuating profile 
        % to middle axis (middle axis is on location of pitch axis)
        xp_c =  R_t*(x_airf_c - xhi_pnO*c) + xcO ;
        xp_0 = xcO;
                
        % global node indizes for assembling and meshing
        inz_global_w = [(i_l-1)*(2*my+1)+1:i_l*(2*my+1)];
        inz_global_c = [(i_l-1)*(my+1)+1:i_l*(my+1)];

        % writing into global coordinates
        uvlm_ob.X_W(inz_global_w,1:3) = xp_w';
        uvlm_ob.X_C(inz_global_c,1:3) = xp_c';
        % this is only for assigning chord-length each node for fsi radius, see line 236
        uvlm_ob.X_0(inz_global_c,1:3)   = xp_0';
        uvlm_ob.X_0_W(inz_global_w,1:3) = [xp_0(1,1)*ones(2*my+1,1),xp_0(2,1)*ones(2*my+1,1),xp_0(3,1)*ones(2*my+1,1)];
        
%         plot3(xp_w(1,:),xp_w(2,:),xp_w(3,:))
%         text(xp_w(1,:),xp_w(2,:),xp_w(3,:),num2str([1:size(xp_w,2)]'))
    end

    % creating camber surface connectivities
    %       2 ___________ 1
    %        |           |
    %        |    BE     |   
    %        |___________|
    %        3           4
    surfaces = [1:mx*(my)]; connectivity_c = [];
    for j = 1:mx
        n2 = surfaces((j-1)*my+1:(j-1)*my+my) + my + j;
        n1 = n2+1;
        n3 = n2-(my+1);
        n4 = n2-my;
        connectivity_c = [connectivity_c; [n1',n2',n3',n4']];
    end
    
    figure; hold on; grid on;
    for i_air = 1:size(connectivity_c,1)
        surf = fill3(uvlm_ob.X_C(connectivity_c(i_air,:),1), uvlm_ob.X_C(connectivity_c(i_air,:),2), uvlm_ob.X_C(connectivity_c(i_air,:),3),[0.00,0.60,0.60],'facealpha',1.0,'edgealpha',0.1);
    end
    text(uvlm_ob.X_C(:,1),uvlm_ob.X_C(:,2),uvlm_ob.X_C(:,3),num2str([1:size(uvlm_ob.X_C,1)]'));
    
    % creating surface connectivities of whole cross-section
    surfaces = [1:mx*(2*my)]; connectivity_w = [];
    for j = 1:mx
        n2 = surfaces((j-1)*2*my+1:(j-1)*2*my+2*my) + 2*my + j;
        n1 = n2+1;
        n3 = n2-(2*my+1);
        n4 = n2-2*my;
        connectivity_w = [connectivity_w; [n1',n2',n3',n4']];
    end

%     figure; hold on; grid on;
%     for i_air = 1:size(connectivity_w,1)
%         surf = fill3(uvlm_ob.X_W(connectivity_w(i_air,:),1), uvlm_ob.X_W(connectivity_w(i_air,:),2), uvlm_ob.X_W(connectivity_w(i_air,:),3),[0.60,0.60,0.60],'facealpha',0.1,'edgealpha',0.1);
%     end
%     for i = 1:size(uvlm_ob.X_W,1)
%         text(uvlm_ob.X_W(i,1),uvlm_ob.X_W(i,2),uvlm_ob.X_W(i,3),num2str(i));
%     end

    % assigning ring radius for whole cross-section, depending on chord
    uvlm_ob.arr_node_fsi_radius_w = [];
    for i = 1:(mx+1)*(2*my+1)
        [inz] = find(abs(uvlm_ob.X_0_W(i,3)) >= uvlm_ob.arr_xre_z);
        uvlm_ob.arr_node_fsi_radius_w(i) = uvlm_ob.arr_c(inz(end));
    end
    
    % assigning ring radius for (chamber) surface, depending on chord
    uvlm_ob.arr_node_fsi_radius = [];
    for i = 1:(mx+1)*(my+1)
        [inz] = find(abs(uvlm_ob.X_0(i,3)) >= uvlm_ob.arr_xre_z);
        uvlm_ob.arr_node_fsi_radius(i) = uvlm_ob.arr_c(inz(end));
    end
%     figure; plot(uvlm_ob.arr_node_fsi_radius);
    
    uvlm_ob.connectivity_c = connectivity_c;
    uvlm_ob.connectivity_w = connectivity_w;
    uvlm_ob.N = my;
    uvlm_ob.M = mx;
% =================================================================================================================
return