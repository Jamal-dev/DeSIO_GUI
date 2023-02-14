% =================================================================================================================
% Function to extract uvlm data from WindIO for aerodynamical grid in DeSiO-Format.
% 
% Author: Christian Hente
% Date: 05.05.2022
%
% input:
%   strName - type of surface of (airfoil) cross-section
%   model - component uvlm object specified in WindIO
%   airfoils - airfoil data specified in WindIO
% output:
%     uvlm_ob - uvlm object containing coordinates and connectivity for creating aerodynamic grid
% =================================================================================================================
function uvlm_obj = fun_extract_uvlm_data(strName,model_uvlm,airfoils)
% =================================================================================================================
    % span-wise and chord-wise discretization
    uvlm_obj.M =  model_uvlm.M_aero; % span-wise
    uvlm_obj.N =  model_uvlm.N_aero; % chord-wise
    
    % natural coordinates of span-and chord-wise discretization
    xhi_x           = [0:1/uvlm_obj.M:1];   % span-wise
    xhi_y_c         = [0:1/(uvlm_obj.N):1]; % chord-wise
    
    nbr             = 0;                    % number of blade root cross-section
    arr_inz_airfoil = [];                   % index array for locating airfoils in span-wise direction
    
    % interpolating coordinates of reference axis according to span-wise discretization
    uvlm_obj.arr_xre_x = fun_lin_interpolation([model_uvlm.reference_axis.x.grid',model_uvlm.reference_axis.x.values'],xhi_x);
    uvlm_obj.arr_xre_y = fun_lin_interpolation([model_uvlm.reference_axis.y.grid',model_uvlm.reference_axis.y.values'],xhi_x);
    uvlm_obj.arr_xre_z = fun_lin_interpolation([model_uvlm.reference_axis.z.grid',model_uvlm.reference_axis.z.values'],xhi_x);
    
    % if-statement according to type of surface cross-section
    if strncmp(strName,'blade',length(strName))
        uvlm_obj.airfoil = model_uvlm.airfoil_position;
        nbr              = model_uvlm.blade_root_position;
        
        % interpolating chord length in span-wise direction
        uvlm_obj.arr_c     = fun_lin_interpolation([model_uvlm.chord.grid',model_uvlm.chord.values'],xhi_x);
        
        % interpolating twist angle according to span-wise discretization
        uvlm_obj.arr_twist = ones(uvlm_obj.M+1,1)*0.0;
        if isfield(model_uvlm,'twist')
            uvlm_obj.arr_twist = fun_lin_interpolation([model_uvlm.twist.grid',model_uvlm.twist.values'],xhi_x);
        end
        % interpolating airfoil in span-wise direction. This is important
        % to identify the blade root and blade's lifting surfaces
        airfoil_grid    = uvlm_obj.airfoil.grid;
        airfoil_values  = [1:length(uvlm_obj.airfoil.grid)];
        arr_inz_airfoil = fix((fun_lin_interpolation([airfoil_grid',airfoil_values'],xhi_x)));
        
        arr_pitch_ax = ones(uvlm_obj.M+1,1)*1.0;
        if isfield(model_uvlm,'pitch_axis')
            uvlm_obj.arr_pitch_ax = fun_lin_interpolation([model_uvlm.pitch_axis.grid',model_uvlm.pitch_axis.values'],xhi_x);
        end
        % natural coordinates of chord-wise discretization for whole
        % surface of cross-section
        xhi_y_w = xhi_y_c;
        
    elseif strncmp(strName,'pipe',length(strName))
        % interpolating chord length in span-wise direction
        uvlm_obj.arr_c          = fun_lin_interpolation([model_uvlm.outer_diameter.grid',model_uvlm.outer_diameter.values'-model_uvlm.thickness.values'],xhi_x);
        uvlm_obj.arr_twist      = ones(uvlm_obj.M+1,1)*0.0; % zero twist
        uvlm_obj.arr_pitch_ax   = ones(uvlm_obj.M+1,1)*0.5; % location of pitch axis
        uvlm_obj.airfoil.labels = {'circular';'circular'};  % artificial "airfoil" sections for pipe
        uvlm_obj.airfoil.grid   = [0.0, 1.0];               % artificial "airfoil" sections grid for pipe
        
        % mesh discretization around circular cross-section
        xhi_y_w = 0.5*(1-cos([0:(pi)/(uvlm_obj.N):pi]));
    end
    
    % loop over airfoils in the model to calculate coordinates of
    % cross-section surfaces
    for i = 1:size(uvlm_obj.airfoil.labels,1)
        airfoil_name = uvlm_obj.airfoil.labels(i);

        % searching current airfoil from airfoil-list
        for j = 1:size(airfoils,1)
            if size(airfoils,1) == 1
                if strcmp(airfoils.name,airfoil_name)
                    airfoilj = airfoils;
                    break
                end
            else
                if strcmp(airfoils{j}.name,airfoil_name)
                    airfoilj = airfoils{j};
                    break
                end
            end
        end

        % extracting coordinates for upper and lower airfoil 
        airfoil_coord = [airfoilj.coordinates.x',airfoilj.coordinates.y'];
        [inz0] = find(airfoil_coord(:,1)==0);
        airfoil_coord_u = airfoil_coord(1:inz0(1),:);   [val,inzsort] = sort(airfoil_coord_u(:,1)); airfoil_coord_u = airfoil_coord_u(inzsort,:);
        airfoil_coord_l = airfoil_coord(inz0(1):end,:); [val,inzsort] = sort(airfoil_coord_l(:,1)); airfoil_coord_l = airfoil_coord_l(inzsort,:);
        
        % interpolating values in airfoil thickness direction according to
        % discretization in chord-wise direction to calculate camber
        % surface coordinates
        xhi_airf_u = fun_lin_interpolation(airfoil_coord_u,xhi_y_c);
        xhi_airf_l = fun_lin_interpolation(airfoil_coord_l,xhi_y_c);
        % coordinates of camber surface
        xhi_airf_c = (xhi_airf_u(:,1) + xhi_airf_l(:,1))/2;
        
        % interpolating values in airfoil thickness direction according to
        % discretization in chord-wise direction to calculate whole
        % cross-section surface coordinates
        xhi_airf_u = fun_lin_interpolation(airfoil_coord_u,xhi_y_w);
        xhi_airf_l = fun_lin_interpolation(airfoil_coord_l,xhi_y_w);
        % coordinates of whole cross-section surface
        xhi_airf_w = [xhi_airf_u(1:end-1);(xhi_airf_u(end)+xhi_airf_l(end))/2;xhi_airf_l(end-1:-1:1)];

        % if-statement to detect, if blade root ends. In case blade root
        % end, switching from whole surface to camber surface
        if i == nbr+1
            if nbr ~= 0
                xhi_airf_w = [xhi_airf_c(1:end-1);xhi_airf_c(end:-1:1)];
            end
        end

        arr_val_c(:,i)   = xhi_airf_c;      % camber surface coordinates 
        arr_val_w(:,i)   = xhi_airf_w(:,1); % whole surface coordinates
    end

    % interpolating camber surface coordinates in span-wise direction 
    for i = 1:size(arr_val_c,1);
        arr_val = fun_lin_interpolation([uvlm_obj.airfoil.grid',arr_val_c(i,:)'],xhi_x);
        arr_xhi_airf_c(:,i) = arr_val(:,1);
    end
    
    % interpolating whole surface coordinates in span-wise direction
    for i = 1:size(arr_val_w,1);
        arr_val = fun_lin_interpolation([uvlm_obj.airfoil.grid',arr_val_w(i,:)'],xhi_x);
        arr_xhi_airf_w(:,i) = arr_val(:,1);
    end
    
    if nbr~=0
        inz  = find(arr_inz_airfoil<=nbr);
        temp_coo = arr_xhi_airf_c(inz(end)+1,:);
        arr_xhi_airf_w(inz(end)+1,:) = [temp_coo(1:end),temp_coo(end-1:-1:1)];
    end
    
%     figure(); hold on; grid on;
%     for i = 1:size(arr_xhi_airf_w,1)
%         plot(arr_xhi_y_w(:,1),arr_xhi_airf_w(i,:));
%     end
    
    uvlm_obj.arr_xhi_airf_c  = arr_xhi_airf_c;
    uvlm_obj.arr_xhi_airf_w  = arr_xhi_airf_w;
    uvlm_obj.xhi_y_w         = xhi_y_w;
    uvlm_obj.xhi_y_c         = xhi_y_c;
    uvlm_obj.arr_xhi_x       = xhi_x;
    uvlm_obj.nbr             = nbr;
    uvlm_obj.arr_inz_airfoil = arr_inz_airfoil;
    
return
