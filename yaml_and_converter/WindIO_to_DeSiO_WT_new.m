% =================================================================================================================
% Function to transform WindIO to DeSiO-input for modeling wind energy
% converter.
% 
% Author: Christian Hente
% Date: 05.05.2022
% updated: 27.10.2022
% =================================================================================================================
function WindIO_to_DeSiO_WT()
% =================================================================================================================
clc; clear all; close all; format short;
i1 = [1,0,0]; i2 = [0,1,0]; i3 = [0,0,1];

% path for necessary functions
addpath(genpath('Pre_post_processing'));
currDir = cd;

%% =================================================================================================================
    path_DeSiO = 'C:\Users\hente\Desktop\DeSiO_compile\DeSiO\DeSiO\x64\Release';
    path_mtee  = 'C:\Users\hente\AppData\Local\Programs\DeSiO';
    %
    str_dir_yaml = 'IEA-15-240-RWT_new.yaml';                             % filename of .yaml input file
% =================================================================================================================
    scaling_blade = 1; 
    scaling_tower = 1;

    cos_msl     = [0,0,0];                                                 % MSL position in global coordinate system
    
    n_yaw       = [0,0,1];                                                 % yaw rotation axis
    n_tilt      = [0,1,0];                                                 % tilt rotation axis
    
    flag_tower       = 1;                                                  % - tower on(1)/off(0)
    flag_tower_aero  = 0;                                                  % - tower aero grid on(1)/off(0)
    flag_tower_struc = 0;                                                  % - tower aero grid on(1)/off(0)

    flag_blade       = 0;                                                  % - blades on(1)/off(0)
    flag_blade_aero  = 0;                                                  % - blades aero grid on(1)/off(0)
    flag_blade_struc = 0;                                                  % - blades structure mesh on(1)/off(0)
    
    flag_foundation       = 0;                                             % - foundation on(1)/off(0)
    flag_foundation_aero  = 0;                                             % - foundation aero grid on(1)/off(0)
    flag_foundation_struc = 0;                                             % - foundation structure mesh on(1)/off(0)
    
    flag_hub     = 1;                                                      % - hub on(1)/off(0)
    flag_nacelle = 1;                                                      % - nacelle on(1)/off(0)
% =================================================================================================================
%% =================================================================================================================
%                   INITIALIZING OBJECTS
% =================================================================================================================
    simu_struct = []; simu_aero = []; simu_fsi = []; 
    mesh_struc  = []; mesh_aero = []; simu_fsi.data = {}; simu_fsi.data_input = []; wake = [];
% =================================================================================================================
%% =================================================================================================================
%                   READING WINDIO-FORMAT and DEFINING WT-OBJECTS
% =================================================================================================================
    model = YAML.read(str_dir_yaml);                                       % read .yaml file
    if isfield(model,'jobname'); strfilename  = model.jobname; end
    strfilename = [strfilename '_' 'pitch' ...
            num2str(model.environment.pitch_angle)...
            '_vel' ...
            num2str(model.environment.vinf) ...
            '_nspan' num2str(model.components.blade.DeSiO.uvlm.M_aero)];
    strfilename = [strfilename '_test'];
% =================================================================================================================
% enviroment conditions
% =================================================================================================================
    yaw_angle = 0.0;
    pitch_angle = 0.0;
    if isfield(model,'environment')
       if isfield(model.environment,'yaw_angle'); yaw_angle = model.environment.yaw_angle; end;
       if isfield(model.environment,'pitch_angle'); pitch_angle = model.environment.pitch_angle; end;
    end
    
    fsi_radius_rbf = 0;
    if isfield(model.simulationparamter,'fsi_radius_rbf')
        fsi_radius_rbf = model.simulationparamter.fsi_radius_rbf;
    end
        
% =================================================================================================================
% blades
% =================================================================================================================
    n_blades = 0;
    if isfield(model.assembly,'number_of_blades')
        n_blades = model.assembly.number_of_blades;
    end
    if isfield(model.components,'blade') && flag_blade == 1
        blade_Aero  = fun_get_uvlm_geometry('blade',model.components.blade.DeSiO.uvlm,model.airfoils,scaling_blade);
        blade_Struc = fun_get_beam_model('blade',model.components.blade.DeSiO.beam,model.materials,scaling_blade);
        wake_diffusion = 0.0;
        if isfield(model.simulationparamter,'wake_diffusion'); 
            wake_diffusion = model.simulationparamter.wake_diffusion;
        end
    end
    
% =================================================================================================================
% tower
% =================================================================================================================
    tower_top = 0;
    tower_bot = 0;
    if isfield(model.components,'tower') && flag_tower == 1
        tower_Aero   = fun_get_uvlm_geometry('pipe',model.components.tower.DeSiO.uvlm,model.airfoils,scaling_tower);    % create 2d mesh for aero grid
        tower_Struc  = fun_get_beam_model('pipe',model.components.tower.DeSiO.beam,model.materials,scaling_tower);      % create structural mesh
        tower_top    = max(model.components.tower.DeSiO.beam.reference_axis.z.values)*scaling_tower;
        tower_bot    = min(model.components.tower.DeSiO.beam.reference_axis.z.values)*scaling_tower;
    end
    
% =================================================================================================================
% monopile
% =================================================================================================================
    tp_mass = 0;
    if isfield(model.components,'monopile')
        monopile_Aero  = fun_get_uvlm_geometry('pipe',model.components.monopile.DeSiO.uvlm,model.airfoils,scaling_tower);
        monopile_Struc = fun_get_beam_model('pipe',model.components.monopile.DeSiO.beam,model.materials,scaling_tower);
        tp_mass        = model.components.monopile.transition_piece_mass;
    end
    
% =================================================================================================================
% nacelle
% =================================================================================================================
    uptilt_angle       = 0.0;
    overhang           = 0.0;
    Twr2Shft           = 0.0;
    e_shaft            = [0,0,0];
    nacelle_centerm_tt = [0,0,0];
    yaw_mass           = 0.0;
    yaw_center         = [0;0;0];
    if isfield(model.components,'nacelle') && flag_nacelle == 1
        uptilt_angle       = model.components.nacelle.DeSiO.drivetrain.uptilt; % in rad
        overhang           = model.components.nacelle.DeSiO.drivetrain.overhang;
        Twr2Shft           = model.components.nacelle.DeSiO.drivetrain.Twr2Shft;
        yaw_mass           = model.components.nacelle.DeSiO.drivetrain.yaw_mass;
        e_shaft            = -[(cos(uptilt_angle)*eye(3) + sin(uptilt_angle)*skew(n_tilt)+(1-cos(uptilt_angle))*n_tilt'*n_tilt)*i1']';
        nacelle_centerm_tt = model.components.nacelle.DeSiO.elastic_properties_mb.center_mass;
    end
    shaft_rb1_tt = nacelle_centerm_tt;
    
% =================================================================================================================
% hub
% =================================================================================================================
    hub_centerm_tt = [0,0,0];
    hub_2Apex_tt   = 0.0;
    hub_diameter   = 0.0;
    if isfield(model.components,'hub') && flag_hub == 1
        hub_diameter   = model.components.hub.DeSiO.diameter;
        hub_2Apex_tt   = model.components.hub.DeSiO.Hub2Apex;
        hub_centerm_tt = model.components.hub.DeSiO.elastic_properties_mb.center_mass;
    end
    shaft_rb2_tt = hub_centerm_tt;
    
    if isfield(model.components,'nacelle') && flag_nacelle == 1
        if isfield(model.components,'hub') && flag_hub == 1
            shaft_rb1_tt   = [0.0, 0.0, Twr2Shft];
            shaft_rb2_tt   = shaft_rb1_tt + e_shaft*overhang;
            hub_centerm_tt = shaft_rb1_tt + e_shaft*(overhang+hub_2Apex_tt);
        end
    end
% =================================================================================================================
%% ================================================================================================================
%                       DEFINING COORDAINTES IN GLOBAL COS (0.0m - MSL)
% =================================================================================================================
    cos_tt = cos_msl + [0,0,tower_top];                                    % global coordinates of tower top
    cos_tb = cos_msl + [0,0,tower_bot];                                    % global coordinates of tower bottom
    cos_hc = cos_tt + hub_centerm_tt;                                      % global coordinates of hub center of mass
    cos_na = cos_tt + nacelle_centerm_tt;                                  % global coordinates nacelle center of mass
% =================================================================================================================
%% =================================================================================================================
%                           NACELLE AND HUB
% =================================================================================================================
    nrb = 0; nhub = 0; nnac = 0;
    if isfield(model.components,'nacelle') && flag_nacelle == 1
        nrb  = nrb + 1;
        nnac = 1;
        simu_struct.rb(nnac).type = 'nacelle';
        simu_struct.rb(nnac).center_mass = cos_na;
        simu_struct.rb(nnac).D1 = [1,0,0]; 
        simu_struct.rb(nnac).D2 = [0,1,0];
        simu_struct.rb(nnac).D3 = [0,0,1];
        simu_struct.rb(nnac).mass_matrix = zeros(1,10);
        if isfield(model.components.nacelle.DeSiO.elastic_properties_mb,'mass_matrix')
            simu_struct.rb(nnac).mass_matrix = model.components.nacelle.DeSiO.elastic_properties_mb.mass_matrix;
        end
    end
% =================================================================================================================    
    if isfield(model.components,'hub') && flag_hub == 1
        nrb  = nrb + 1;
        nhub = nnac + 1;
        
        alpha_hub = 45;
        R_hub     = cos(alpha_hub*pi/180)*eye(3) + sin(alpha_hub*pi/180)*skew(n_yaw)+(1-cos(alpha_hub*pi/180))*n_yaw'*n_yaw;

        simu_struct.rb(nhub).type = 'hub';
        simu_struct.rb(nhub).center_mass = cos_hc;
        simu_struct.rb(nhub).D1 = (R_hub*[1;0;0])';
        simu_struct.rb(nhub).D2 = (R_hub*[0;1;0])';
        simu_struct.rb(nhub).D3 = (R_hub*[0;0;1])';
        simu_struct.rb(nhub).mass_matrix = zeros(1,10);
        if isfield(model.components.hub.DeSiO.elastic_properties_mb,'mass_matrix')
            simu_struct.rb(nhub).mass_matrix = model.components.hub.DeSiO.elastic_properties_mb.mass_matrix;
        end
    end
% =================================================================================================================
%% ================================================================================================================
%                                 FOUNDATION STRUCTURE (right now only monopile)
% =================================================================================================================
    imat0    = 0; % couner for material
    nn12     = 0; % global beam node12 counter
    n_struc  = size(mesh_struc,2); 
    n_surf   = size(mesh_aero,2);
    monopile = [];
    if isfield(model.components,'monopile') && flag_foundation == 1
        % creating structural mesh in DeSiO-Format
        if flag_foundation_struc == 1; 
            X_R = [ones(monopile_Struc.M+1,1)*cos_msl(1),ones(monopile_Struc.M+1,1)*cos_msl(2),ones(monopile_Struc.M+1,1)*cos_msl(3)];
            monopile.grid.struc.X_RE         = [monopile_Struc.arr_coordinates(:,1:3)] + X_R;
            monopile.grid.struc.D1           = [monopile_Struc.arr_coordinates(:,4:6)];
            monopile.grid.struc.D2           = [monopile_Struc.arr_coordinates(:,7:9)];
            monopile.grid.struc.D3           = [monopile_Struc.arr_coordinates(:,10:12)];
            monopile.grid.struc.M            = monopile_Struc.M;
            monopile.grid.struc.connectivity = monopile_Struc.connectivity;

            monopile.grid.struc.matBeam.n_gl = nrb + nn12 + [1:monopile_Struc.M + 1];
            monopile.grid.struc.matBeam.inz  = [1:monopile_Struc.M];
            monopile.grid.struc.matBeam.Cmat = monopile_Struc.arr_stiff_matrix;
            monopile.grid.struc.matBeam.mmat = monopile_Struc.arr_mass_matrix;
            monopile.grid.struc.matBeam.diss = monopile_Struc.dissipation;

            mesh_struc                             = fun_set_struc_mesh(mesh_struc,monopile);
            mesh_struc(size(mesh_struc,2)).strname = 'monopile';
            imat0                                  = mesh_struc(size(mesh_struc,2)).imat;
            n_struc                                = n_struc + 1;
            nn12                                   = monopile.grid.struc.matBeam.n_gl(end)-nrb;
        end
        
        % creating aero grid in DeSiO-Format
        if flag_foundation_aero == 1; 
            X_R                                = [ones((monopile_Aero.M+1)*(2*monopile_Aero.N+1),1)*cos_msl(1),ones((monopile_Aero.M+1)*(2*monopile_Aero.N+1),1)*cos_msl(2),ones((monopile_Aero.M+1)*(2*monopile_Aero.N+1),1)*cos_msl(3)];
            monopile.grid.aero.X               = monopile_Aero.X_W + X_R;
            monopile.grid.aero.M               = monopile_Aero.M; monopile(1).grid.aero.N = 2*monopile_Aero.N;
            monopile.grid.aero.connectivity    = monopile_Aero.connectivity_w;
            monopile.grid.aero.node_fsi_radius = monopile_Aero.arr_node_fsi_radius;
            
            mesh_aero = fun_set_aero_mesh(mesh_aero,monopile); 
            n_surf    = n_surf + 1;
            fsys      = figure(); hold on; grid on; axis equal; view(3); title('monopile');
            fun_plot_3Dmesh(fsys,'2D',mesh_aero(1));
        end
        
        % set input for FSI
        if flag_foundation_aero == 1 && flag_foundation_struc == 1
            if isempty(monopile.grid.aero.node_fsi_radius)
                simu_fsi.data(end+1,:) = {'beam',n_struc,n_surf, fsi_radius_rbf};
            else
                fsi_radius_filename = ['file_fsi_radius_' num2str(n_surf) '_input.txt'];
                simu_fsi.data(end+1,:) = {'beam',n_struc,n_surf, fsi_radius_filename};
                simu_fsi.data_input(end+1).node_fsi_radius = monopile.grid.aero.node_fsi_radius;
            end
        end
        
    end
% =================================================================================================================
%% ================================================================================================================
%                                 TOWER
% =================================================================================================================
    n_struc = size(mesh_struc,2); 
    n_surf  = size(mesh_aero,2);
    tower = [];
    if isfield(model.components,'tower') && flag_tower == 1
        % creating structural mesh in DeSiO-Format
        if flag_tower_struc == 1; 
            X_R                           = [ones(tower_Struc.M+1,1)*cos_msl(1),ones(tower_Struc.M+1,1)*cos_msl(2),ones(tower_Struc.M+1,1)*cos_msl(3)];
            tower.grid.struc.X_RE         = [tower_Struc.arr_coordinates(:,1:3)] + X_R;
            tower.grid.struc.D1           = [tower_Struc.arr_coordinates(:,4:6)];
            tower.grid.struc.D2           = [tower_Struc.arr_coordinates(:,7:9)];
            tower.grid.struc.D3           = [tower_Struc.arr_coordinates(:,10:12)];
            tower.grid.struc.M            = tower_Struc.M;
            tower.grid.struc.connectivity = tower_Struc.connectivity;

            tower.grid.struc.matBeam.n_gl = nrb + nn12 + [1:tower_Struc.M + 1];
            tower.grid.struc.matBeam.inz  = [1:tower_Struc.M];
            tower.grid.struc.matBeam.Cmat = tower_Struc.arr_stiff_matrix;
            tower.grid.struc.matBeam.mmat = tower_Struc.arr_mass_matrix;
            tower.grid.struc.matBeam.diss = tower_Struc.dissipation;
            
            mesh_struc                             = fun_set_struc_mesh(mesh_struc,tower);
            mesh_struc(size(mesh_struc,2)).strname = 'tower';
            imat0                                  = mesh_struc(size(mesh_struc,2)).imat;
            n_struc                                = n_struc + 1;
            nn12                                   = tower.grid.struc.matBeam.n_gl(end)-nrb;
        end
        
        % creating aero grid in DeSiO-Format
        if flag_tower_aero == 1; 
            X_R                             = [ones((tower_Aero.M+1)*(2*tower_Aero.N+1),1)*cos_msl(1),ones((tower_Aero.M+1)*(2*tower_Aero.N+1),1)*cos_msl(2),ones((tower_Aero.M+1)*(2*tower_Aero.N+1),1)*cos_msl(3)];
            tower.grid.aero.X               = tower_Aero.X_W + X_R;
            tower.grid.aero.M               = tower_Aero.M; tower(1).grid.aero.N = 2*tower_Aero.N;
            tower.grid.aero.connectivity    = tower_Aero.connectivity_w;
            tower.grid.aero.node_fsi_radius = tower_Aero.arr_node_fsi_radius;
            
            mesh_aero = fun_set_aero_mesh(mesh_aero,tower); 
            n_surf    = n_surf + 1;
            fsys      = figure(); hold on; grid on; axis equal; view(3); title('tower');
            fun_plot_3Dmesh(fsys,'2D',mesh_aero(n_surf));
        end
        
        % set input for FSI
        if flag_tower_aero == 1 && flag_tower_struc == 1
            if isempty(tower.grid.aero.node_fsi_radius)
                simu_fsi.data(end+1,:) = {'beam',n_struc,n_surf, fsi_radius_rbf};
            else
                fsi_radius_filename = ['file_fsi_radius_' num2str(n_surf) '_input.txt'];
                simu_fsi.data(end+1,:) = {'beam',n_struc,n_surf, fsi_radius_filename};
                simu_fsi.data_input(end+1).node_fsi_radius = tower.grid.aero.node_fsi_radius;
            end            
        end
        
    end
% =================================================================================================================
%% ================================================================================================================
%                                 BLADES
% =================================================================================================================
    nbl = 0; imatb = imat0;
    blades = [];
    if isfield(model.components,'blade') && flag_blade == 1
        ir1    = -i2; 
        ir2    = i1;
        ir3    = i3;
        
        % Blade root and rest of the blade
        nn = 0; 
        ne = 0; 
        inz_blade_root = [];
        if blade_Aero.nbr~=0
            % getting indices for extracting blade root from arrays
            arr_inz_airfoil = blade_Aero.arr_inz_airfoil;
            inz = find(arr_inz_airfoil<=blade_Aero.nbr);
            inz_blade_root = [1:inz(end)+1];

            M_br = length(inz_blade_root)-1;
            N_br = 2*blade_Aero.N;
           
            nn_br = (M_br+1)*(blade_Aero.N+1);
            ne_br = M_br*blade_Aero.N;
        end
        
        for i = 1:n_blades
            n_struc = size(mesh_struc,2);
            n_surf  = size(mesh_aero,2);
            
            % position angle of blade i in rotor
            alpha_blades = (i-1)*2*pi/n_blades*180/pi;
            
            % positioning blade in rotor
            blades(i).grid = fun_blades2rotor(blade_Aero, blade_Struc, ir3, -pitch_angle, ir2, alpha_blades, n_tilt, uptilt_angle*180/pi, [0,0,hub_diameter/2]); 

            % creating structural mesh in DeSiO-Format
            if flag_blade_struc == 1;
                nbl = blade_Struc.M+1;
                X_R = [ones(blade_Struc.M+1,1)*cos_hc(1),ones(blade_Struc.M+1,1)*cos_hc(2),ones(blade_Struc.M+1,1)*cos_hc(3)];
                
                blades(i).grid.struc.X_RE = blades(i).grid.struc.X_RE + X_R;
                blades(i).grid.struc.M = blade_Struc.M;
                blades(i).grid.struc.connectivity = blade_Struc.connectivity;
                blades(i).grid.struc.matBeam.n_gl = nrb + nn12 + (i-1)*nbl + [1:nbl];
                blades(i).grid.struc.matBeam.inz  = [1:blade_Struc.M];
                blades(i).grid.struc.matBeam.Cmat = blade_Struc.arr_stiff_matrix;
                blades(i).grid.struc.matBeam.mmat = blade_Struc.arr_mass_matrix;
                blades(i).grid.struc.matBeam.diss = blade_Struc.dissipation;

                mesh_struc = fun_set_struc_mesh(mesh_struc,blades(i));
                mesh_struc(size(mesh_struc,2)).strname = ['blade ' num2str(i)];
                mesh_struc(size(mesh_struc,2)).imat = imatb;
                
                n_struc = n_struc + 1;
                nbl = n_blades*nbl;
            end
            
            % creating aero grid in DeSiO-Format
            if flag_blade_aero == 1; 
                X_R = [ones((blade_Aero.M+1)*(blade_Aero.N+1),1)*cos_hc(1),ones((blade_Aero.M+1)*(blade_Aero.N+1),1)*cos_hc(2),ones((blade_Aero.M+1)*(blade_Aero.N+1),1)*cos_hc(3)];
                X_C = blades(i).grid.aero.X_C  + X_R;

                X_R = [ones((blade_Aero.M+1)*(2*blade_Aero.N+1),1)*cos_hc(1),ones((blade_Aero.M+1)*(2*blade_Aero.N+1),1)*cos_hc(2),ones((blade_Aero.M+1)*(2*blade_Aero.N+1),1)*cos_hc(3)];
                X_W = blades(i).grid.aero.X_W  + X_R;

                % Blade root and rest of the blade
                if blade_Aero.nbr~=0
                    % extracting blade root from arrays
                    blade_obj(1).grid.aero.N = N_br;
                    blade_obj(1).grid.aero.M = M_br;
                    blade_obj(1).grid.aero.X = X_W([1:(N_br+1)*(M_br+1)],:);
                    blade_obj(1).grid.aero.node_fsi_radius = blade_Aero.arr_node_fsi_radius_w([1:(N_br+1)*(M_br+1)]);
                    blade_obj(1).grid.aero.connectivity = blade_Aero.connectivity_w([1:M_br*N_br],:);

                    % extracting rest of blade from arrays
                    blade_obj(2).grid.aero.N = blade_Aero.N;
                    blade_obj(2).grid.aero.M = blade_Aero.M-M_br;
                    blade_obj(2).grid.aero.X = X_C(nn_br-blade_Aero.N:end,:);
                    blade_obj(2).grid.aero.node_fsi_radius = blade_Aero.arr_node_fsi_radius(nn_br-blade_Aero.N:end);
                    blade_obj(2).grid.aero.connectivity = blade_Aero.connectivity_c(ne_br+1:end,:)-(nn_br-blade_Aero.N)+1;
                else
                    blade_obj(1).grid.aero.N = blade_Aero.N;
                    blade_obj(1).grid.aero.M = blade_Aero.M;
                    blade_obj(1).grid.aero.X = X_C;
                    blade_obj(1).grid.aero.node_fsi_radius = blade_Aero.arr_node_fsi_radius;
                    blade_obj(1).grid.aero.connectivity = blade_Aero.connectivity_c;
                end
            
                % creating mesh in DeSiO-Format
                mesh_aero = fun_set_aero_mesh(mesh_aero,blade_obj);
                
                if i == 1; fsys = figure(); hold on; grid on; axis equal; title('blades'); end;
                for j = 1:size(blade_obj,2)
                    fun_plot_3Dmesh(fsys,'2D',mesh_aero(n_surf+j));
                end
                wake = fun_blade_wake(wake,n_surf + size(blade_obj,2),mesh_aero(n_surf + size(blade_obj,2)),model.simulationparamter.nwakerows,model.simulationparamter.nwakerows,wake_diffusion);
            end
            
            % setting input for FSI
            if flag_blade_aero == 1 && flag_blade_struc == 1
                for j = 1:size(blade_obj,2)
                    if isempty(blade_obj(j).grid.aero.node_fsi_radius)
                        simu_fsi.data(end+1,:) = {'beam',n_struc,n_surf+j, fsi_radius_rbf};
                    else
                        fsi_radius_filename = ['file_fsi_radius_' num2str(n_surf+j) '_input.txt'];
                        simu_fsi.data(end+1,:) = {'beam',n_struc,n_surf+j, fsi_radius_filename};
                        simu_fsi.data_input(end+1).node_fsi_radius = blade_obj(j).grid.aero.node_fsi_radius;
                    end
                end
            end
        end
    end
% =================================================================================================================
%% ================================================================================================================
%                       ADDITIONAL MASSES
% =================================================================================================================
    % additional mass for transition piece at monopile top
    nam = 0;
    if not(isempty(monopile)) && flag_foundation_struc == 1
        nam = nam + 1;
        simu_struct.pointmass12(nam).mass = tp_mass;
        simu_struct.pointmass12(nam).node = monopile.grid.struc.matBeam.n_gl(end);
    end
    % additional mass for yaw bearing at tower top
    if isfield(model.components,'nacelle') && flag_nacelle == 1
        if not(isempty(tower)) && flag_tower_struc == 1
            nam = nam + 1;
            simu_struct.pointmass12(nam).mass = yaw_mass;
            simu_struct.pointmass12(nam).node = tower.grid.struc.matBeam.n_gl(end);
        end
    end
% =================================================================================================================
%% ================================================================================================================
%                           CONSTRAINTS
% =================================================================================================================
    simu_struct.constraints = [];
    remo_internal_const = [];
    nco = 0;
    if not(isempty(tower)) && flag_tower_struc == 1
         if not(isempty(monopile)) && flag_foundation_struc == 1
            nco = nco + 1;
            simu_struct.constraints(nco).type  = 'rigidsupport';
            simu_struct.constraints(nco).nodes = [monopile.grid.struc.matBeam.n_gl(1) 0];
            remo_internal_const(end+1) = monopile.grid.struc.matBeam.n_gl(1);
            
            nco = nco + 1;
            simu_struct.constraints(nco).type  = 'rigidconnection';
            simu_struct.constraints(nco).nodes = [monopile.grid.struc.matBeam.n_gl(end) tower.grid.struc.matBeam.n_gl(1)];
            remo_internal_const(end+1) = monopile.grid.struc.matBeam.n_gl(end);
        else
            nco = nco + 1;
            simu_struct.constraints(nco).type  = 'rigidsupport';
            simu_struct.constraints(nco).nodes = [tower.grid.struc.matBeam.n_gl(1) 0];
            remo_internal_const(end+1) = tower.grid.struc.matBeam.n_gl(1);
        end
        if isfield(model.components,'nacelle') && flag_nacelle == 1
            nco = nco + 1;
            simu_struct.constraints(nco).type  = 'rigidconnection';
            simu_struct.constraints(nco).nodes = [tower.grid.struc.matBeam.n_gl(end) nnac];
            remo_internal_const(end+1) = tower.grid.struc.matBeam.n_gl(end);
        end
    end
    if isfield(model.components,'hub') && flag_hub == 1
        if isfield(model.components,'nacelle') && flag_nacelle == 1
            nco  = nco + 1;
            nhub = nnac + 1;
            
            simu_struct.constraints(nco).type  = 'revolutejoint';
            simu_struct.constraints(nco).nodes = [nnac nhub];
            
            phi1 = [cos_tt + shaft_rb1_tt] - simu_struct.rb(nnac).center_mass ;
            simu_struct.constraints(nco).phi1 = [phi1*simu_struct.rb(nnac).D1', phi1*simu_struct.rb(nnac).D2', phi1*simu_struct.rb(nnac).D3'];

            phi2 = [cos_tt + shaft_rb2_tt] - simu_struct.rb(nhub).center_mass;
            simu_struct.constraints(nco).phi2 = [phi2*simu_struct.rb(nhub).D1', phi2*simu_struct.rb(nhub).D2', phi2*simu_struct.rb(nhub).D3'];
            
            % rotor rotation axis in global cos and local hub cos:
            dir_gl = simu_struct.rb(nhub).center_mass - (simu_struct.rb(nnac).center_mass+phi1);
            dir_gl = dir_gl/norm(dir_gl);
            dir_lo = [dir_gl*simu_struct.rb(nhub).D1', dir_gl*simu_struct.rb(nhub).D2', dir_gl*simu_struct.rb(nhub).D3'];
        end
    end
    
    if not(isempty(blades)) && flag_blade_struc == 1
        if isfield(model.components,'hub') && flag_hub == 1
            nco0 = nco;
            for i = 1:n_blades
                nco  = nco0 + i;
                nhub = nnac + 1;
                
                simu_struct.constraints(nco).type  = 'rigidconnection';
                simu_struct.constraints(nco).nodes = [nhub blades(i).grid.struc.matBeam.n_gl(1)];
                simu_struct.constraints(nco).phi1  = [0,0,0];

                remo_internal_const(end+1) = blades(i).grid.struc.matBeam.n_gl(1);
            end
            if (isempty(tower)) || flag_tower_struc == 0
                nco = nco + 1;
                simu_struct.constraints(nco).type  = 'sphericalsupport';
                simu_struct.constraints(nco).nodes = [nhub 0];
                nco = nco + 1;
                simu_struct.constraints(nco).type  = 'rotation_global';
                simu_struct.constraints(nco).nodes = [nhub 0];
                simu_struct.constraints(nco).dir   = [0.00,1.00,0.00];
                nco = nco + 1;
                simu_struct.constraints(nco).type  = 'rotation_global';
                simu_struct.constraints(nco).nodes = [nhub 0];
                simu_struct.constraints(nco).dir   = [0.00,0.00,1.00];
            end
        end
    end
    
    nn = nrb + nn12 + nbl; nodes = 1:nn;
    ai = size(simu_struct.constraints,2);
    for i = 1:nn
        if all(remo_internal_const-nodes(i))
            ai = ai + 1;
            simu_struct.constraints(ai).type  = 'internal';
            simu_struct.constraints(ai).nodes = [nodes(i) 0];
        end
    end
% =================================================================================================================
%% ================================================================================================================
%                           SIMULATION SETTINGS
% =================================================================================================================
    simu_aero.currDir     = currDir;
    simu_aero.strfilename = [strfilename];
    simu_aero.time        = model.simulationparamter.time;
    simu_aero.deltat      = model.simulationparamter.deltat_aero;
    simu_aero.cutoff      = model.simulationparamter.cutoff;
    simu_aero.density     = model.environment.air_density;
    simu_aero.i_vinf      = model.environment.vinf;
    simu_aero.d_vinf      = model.environment.winddir;
    
    simu_struct.currDir     = currDir;
    simu_struct.strfilename = [strfilename];
    simu_struct.time        = model.simulationparamter.time;
    simu_struct.deltat      = model.simulationparamter.deltat_struc;
    simu_struct.tol         = model.simulationparamter.tolerance;
    simu_struct.niter       = model.simulationparamter.niter;
    simu_struct.grav        = model.simulationparamter.gravity;
    simu_struct.grav_vec    = model.simulationparamter.grav_vector;
    simu_struct.loads       = [];
    
    simu_fsi.lf_type = 'constant';
    simu_fsi.lf_duration = 10.0;
    
    % flags for pre-calculations
    flag_init_self_weight = 0;
    if isfield(model.simulationparamter,'flag_init_self_weight') 
        flag_init_self_weight = model.simulationparamter.flag_init_self_weight;
    end
    simu_struct.grav = flag_init_self_weight;
    
    flag_init_rotor_velocity = 0;
    if isfield(model.simulationparamter,'flag_init_rotor_velocity')
        flag_init_rotor_velocity = model.simulationparamter.flag_init_rotor_velocity;
    end
    
    % simulation settings regarding DeSiO executeable
    os = {'windows'};
    if isfield(model.simulationparamter,'os')
        os = model.simulationparamter.os;
    end
    
    path_DeSiO = {''};
    if isfield(model.simulationparamter,'path_DeSiO')
        path_DeSiO = model.simulationparamter.path_DeSiO;
    end
    
    ifort_version = {''};
    if isfield(model.simulationparamter,'ifort_version')
        ifort_version = model.simulationparamter.ifort_version;
    end
    
    % fsi settings
    simu_fsi.currDir     = currDir;
    simu_fsi.strfilename = [strfilename];
    
    simu_fsi.radius_rbf  = 0;
    if isfield(model.simulationparamter,'fsi_radius_rbf')
        simu_fsi.radius_rbf  = model.simulationparamter.fsi_radius_rbf;
    end
    
    % external flow field data
    simu_aero.sort  = 'constant';
    if isfield(model.environment,'wind_field_file')
        wind_field_file = model.environment.wind_field_file;
        simu_aero.sort  = 'file';
        simu_aero.wind_field_file = wind_field_file;
        simu_aero.grid_center = simu_struct.rb(nhub).center_mass;
    end
    
    % beam element settings for elasticity matrix
    simu_struct.flag_matbeam = 1.0;
    imat = 0;
    if not(isempty(monopile)) && flag_foundation_struc == 1
        for i = 1:size(monopile.grid.struc.matBeam.Cmat,1)
            simu_struct.matbeam(i).cmat    = monopile.grid.struc.matBeam.Cmat(i,:);
            simu_struct.matbeam(i).mmat    = monopile.grid.struc.matBeam.mmat(i,:);
            simu_struct.matbeam(i).diss    = monopile.grid.struc.matBeam.diss;
            simu_struct.matbeam(i).strname = 'monopile';
        end
        imat = i;
    end
    if not(isempty(tower)) && flag_tower_struc == 1
        for i = 1:size(tower.grid.struc.matBeam.Cmat,1)
            simu_struct.matbeam(i+imat).cmat    = tower.grid.struc.matBeam.Cmat(i,:);
            simu_struct.matbeam(i+imat).mmat    = tower.grid.struc.matBeam.mmat(i,:);
            simu_struct.matbeam(i+imat).diss    = tower.grid.struc.matBeam.diss;
            simu_struct.matbeam(i+imat).strname = 'tower';
        end
        imat = i + imat;
    end
    if not(isempty(blades)) && flag_blade_struc == 1      
        for i = 1:size(blades(1).grid.struc.matBeam.Cmat,1)
            simu_struct.matbeam(i+imat).cmat = blades(1).grid.struc.matBeam.Cmat(i,:);
            simu_struct.matbeam(i+imat).mmat = blades(1).grid.struc.matBeam.mmat(i,:);
            simu_struct.matbeam(i+imat).diss = blades(1).grid.struc.matBeam.diss;
            simu_struct.matbeam(i+imat).strname = 'blade';
        end
    end
% =================================================================================================================
%% ================================================================================================================
%                           WRITING DESIO-INPUT FILES
% =================================================================================================================
    if size(simu_aero,2) ~= 0 && size(mesh_aero,2) ~= 0
        % write DeSiO-Aero input
        fun_writeDeSiOAeroInput(simu_aero,mesh_aero,wake);
    end
    
    if size(simu_struct,2) ~= 0 && size(mesh_struc,2) ~= 0
        % initial DeSiO-Structure input
        path_struc = [simu_struct.currDir '\' simu_struct.strfilename '\DeSiO-Structure'];
        simu_struct.jobname  = 'initial_data';
        fun_writeDeSiOStructureInput(simu_struct,mesh_struc);
        
        % writing structural input file for static pre-analysis for self-weight
        simu_struct_static = simu_struct;
        simu_struct_static.type     = 'static';
        simu_struct_static.jobname  = 'self_weight';
        simu_struct_static.settings = [1.0 1.0e-1 1.0e-6 50 1];
        % adding constraint to fix rotor-axis
        simu_struct_static.constraints(end+1).type  = 'rotation_local';
        simu_struct_static.constraints(end).nodes = [nhub 0];
        simu_struct_static.constraints(end).dir   = [dir_lo];
        fun_writeDeSiOStructureInput(simu_struct_static,mesh_struc);
        
        % writing structural input file for dynamic pre-analysis for initial rotor speed
        simu_struct_dynamic = simu_struct;
        if flag_init_rotor_velocity == 1
            simu_struct_dynamic.type     = 'dynamic';
            simu_struct_dynamic.jobname  = 'init_rotor_speed';
            if flag_init_self_weight == 1
                simu_struct_dynamic.settings = [10.0 1.0e-1 1.0e-6 50 1];
            else
                simu_struct_dynamic.settings = [10.0 1.0e-1 1.0e-6 50 0];
            end
            % asign simulation parameter
            simu_struct_dynamic.constraints(end+1).type  = 'angularvelocity_local';
            simu_struct_dynamic.constraints(end).nodes   = [nhub 0];
            simu_struct_dynamic.constraints(end).dir     = [dir_lo];
            init_rotor_velocity = 0;
            if isfield(model.environment,'init_rotor_velocity'); 
                init_rotor_velocity = -model.environment.init_rotor_velocity*2*pi/60;
            end        
            simu_struct_dynamic.boundary12(1).sort         = 'linearC';
            simu_struct_dynamic.boundary12(1).intensity    = init_rotor_velocity;
            simu_struct_dynamic.boundary12(1).time         = 10.0;
            simu_struct_dynamic.boundary12(1).constraintID = size(simu_struct_dynamic.constraints,2);
            simu_struct_dynamic.boundary12(1).file         = 'none';
            % create files but statement to ask if self weight is on or not
            if flag_init_self_weight == 1
                simu_struct_dynamic.prevjobname = simu_struct_static.jobname;
            end
            fun_writeDeSiOStructureInput(simu_struct_dynamic,mesh_struc);
        end
        
        % writing structural input file for modal analysis
        simu_struct_modal = simu_struct;
        simu_struct_modal.type     = 'modal';
        simu_struct_modal.jobname  = 'modalanalysis';
        simu_struct_modal.settings = [10 1.0e-6 0.0 0.0 0];
        fun_writeDeSiOStructureInput(simu_struct_modal,mesh_struc);

    end
    if size(simu_fsi,2) ~= 0 
        if size(simu_aero,2) ~= 0 && size(simu_struct,2) ~= 0
            path_fsi   = [simu_fsi.currDir '\' simu_fsi.strfilename '\DeSiO'];
            path_aero  = [simu_aero.currDir '\' simu_aero.strfilename '\DeSiO-Aero'];
            path_struc = [simu_struct.currDir '\' simu_struct.strfilename '\DeSiO-Structure'];
            mkdir(path_fsi);
            filePattern = fullfile(path_aero,'*.txt'); dd = dir(filePattern);
            for j = 1:size(dd,1)
                if isempty(strfind(dd(j).name,'simulationinput'))
                    copyfile([path_aero '\' dd(j).name], path_fsi);
                end
            end
            filePattern = fullfile([path_struc '\initial_data'],'*.txt'); dd = dir(filePattern);
            for j = 1:size(dd,1)
                if isempty(strfind(dd(j).name,'simulationinput'))
                    copyfile([path_struc '\initial_data\' dd(j).name], path_fsi);
                end
            end
        end
        
        % create fsi input files
        if flag_init_rotor_velocity == 1 % set jobname for fsi calculation
            simu_fsi.structprevjobname = simu_struct_dynamic.jobname;
            simu_fsi.aeroprevjobname  = 'none';
        elseif flag_init_rotor_velocity == 0 && flag_init_self_weight == 1
            simu_fsi.structprevjobname = simu_struct_static.jobname;
            simu_fsi.aeroprevjobname  = 'none';
        end
        fun_writeDeSiOFSIInput(simu_fsi,simu_struct,simu_aero);
        
        % create batch or bash files for os windows/linux
        if strncmp(os,'windows',length(os))
            fun_create_batchrunfile(path_DeSiO, path_mtee, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic, flag_init_self_weight, flag_init_rotor_velocity);
        else
            fun_create_bashrunfile(path_DeSiO, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic, flag_init_self_weight, flag_init_rotor_velocity);
        end
        
    end
    cd(currDir);
    close all;
return
% =================================================================================================================
function fun_create_bashrunfile(path_DeSiO, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic, flag_init_self_weight, flag_init_rotor_velocity)
% write batch file for fsi calculation
    cd([simu_fsi.currDir '\' simu_fsi.strfilename]);
    if not(exist([simu_fsi.currDir '\my_slurm_set.txt'], 'file'))
        % create file with default settings for Luis Cluster at LUH
        fid = fopen(['my_slurm_set.txt'],'w');
            fprintf(fid,'%s\n',['#!/bin/bash -l']);
            fprintf(fid,'%s\n',['#SBATCH --ntasks=10']);
            fprintf(fid,'%s\n',['#SBATCH --nodes=1']);
            fprintf(fid,'%s\n',['#SBATCH --ntasks-per-node=1']);
            fprintf(fid,'%s\n',['#SBATCH --mem=20G']);
            fprintf(fid,'%s\n',['#SBATCH --time=12:00:00']);
            fprintf(fid,'%s\n',['#SBATCH --partition=amo']);
            fprintf(fid,'%s\n',['#SBATCH --output output.out']);
            fprintf(fid,'%s\n',['#SBATCH --error error.err']);
        fclose(fid);
    else
        copyfile([simu_fsi.currDir '\my_slurm_set.txt'],[simu_fsi.currDir '\' simu_fsi.strfilename '\my_slurm_set.txt'])
    end
    fid = fopen('DeSiO.sh','w');
        % set path variables
        fprintf(fid,'%s\n',['#!/bin/bash -l']);
        fprintf(fid,'%s\n',['# set path directories']);
        % set directories
        fprintf(fid,'%s\n',['mypath="$(pwd)"']);
        fprintf(fid,'%s\n',['path_desio="' path_DeSiO '"']);
        fprintf(fid,'%s\n',['path_fsi="$mypath/DeSiO"' ' # to fsi simulation']);
        fprintf(fid,'%s\n',['path_struc="$mypath/DeSiO-Structure"' '  # to structural simulation']);
        % set jobnames
        fprintf(fid,'%s\n',['#']);
        fprintf(fid,'%s\n',['# set jobnames']);
        fprintf(fid,'%s\n',['jobname_fsi="' simu_fsi.strfilename '"']);
        if flag_init_self_weight == 1
            fprintf(fid,'%s\n',['jobname_structure_self_weight="' simu_struct_static.jobname '"']);
        end
        if flag_init_rotor_velocity == 1
            fprintf(fid,'%s\n',['jobname_structure_init_vel="' simu_struct_dynamic.jobname '"']);
        end
        % set job settings in fsi, e.g. for parallelization. create runDeSiO.sh
        fprintf(fid,'%s\n',['#']);
        fprintf(fid,'%s\n',['# copy my_slurm settings to files']);
        fprintf(fid,'%s\n',['cp "$mypath/my_slurm_set.txt" "$path_fsi/runDeSiO.sh"']);
        fprintf(fid,'%s\n',['echo "#SBATCH --job-name=$jobname_fsi" >> "$path_fsi/runDeSiO.sh"']);
        fprintf(fid,'%s\n',['echo "module load intel/2021a" >> "$path_fsi/runDeSiO.sh"']);
        fprintf(fid,'%s\n',['echo "$path_desio" >> "$path_fsi/runDeSiO.sh"']);
        %fprintf(fid,'%s\n',['cat "$mypath/version_settings.txt" >> "$path_fsi/runDeSiO.sh"']);
        fprintf(fid,'%s\n',['#']);
        if flag_init_self_weight == 1
            % set job settings in structure self-weight, e.g. for parallelization. create runDeSiO.sh
            fprintf(fid,'%s\n',['cp "$mypath/my_slurm_set.txt" "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['echo "#SBATCH --job-name=$jobname_structure_self_weight" >> "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['echo "module load intel/2021a" >> "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['echo "$path_desio" >> "$path_struc/$jobname_structure_self_weight/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['#']);
        end
        if flag_init_rotor_velocity == 1
            % set job settings in structure initial velocity, e.g. for parallelization. create runDeSiO.sh
            fprintf(fid,'%s\n',['cp "$mypath/my_slurm_set.txt" "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['echo "#SBATCH --job-name=$jobname_structure_init_vel" >> "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['echo "module load intel/2021a" >> "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['echo "$path_desio" >> "$path_struc/$jobname_structure_init_vel/runDeSiO.sh"']);
            fprintf(fid,'%s\n',['#']);
        end
        if flag_init_self_weight == 1
            % changing dir for static self-weight calculation
            fprintf(fid,'%s\n',['# switching to static self-weight calculations']);
            fprintf(fid,'%s\n',['cd "$path_struc/$jobname_structure_self_weight"']);
            fprintf(fid,'%s\n',['rm "check.log"']);
            fprintf(fid,'%s\n',['#sh runDeSiO.sh']);
            fprintf(fid,'%s\n',['sbatch runDeSiO.sh']);
            fprintf(fid,'%s\n',['until [ -f "check.log" ]']);
            fprintf(fid,'%s\n',['do']);
             fprintf(fid,'%s\n',['   sleep 5']);
            fprintf(fid,'%s\n',['done']);
            if flag_init_rotor_velocity == 1
                fprintf(fid,'%s\n',['# copying result files to initial rotor velocity']);
                fprintf(fid,'%s\n',['cp "${jobname_structure_self_weight}_q.dres" "$path_struc/$jobname_structure_init_vel/${jobname_structure_self_weight}_qs.dres"']);
                fprintf(fid,'%s\n',['cp "${jobname_structure_self_weight}_v.dres" "$path_struc/$jobname_structure_init_vel/${jobname_structure_self_weight}_vs.dres"']);
                fprintf(fid,'%s\n',['cp "${jobname_structure_self_weight}_lambda.dres" "$path_struc/$jobname_structure_init_vel/${jobname_structure_self_weight}_lambdas.dres"']);
                fprintf(fid,'%s\n',['#']);
            else
                fprintf(fid,'%s\n',['# copying result files to initial rotor velocity']);
                fprintf(fid,'%s\n',['cp "${jobname_structure_self_weight}_q.dres" "$path_fsi/${jobname_structure_self_weight}_qs.dres"']);
                fprintf(fid,'%s\n',['cp "${jobname_structure_self_weight}_v.dres" "$path_fsi/${jobname_structure_self_weight}_vs.dres"']);
                fprintf(fid,'%s\n',['cp "${jobname_structure_self_weight}_lambda.dres" "$path_fsi/${jobname_structure_self_weight}_lambdas.dres"']);
                fprintf(fid,'%s\n',['#']);
            end
        end
        if flag_init_rotor_velocity == 1
            % changing dir for dynamic initial velocity calculation
            fprintf(fid,'%s\n',['# switching to dynamic initial rotor velocity calculations']);
            fprintf(fid,'%s\n',['cd "$path_struc/$jobname_structure_init_vel"']);
            fprintf(fid,'%s\n',['rm "check.log"']);
            fprintf(fid,'%s\n',['#sh runDeSiO.sh']);
            fprintf(fid,'%s\n',['sbatch runDeSiO.sh']);
            fprintf(fid,'%s\n',['until [ -f "check.log" ]']);
            fprintf(fid,'%s\n',['do']);
            fprintf(fid,'%s\n',['sleep 5']);
            fprintf(fid,'%s\n',['done']);
            fprintf(fid,'%s\n',['# copying result files to initial rotor velocity']);
            fprintf(fid,'%s\n',['cp "${jobname_structure_init_vel}_q.dres" "$path_fsi/${jobname_structure_init_vel}_qs.dres"']);
            fprintf(fid,'%s\n',['cp "${jobname_structure_init_vel}_v.dres" "$path_fsi/${jobname_structure_init_vel}_vs.dres"']);
            fprintf(fid,'%s\n',['cp "${jobname_structure_init_vel}_lambda.dres" "$path_fsi/${jobname_structure_init_vel}_lambdas.dres"']);
            fprintf(fid,'%s\n',['#']);
        end
        % changing dir for fsi calculation
        fprintf(fid,'%s\n',['#switching to fsi calculation']);
        fprintf(fid,'%s\n',['cd "$path_fsi"']);
        fprintf(fid,'%s\n',['sbatch runDeSiO.sh']);
        fprintf(fid,'%s\n',['#sh runDeSiO.sh']);
    fclose(fid);
% =================================================================================================================
function fun_create_batchrunfile(path_DeSiO, path_mtee, path_fsi, path_struc, simu_fsi, simu_struct_static, simu_struct_dynamic,flag_init_self_weight,flag_init_rotor_velocity)
% write batch file for fsi calculation
    cd(path_fsi)
    fid = fopen('DeSiO.bat','w');
    % set path variables
    fprintf(fid,'%s\n',['set path_DeSiO="' path_DeSiO '"']);
    fprintf(fid,'%s\n',['set path_mtee="' path_mtee '"']);
    fprintf(fid,'%s\n',['REM']);
    fprintf(fid,'%s\n',['set path_fsi="' path_fsi '"']);
    fprintf(fid,'%s\n',['set path_struc="' path_struc '"']);
    % if flag_init_self_weight == 1 run self weight simulation
    if flag_init_self_weight == 1
        fprintf(fid,'%s\n',['REM Start static self-weight analysis with DeSiO']);
        fprintf(fid,'%s\n',['cd %path_struc%' '\' simu_struct_static.jobname]);
        fprintf(fid,'%s\n',['"' path_DeSiO '\DeSiO.exe' '"' ' | ' '"' path_mtee '\mtee.exe' '"' ' "' '%path_struc%' '\' simu_struct_static.jobname '\output.txt"']);
        if flag_init_rotor_velocity == 1 % then copy self weight dres files to init_rotor_vel directory
        fprintf(fid,'%s\n',['REM Copying results files of self weight to dynamic init_rotor_vel simulation']);
            fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_static.jobname '\' simu_struct_static.jobname '_q.dres' ' ' ... 
                        '%path_struc%' '\' simu_struct_dynamic.jobname '\' simu_struct_static.jobname '_qs.dres']);
            fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_static.jobname '\' simu_struct_static.jobname '_v.dres' ' ' ... 
                        '%path_struc%' '\' simu_struct_dynamic.jobname '\' simu_struct_static.jobname '_vs.dres']);
            fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_static.jobname '\' simu_struct_static.jobname '_lambda.dres' ' ' ... 
                        '%path_struc%' '\' simu_struct_dynamic.jobname '\' simu_struct_static.jobname '_lambdas.dres']);
        else % then copy self weight dres files to fsi directory
            fprintf(fid,'%s\n',['REM Copying results files of self weight to fsi simulation']);
                fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_static.jobname '\' simu_struct_static.jobname '_q.dres' ' ' ... 
                        '%path_fsi%' '\' simu_struct_static.jobname '_qs.dres']);
                fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_static.jobname '\' simu_struct_static.jobname '_v.dres' ' ' ... 
                        '%path_fsi%' '\' simu_struct_static.jobname '_vs.dres']);
                fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_static.jobname '\' simu_struct_static.jobname '_lambda.dres' ' ' ... 
                        '%path_fsi%' '\' simu_struct_static.jobname '_lambdas.dres']);
        end
    end
    % if flag_init_rotor_velocity == 1 run init_rotor_vel simulation
    if flag_init_rotor_velocity == 1
        fprintf(fid,'%s\n',['REM Start dynamic init_rotor_vel analysis with DeSiO']);
        fprintf(fid,'%s\n',['cd %path_struc%' '\' simu_struct_dynamic.jobname]);
        fprintf(fid,'%s\n',['"' path_DeSiO '\DeSiO.exe' '"' ' | ' '"' path_mtee '\mtee.exe' '"' ' "' '%path_struc%' '\' simu_struct_dynamic.jobname '\output.txt"']);
        fprintf(fid,'%s\n',['REM Copying results files of dynamic struc simulation for init_rotor_vel to fsi simulation']);
            fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_dynamic.jobname '\' simu_struct_dynamic.jobname '_q.dres' ' ' ... 
                        '%path_fsi%' '\' simu_struct_dynamic.jobname '_qs.dres']);
            fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_dynamic.jobname '\' simu_struct_dynamic.jobname '_v.dres' ' ' ... 
                        '%path_fsi%' '\' simu_struct_dynamic.jobname '_vs.dres']);
            fprintf(fid,'%s\n',['copy ' '%path_struc%' '\' simu_struct_dynamic.jobname '\' simu_struct_dynamic.jobname '_lambda.dres' ' ' ... 
                        '%path_fsi%' '\' simu_struct_dynamic.jobname '_lambdas.dres']);
    end
    fprintf(fid,'%s\n',['REM Start fsi simulation with DeSiO']);
    fprintf(fid,'%s\n',['cd %path_fsi%']);
    fprintf(fid,'%s\n',['REM Boost thread priority']);
    fprintf(fid,'%s\n',['SET desio_exe=' path_DeSiO '\DeSiO.exe']);
    fprintf(fid,'%s\n',['start "" /REALTIME /B /W  %desio_exe%'  ' | ' '"' path_mtee '\mtee.exe' '"' ' "%path_fsi%\output.txt"']);
    %fprintf(fid,'%s\n',['%desio_exe%'  ' | ' '"' path_mtee '\mtee.exe' '"' ' "%path_fsi%\output.txt"']);
    fprintf(fid,'%s\n',['pause']);
    fclose(fid);
% =================================================================================================================
function mesh = fun_set_struc_mesh(mesh,struct_var)
    nmesh = size(mesh,2);
    imat  = 0;
    if nmesh ~= 0
        imat = mesh(nmesh).imat;
    end
    mesh(nmesh+1).nodes               = [struct_var.grid.struc.X_RE,struct_var.grid.struc.D1,struct_var.grid.struc.D2,struct_var.grid.struc.D3];
    mesh(nmesh+1).mx                  = struct_var.grid.struc.M; 
    mesh(nmesh+1).nn                  = (struct_var.grid.struc.M+1);
    mesh(nmesh+1).connectivity(:,1:2) = struct_var.grid.struc.connectivity(:,1:2);
    % element connectivity and material
    mesh(nmesh+1).connectivity(:,3)   = struct_var.grid.struc.matBeam.inz + imat;
    mesh(nmesh+1).imat                = max(mesh(nmesh+1).connectivity(:,3));
return
% =================================================================================================================
function mesh = fun_set_aero_mesh(mesh,struct_var)
    nmesh = size(mesh,2);
    for i = 1:size(struct_var,2)
        mesh(nmesh+i).nodes        = struct_var(i).grid.aero.X; 
        mesh(nmesh+i).connectivity = struct_var(i).grid.aero.connectivity; 
        mesh(nmesh+i).mx           = struct_var(i).grid.aero.N; 
        mesh(nmesh+i).my           = struct_var(i).grid.aero.M;
        mesh(nmesh+i).nn           = (mesh(nmesh+i).mx+1)*(mesh(nmesh+i).my+1);
    end
return
% =================================================================================================================
function fun_plot_3Dmesh(fig,str_b,mesh)
    % plot system
    figure(fig);
    for i = 1:size(mesh,2)
        for i_air = 1:size(mesh(i).connectivity,1)
            if str_b == '2d' | str_b == '2D'
                surf = fill3(mesh(i).nodes(mesh(i).connectivity(i_air,:),1), mesh(i).nodes(mesh(i).connectivity(i_air,:),2), mesh(i).nodes(mesh(i).connectivity(i_air,:),3),'g','facealpha',0.3,'edgealpha',0.5);
            end
        end
%         text(mesh(i).nodes(:,1),mesh(i).nodes(:,2),mesh(i).nodes(:,3),num2str([1:size(mesh(i).nodes,1)]'));
    end
    
return
% =================================================================================================================
function blade_ro_cos = fun_blades2rotor(blade_aero, blade_struc, n_pitch, a_pitch, n_rotor, a_rotor, n_tilt, a_tilt, x_r)
% translating, pitching and rotating of blade to bring in rotor position
    % rotation matrix around pitch axis
    Rp = cos(a_pitch*pi/180)*eye(3) + sin(a_pitch*pi/180)*skew(n_pitch)+(1-cos(a_pitch*pi/180))*n_pitch'*n_pitch;
    % rotation matrix around rotor axis
    Rr = cos(a_rotor*pi/180)*eye(3) + sin(a_rotor*pi/180)*skew(n_rotor)+(1-cos(a_rotor*pi/180))*n_rotor'*n_rotor;
    % rotation matrix around tilt
    Rt = cos(a_tilt*pi/180)*eye(3) + sin(a_tilt*pi/180)*skew(n_tilt)+(1-cos(a_tilt*pi/180))*n_tilt'*n_tilt;
    
    % calculate new position vector for aero grid in global cos
    nnodes = (blade_aero.M+1)*(blade_aero.N+1);
    X_R = [ones(nnodes,1)*x_r(1),ones(nnodes,1)*x_r(2),ones(nnodes,1)*x_r(3)];
    blade_ro_cos.aero.X_C  = (Rt*Rr*(Rp*blade_aero.X_C' + X_R'))';
    
    nnodes = (blade_aero.M+1)*(2*blade_aero.N+1);
    X_R = [ones(nnodes,1)*x_r(1),ones(nnodes,1)*x_r(2),ones(nnodes,1)*x_r(3)];
    blade_ro_cos.aero.X_W  = (Rt*Rr*(Rp*blade_aero.X_W' + X_R'))';
    
    nnodes = (blade_struc.M+1);
    X_R = [ones(nnodes,1)*x_r(1),ones(nnodes,1)*x_r(2),ones(nnodes,1)*x_r(3)];
    blade_ro_cos.struc.X_RE = (Rt*Rr*(Rp*blade_struc.arr_coordinates(:,1:3)' + X_R'))';
    blade_ro_cos.struc.D1 = (Rt*Rr*Rp*blade_struc.arr_coordinates(:, 4:6)')';
    blade_ro_cos.struc.D2 = (Rt*Rr*Rp*blade_struc.arr_coordinates(:, 7:9)')';
    blade_ro_cos.struc.D3 = (Rt*Rr*Rp*blade_struc.arr_coordinates(:,10:12)')';
    return
% =================================================================================================================
function wake = fun_blade_wake(wake, wakesurf, obj, tnrows, nrows, wake_diffusion)
% define seperation edge
    if isempty(wake)
        w = 1;
    else
        w = size(wake.wakes,2)+1;
    end
    
    % wake at trailing edge
    a = 0;
    for j = (obj.mx+1):(obj.mx+1):(obj.mx+1)*(obj.my+1)-1
        a = a + 1;
        [row1,col1,v] = find(obj.connectivity == [j]);
        [row2,col2,v] = find(obj.connectivity(row1,:) == [j+(obj.mx+1)]);
        if isempty(row2); row2 = 1; end
        wake.wakes(w).inf(a,1:3) = [j,j+(obj.mx+1), row1(row2)] ;
    end
    wake.wakes(w).nsurf     = wakesurf;
    wake.wakes(w).nsegments = a;
    wake.wakes(w).nproperty = w;

    % define wake properties
    wake.property(w,1:3) = [tnrows,nrows,wake_diffusion];

    % wake at blade tip
    a = 0;
    for j = (obj.mx+1)*(obj.my+1):-1:(obj.mx+1)*obj.my+1+1
        a = a + 1;
        [row1,col1,v] = find(obj.connectivity == [j]);
        [row2,col2,v] = find(obj.connectivity(row1,:) == [j-1]);
        if isempty(row2); row2 = 1; end
        wake.wakes(w+1).inf(a,1:3) = [j,j-1, row1(row2)] ;
    end
    wake.wakes(w+1).nsurf     = wakesurf;
    wake.wakes(w+1).nsegments = a;
    wake.wakes(w+1).nproperty = w+1;

    % define wake properties
    wake.property(w+1,1:3) = [tnrows,nrows,wake_diffusion];
return