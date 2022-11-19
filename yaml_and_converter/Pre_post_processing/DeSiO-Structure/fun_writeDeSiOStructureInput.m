% =================================================================================================================
% Function for writing DeSiO-Aero input files
% 
% Author: Christian Hente
% Date: 05.05.2022
% 
% input:
%   simu - struc variable that contains DeSiO model components, i.e. rb,
%   pointmass12,constraints,simulation settings
%   mesh - struc variable that containts grid coordinates and
%   connectivities of structural finite beam mesh
% =================================================================================================================
function fun_writeDeSiOStructureInput(simu,mesh)
% =================================================================================================================
currDir = cd;

% creating structure input files and directories
if isfield(simu,'jobname')
    caseDir = [simu.currDir '\' simu.strfilename '\DeSiO-Structure\' simu.jobname];
else
    caseDir = [simu.currDir '\' simu.strfilename '\DeSiO-Structure\'];
end
mkdir(caseDir);

if exist([currDir '\' 'DeSiO.bat'], 'file')
  copyfile([currDir '\' 'DeSiO.bat'],[caseDir]);
else
  warningMessage = sprintf('Warning: file does not exist:\n%s', [currDir '\' 'DeSiO-FSI.bat']);
end

% writing input files
cd(caseDir);
fid = fopen('beaminput.txt','w');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!number of beams (1), number of cross-section properties (2), flag for property input type (3):\n'); 
    fprintf(fid,'%i\t%i\t%i\n',[size(mesh,2),size(simu.matbeam,2), simu.flag_matbeam]);
    arr_node = []; nnodes = 0;
    for i_s = 1:size(mesh,2)
        strbeamname = '';
        if isfield(mesh(i_s),'strname')
            strbeamname = mesh(i_s).strname;
        end        
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! beam %s: number of nodes (1), number of elements (2):\n',strbeamname);
        fprintf(fid,'%i\t%i\n',[mesh(i_s).nn mesh(i_s).mx]);
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! beam %s: nodal phi (1, 2, 3), nodal d1 (4, 5, 6), nodal d2 (7, 8, 9), nodal d3 (10, 11, 12):\n',strbeamname);
        for i = 1:size(mesh(i_s).nodes,1)
            fprintf(fid,'%20.15fd0\t', [mesh(i_s).nodes(i,1:12)]);
            fprintf(fid,'\n');
        end
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! beam %s: connectivities (1, 2), cross-section property (3):\n',strbeamname);
        for i = 1:size(mesh(i_s).connectivity,1)
            fprintf(fid,'%i\t%i\t%i\n',[mesh(i_s).connectivity(i,:)]);        
        end
        % global nodes
        nnodes = nnodes + mesh(i_s).nn;
        arr_node(end+1:nnodes) = [length(arr_node)+1:nnodes];
    end
    if size(simu.matbeam) ~=0 
        for i_mat = 1:size(simu.matbeam,2)
            strbeamname = '';
            if isfield(simu.matbeam(i_mat),'strname')
                strbeamname = simu.matbeam(i_mat).strname;
            end        
            if simu.flag_matbeam == 1 % general DeSiO input-format (Voigt notation)
                fprintf(fid,'!! cross-section property %i (%s) (Voigt notation):\n',i_mat,strbeamname); 
                fprintf(fid,'!! row1: cbeam\n');
                fprintf(fid,'!! row2: cmass\n');
            else % isotropic material DeSiO input-format 
                fprintf(fid,'!! cross-section property %i (%s):\n',i_mat,strbeamname);
                fprintf(fid,'!! row1: EA(1), GA1(2), GA2(3), EI1(4), EI2(5), GI3(6), ES1(7), ES2(8), GS1(9), GS2(10), EI12(11)\n');
                fprintf(fid,'!! row2: rhoA3(1), rhoI1(2), rhoI2(3), rhoS1(4), rhoS2(5), rhoI12(6)\n');
            end
            str = strrep(sprintf('%10.5e\t',simu.matbeam(i_mat).cmat),'e','d');
            fprintf(fid,'%s',str); fprintf(fid,'\n');
            
            str = strrep(sprintf('%10.5e\t',simu.matbeam(i_mat).mmat),'e','d');
            fprintf(fid,'%s',str); fprintf(fid,'\n');
            
            str = strrep(sprintf('%10.5e\t',simu.matbeam(i_mat).diss),'e','d');
            fprintf(fid,'%s',str); fprintf(fid,'\n');
        end
    end
fclose(fid);

if isfield(simu,'pointmass12')
    fid = fopen('pointmass12input.txt','w');
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! number of point masses\n');
        fprintf(fid,'%i\n',size(simu.pointmass12,2));
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! (1) node, (2) point mass:\n');
        for ipm12 = 1:size(simu.pointmass12,2)
            fprintf(fid,'%i\t%10.5fd0\n',simu.pointmass12(ipm12).node,simu.pointmass12(ipm12).mass);
        end
    fclose(fid);
end

if isfield(simu,'constraints')
    fid = fopen('constraint12input.txt','w');
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! number of constraints for nodes with 12 coordinates:\n');
        fprintf(fid,'%i\n',size(simu.constraints,2));
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! constraints for nodes with 12 coordinates: sort (1), nodes (2, 3), phi1 (4, 5, 6), phi2 (7, 8, 9), dir (10, 11, 12):\n');
        for ic = 1:size(simu.constraints,2)
            phi1 = [0,0,0]; phi2 = [0,0,0]; dir  = [0,0,0];
            fprintf(fid,'%s\t',simu.constraints(ic).type);
            fprintf(fid,'%i\t',simu.constraints(ic).nodes);
            if isfield(simu.constraints(ic),'phi1')
                if not(isempty(simu.constraints(ic).phi1)); phi1 = simu.constraints(ic).phi1; end
            end
            if isfield(simu.constraints(ic),'phi2')
                if not(isempty(simu.constraints(ic).phi2)); phi2 = simu.constraints(ic).phi2; end
            end
            if isfield(simu.constraints(ic),'dir')
                if not(isempty(simu.constraints(ic).dir));  dir  = simu.constraints(ic).dir; end
            end
            
            str = strrep(sprintf('%10.8e\t',phi1),'e','d');
            fprintf(fid,'%s',str);
            
            str = strrep(sprintf('%10.8e\t',phi2),'e','d');
            fprintf(fid,'%s',str);
            
            str = strrep(sprintf('%10.8e\t',dir),'e','d');
            fprintf(fid,'%s',str);
            
            fprintf(fid,'\n');
        end
    fclose(fid);
end

if isfield(simu,'rb')
    fid = fopen('rigidbodyinput.txt','w');
        fprintf(fid,'!! rigid body input\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! number of rigid bodies (1), number of body properties (2):\n');
        fprintf(fid,'%i\t%i\n',[size(simu.rb,2), size(simu.rb,2)]);
        for irb = 1:size(simu.rb,2)
            type = '';
            fprintf(fid,'!!\n');
            fprintf(fid,'!!\n');
            if isfield(simu.rb(irb),'type'); type = simu.rb(irb).type; end;
            fprintf(fid,'!! rigid body %i, %s phi(1:3), d1(4:6) d2(7:9), d3(10:12)\n',irb,type);
            fprintf(fid,'%10.8fd0\t',simu.rb(irb).center_mass);
            fprintf(fid,'%10.8fd0\t',simu.rb(irb).D1);
            fprintf(fid,'%10.8fd0\t',simu.rb(irb).D2);
            fprintf(fid,'%10.8fd0\t',simu.rb(irb).D3);
            fprintf(fid,'\n');
            fprintf(fid,'!!\n');
            fprintf(fid,'!!\n');
            fprintf(fid,'!! rigid body %i, %s property\n',irb,type);
            fprintf(fid,'%i\n',irb);
        end
        for irb = 1:size(simu.rb,2)
            type = '';
            fprintf(fid,'!!\n');
            fprintf(fid,'!!\n');
            if isfield(simu.rb(irb),'type'); type = simu.rb(irb).type; end;
            fprintf(fid,'!! rigid body %i property, %s hub mass(1), J11(2), J22(3), J33(4), J23(5), J13(6), rhoxhi3(7), rhoxhi2(8), rhoxhi1(9), J12(10)\n',irb,type);
            fprintf(fid,'%10.8fd0\t',simu.rb(irb).mass_matrix);
            fprintf(fid,'\n');
        end
    fclose(fid);    
end

strtypecomments = {'row 1: simulation settings: totalt (1), deltat (2), tolerance (3), iteration limit (4), gravity flag (5)',...   % dynamic and static solver
                   'row 1: simulation settings: number of EV (1), tolerance (2), emin, emax, flag (0 - dense; 1 - sparse)',...      % modal and linear buckling solver
                   'row 1: simulation settings: totalt (1), deltat (2), tolerance (3), iteration limit (4), arc length method (5), desired iterations (6)'}; % static post-buckling solver

if isfield(simu,'type')
    strjobname     = 'solution'; 
    strprevjobname = 'none';
    if strncmp(simu.type,'dynamic',7)
        if isfield(simu,'jobname'); strjobname = simu.jobname; end
        if isfield(simu,'prevjobname'); strprevjobname = simu.prevjobname; end
        strtype                  = 'dynamic';
        strtypecomment           = strtypecomments{1};
        str_simu_settings_format = ['%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\t%i\n'];
        simu_settings            = [1.0e0 1.0e-1 1.0e-6 50 0];
        if isfield(simu,'settings'); simu_settings = simu.settings; end
        simu_grav_vec = [0,0,0]; 
        if isfield(simu,'grav_vec'); simu_grav_vec = simu.grav_vec; end
    elseif strncmp(simu.type,'static',6)
        if isfield(simu,'jobname'); strjobname = simu.jobname; end
        if isfield(simu,'prevjobname'); strprevjobname = simu.prevjobname; end
        strtype                  = 'static';
        strtypecomment           = strtypecomments{1};
        str_simu_settings_format = ['%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\t%i\n'];
        simu_settings            = [1.0e0 1.0e-1 1.0e-6 50 0];
        if isfield(simu,'settings'); simu_settings = simu.settings; end
        simu_grav_vec = [0,0,0]; 
        if isfield(simu,'grav_vec'); simu_grav_vec = simu.grav_vec; end    
    elseif strncmp(simu.type,'modal',5)
        if isfield(simu,'jobname'); strjobname = simu.jobname; end
        if isfield(simu,'prevjobname'); strprevjobname = simu.prevjobname; end
        strtype                  = 'modal';
        strtypecomment           = strtypecomments{2};
        str_simu_settings_format = ['%i\t%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\n'];
        simu_settings            = [10 1.0e-6 0.0 1.0e0 0];
        if isfield(simu,'settings'); simu_settings = simu.settings; end
        simu_grav_vec = [0,0,0]; 
    elseif strncmp(simu.type,'buckling',8)
        if isfield(simu,'jobname'); strjobname = simu.jobname; end
        if isfield(simu,'prevjobname'); strprevjobname = simu.prevjobname; end
        strtype                  = 'buckling';
        strtypecomment           = strtypecomments{2};
        str_simu_settings_format = ['%i\t%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\n'];
        simu_settings            = [10 1.0e-6 0.0 1.0e0 0];
        if isfield(simu,'settings'); simu_settings = simu.settings; end
        simu_grav_vec = [0,0,0]; 
        if isfield(simu,'grav_vec'); simu_grav_vec = simu.grav_vec; end
    elseif strncmp(simu.type,'static_arc',10)
        if isfield(simu,'jobname'); strjobname = simu.jobname; end
        if isfield(simu,'prevjobname'); strprevjobname = simu.prevjobname; end
        strtype                  = 'static_arc';
        strtypecomment           = strtypecomments{3};
        str_simu_settings_format = ['%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\t%i\n'];
        simu_settings            = [1.0d0 1.0d0 1.0d-6 50 2 4];
        if isfield(simu,'settings'); simu_settings = simu.settings; end
        simu_grav_vec = [0,0,0]; 
        if isfield(simu,'grav_vec'); simu_grav_vec = simu.grav_vec; end
    end
else
        strjobname               = 'solution';
        strprevjobname          = 'none';
        strtype                  = 'modal';
        strtypecomment           = strtypecomments{2};
        str_simu_settings_format = ['%i\t%10.8fd0\t%10.8fd0\t%10.8fd0\t%i\n'];
        simu_settings            = [10 1.0e-6 0.0 1.0e0 0];
        simu_grav_vec            = [0,0,0]; 
end

  fid = fopen('simulationinput_structure.txt','w');
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! filename\n');
        fprintf(fid,'%s\t%s\n',strjobname,strprevjobname);
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! number of simulations:\n');
        fprintf(fid,'%i\n',1);
        fprintf(fid,'!! ... place for comments ...\n');
        fprintf(fid,'!! %s\n',strtypecomment);
        fprintf(fid,'!! row2: flag for writing matrices - on/off = 1/0\n');
        fprintf(fid,'%s\n',strtype);
        fprintf(fid,str_simu_settings_format,simu_settings);
        fprintf(fid,'0\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! gravity vector (1, 2, 3):\n');
        fprintf(fid,'%10.5fd0',simu_grav_vec);
        fprintf(fid,'\n');
    fclose(fid);

if isfield(simu,'loads')
    fid = fopen('load12input.txt','w');
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! number of loads for nodes with 12 coordinates.\n');
        fprintf(fid,'%i\n',size(simu.loads,2));
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! loads for nodes with 12 coordinates: sort (1), intensity (2), duration (3), node (4), spatial (5, 6, 7, 8, 9, 10), material (11, 12, 13, 14, 15, 16).\n');
        for i = 1:size(simu.loads,2)
            fprintf(fid,'%s',simu.loads(i).sort)
            fprintf(fid,'%20.15fd0\t%20.15fd0\t%i\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\n',...
                    [simu.loads(i).intensity,simu.loads(1).duration,simu.loads(1).node,simu.loads(1).spatial,simu.loads(1).material]);
        end
    fclose(fid);
end

if isfield(simu,'boundary12')
    fid = fopen('boundary12input.txt','w');
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! number of inhomogeneous constraints for nodes with 12 coordinates.\n');
        fprintf(fid,'%i\n',size(simu.boundary12,1));
        fprintf(fid,'!!\n');
        fprintf(fid,'!! inhomogeneous constraints for nodes with 12 coordinates: \n');
        fprintf(fid,'!! sort (1), intensity (2), duration (3), constraint ID (4), file (5) (if no, then none).\n');
        for i = 1:size(simu.boundary12,2)
            fprintf(fid,'%s\t%20.15fd0\t%20.15fd0\t%i\t%s\n',...
                    simu.boundary12(i).sort,...
                     simu.boundary12(i).intensity,...
                     simu.boundary12(i).time,...
                     simu.boundary12(i).constraintID,...
                     simu.boundary12(i).file);
        end
    fclose(fid);
end

% =================================================================================================================
% close all;
cd(simu.currDir);
disp('creating DeSiO-Structure input files');
% =================================================================================================================
return