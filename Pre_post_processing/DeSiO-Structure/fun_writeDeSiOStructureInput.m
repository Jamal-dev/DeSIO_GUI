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
caseDir = [simu.currDir '\' simu.strfilename '\DeSiO-Structure\'];
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
            
%             fprintf(fid,'%10.8fd0\t',phi1);
%             fprintf(fid,'%10.8fd0\t',phi2);
%             fprintf(fid,'%10.8fd0\t',dir);
            
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

fid = fopen('simulationinput_structure.txt','w');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! filename\n');
    fprintf(fid,'solution\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! number of simulations:\n');
    fprintf(fid,'%i\n',1);
    fprintf(fid,'!! Settings for nonlinear dynamic solver with constant time stepping\n');
    fprintf(fid,'!! row1: simulation settings: totalt (1), deltat (2), tolerance (3), iteration limit (4), gravity flag (5)\n');
    fprintf(fid,'!! row2: flag for writing matrices - on/off = 1/0\n');
    fprintf(fid,'modal\n');
    fprintf(fid,'10\t1.0d-8\t0.0d0\t1.0d8\t0\n');
    fprintf(fid,'0\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! gravity vector (1, 2, 3):\n');
    fprintf(fid,'%10.5fd0',[simu.grav_vec]);
    fprintf(fid,'\n');
fclose(fid);

if size(simu.loads,1) ~= 0
    fid = fopen('load12input.txt','w');
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! number of loads for nodes with 12 coordinates.\n');
        fprintf(fid,'%i\n',size(simu.loads,1));
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! loads for nodes with 12 coordinates: sort (1), intensity (2), duration (3), node (4), spatial (5, 6, 7, 8, 9, 10), material (11, 12, 13, 14, 15, 16).\n');
        for i = 1:size(simu.loads,1)
            fprintf(fid,'%s',simu.loads(i).sort)
            fprintf(fid,'%20.15fd0\t%20.15fd0\t%i\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\n',...
                    [simu.loads(i).intensity,simu.loads(1).duration,simu.loads(1).node,simu.loads(1).spatial,simu.loads(1).material]);
        end
    fclose(fid);
end

% =================================================================================================================
% close all;
cd(simu.currDir);
disp('creating DeSiO-Structure input files');
% =================================================================================================================
return