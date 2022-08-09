% =================================================================================================================
% Function for writing DeSiO-Aero input files
% 
% Author: Christian Hente
% Date: 05.05.2022
% 
% input:
%   simu - struc variable that contains simulation specified parameter (jobname,
%   deltat, cutoff, density, i_vinf, d_vinf
%   mesh - struc variable that containts grid coordinates and
%   connectivities
%   wake - struc variable that containts wake information data
% =================================================================================================================
function fun_writeDeSiOAeroInput(simu,mesh,wake)
% =================================================================================================================
currDir = cd;

% creating uvlm input files and directories
caseDir = [simu.currDir '\' simu.strfilename '\DeSiO-Aero\'];
mkdir(caseDir);

% checking, if DeSiO.bat exists
if exist([currDir '\' 'DeSiO.bat'], 'file')
  copyfile([currDir '\' 'DeSiO.bat'],[caseDir]);
else
  warningMessage = sprintf('Warning: file does not exist:\n%s', [currDir '\' 'DeSiO-FSI.bat']);
end

% writing input files
cd(caseDir);
fid = fopen('surfaceinput.txt','w');
    fprintf(fid,'!! \n');
    fprintf(fid,'!! \n');
    fprintf(fid,'!! number of surfaces (1):\n');
    fprintf(fid,'%i \n',size(mesh,2));
    for i_s = 1:size(mesh,2)
        fprintf(fid,'!! \n');
        fprintf(fid,'!! \n');
        fprintf(fid,'!! !! surface: number of nodes (1), number of rings (2), number of nodes along first dimension (3), and number of nodes along second direction (4). \n');
        fprintf(fid,'%i %i %i %i \n',[mesh(i_s).nn mesh(i_s).mx*mesh(i_s).my mesh(i_s).mx+1 mesh(i_s).my+1]);
        fprintf(fid,'!! \n');
        fprintf(fid,'!! \n');
        fprintf(fid,'!! surface node coordinates: nodal phi (1, 2, 3). \n');
        for i = 1:size(mesh(i_s).nodes,1)
            fprintf(fid,'%20.10fd0 %20.10fd0 %20.10fd0 \n', [mesh(i_s).nodes(i,1), mesh(i_s).nodes(i,2), mesh(i_s).nodes(i,3)]);
        end
        fprintf(fid,'!! \n');
        fprintf(fid,'!! \n');
        fprintf(fid,'!! surface: connectivities (1, 2, 3, 4). \n');
        for i = 1:size(mesh(i_s).connectivity,1)
            fprintf(fid,'%i %i %i %i \n',mesh(i_s).connectivity(i,:));        
        end
    end
fclose(fid);

if isempty(wake)
    n_wakes = 0;
    nprop   = 0;
else
    n_wakes = size(wake.wakes,2);
    nprop   = size(wake.property,1);
end
  
fid = fopen('wakeinput.txt','w');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! number of wakes (1), number of wake properties (2)\n');
    fprintf(fid,'%i\t%i\n',[n_wakes, nprop]);
    for i = 1:n_wakes
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! wake %i: number of surface (1), number of segments (2), property number (3)\n', i);
        fprintf(fid,'%i\t%i\t%i\n',[wake.wakes(i).nsurf, wake.wakes(i).nsegments, wake.wakes(i).nproperty]);
        fprintf(fid,'!!\n');
        fprintf(fid,'!!\n');
        fprintf(fid,'!! wake %i: nodes (1, 2), ring (3)\n', i);
        for j = 1:size(wake.wakes(i).inf,1)
            fprintf(fid,'%i\t%i\t%i\n',[wake.wakes(i).inf(j,:)]);
        end
    end
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! wake property %i: time cut (1), max. rows (2)\n', i);
    for i = 1:nprop
        fprintf(fid,'%i\t%i\t0.0d0\t\n',wake.property(i,:));
    end
fclose(fid);

fid = fopen('simulationinput_aero.txt','w');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! filename\n');
    fprintf(fid,'solution\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! simulation settings: totalt (1), deltat (2), cutoff (3)\n');
    fprintf(fid,'%20.15fd0\t%20.15fd0\t%20.15fd0\n',[simu.time,simu.deltat,simu.cutoff]);
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! wind data: sort(1), fluid density(2), intensity(3), duration(4), direction vector (5-7)\n');
    fprintf(fid,'constant\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\n',[simu.density simu.i_vinf simu.time simu.d_vinf]);
fclose(fid);
% =================================================================================================================
cd(simu.currDir);
disp('creating DeSiO-Aero input files');
% =================================================================================================================
return