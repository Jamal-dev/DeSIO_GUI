% =================================================================================================================
% Function for writing DeSiO-Aero input files
% 
% Author: Christian Hente
% Date: 05.05.2022
% 
% input:
%   simu - struc variable that contains DeSiO-FSI simulation settings, i.e.
%   jobname, search radius
%   simu_struc - settings for structural simulation
%   simu_aero - settings for aerodynamic simulation
% =================================================================================================================
function fun_writeDeSiOFSIInput(simu,simu_struc,simu_aero)
% =================================================================================================================
currDir = cd;

% creating uvlm input files and directories
caseDir = [simu.currDir '\' simu.strfilename '\DeSiO\'];
mkdir(caseDir);

if exist([currDir '\' 'DeSiO.bat'], 'file')
  copyfile([currDir '\' 'DeSiO.bat'],[caseDir]);
else
  warningMessage = sprintf('Warning: file does not exist:\n%s', [currDir '\' 'DeSiO-FSI.bat']);
end

cd(caseDir);
fid = fopen('simulationinput_fsi.txt','w');
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! filename (1), prev filename aero (2), prev filename structure (3)\n');
    fprintf(fid,'%s\n',simu.strfilename);
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! total time (1) deltat aero (2) cutoff (3) time for load factor (4) type of load factor function (5)\n');
    
    lf_duration = simu_aero.time;
    if isfield(simu,'lf_duration'); lf_duration = simu.lf_duration; end
    
    lf_type = 'constant';
    if isfield(simu,'lf_type'); lf_type = simu.lf_type; end
    
    fprintf(fid,'%20.15fd0\t%20.15fd0\t%20.15fd0\t%20.15fd0\t%s\n',simu_aero.time,simu_aero.deltat,simu_aero.cutoff,lf_duration,lf_type);
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! wind data: sort(1), fluid density(2), intensity(3), duration(4), direction vector (5-7)\n');
    fprintf(fid,'%s\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\t%20.10fd0\n','constant',simu_aero.density,simu_aero.i_vinf,simu_aero.time,simu_aero.d_vinf);
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! deltat (1), tolerance (2), iteration limit (3), gravity flag (4), row2: flag for writing matrices - on/off = 1/0\n');
    
    struc_simType = 'dynamic';
    if isfield(simu_struc,'simType'); struc_simType = simu_struc.simType; end
    
    fprintf(fid,'%s\n',struc_simType);
    fprintf(fid,'%20.10fd0\t%20.10fd0\t%i\t%i\n',simu_struc.deltat,simu_struc.tol,simu_struc.niter,simu_struc.grav);
    
    flag_linearization = 0;
    if isfield(simu,'flag_linearization'); flag_linearization = simu.flag_linearization; end
    
    fprintf(fid,'%i\t%i\n',0,flag_linearization);
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! gravity vector (1, 2, 3)\n');
    fprintf(fid,'%20.10fd0\t%20.10fd0\t%20.10fd0\n',simu_struc.grav_vec);
fclose(fid);

fid = fopen('fsi_input.txt','w');
    nfsi = 0;
    if isfield(simu,'data'); nfsi = size(simu.data,1); end;
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! number of fsi (1)\n');
    fprintf(fid,'%i\n',nfsi);
    fprintf(fid,'!!\n');
    fprintf(fid,'!!\n');
    fprintf(fid,'!! input for fluid-structure interaction: from shell/beaminput (1), surface from surfaceinput (2) local search radius (3)\n');
    for i = 1:nfsi
        fprintf(fid,'%s\t%i\t%i\t%20.10fd0\n',simu.data{i,1},simu.data{i,2},simu.data{i,3},simu.data{i,4});
    end
fclose(fid);

% =================================================================================================================
cd(simu.currDir);
disp('creating DeSiO-FSI input files');
% =================================================================================================================
return