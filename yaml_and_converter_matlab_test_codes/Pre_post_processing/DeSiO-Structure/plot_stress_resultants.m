% =========================================================================
% Plotting DeSiO Stress Resultants
% =========================================================================
clc; clear all; close all

addpath(genpath('U:\PhD-2021-12_Christian_Hente\Pre_post_processing'));

model = structure_readmodel;

% Loading of DeSiO result files
time   = load([model.strSimName '_t.dres']);
stress = load([model.strSimName '_stress.dres']);

% Function to extract dof-solution for node from solution.m file
element = [50+1,50+100+1,50+2*100+1];
for i = 1:length(element)
    F1 = stress(:,6*(element(i)-1)+1); % shear force in local 1
    F2 = stress(:,6*(element(i)-1)+2); % shear force in local 2
    F3 = stress(:,6*(element(i)-1)+3); % normal force in local 3
    M1 = stress(:,6*(element(i)-1)+4); % bending moment around local 1
    M2 = stress(:,6*(element(i)-1)+5); % bending moment around local 2
    M3 = stress(:,6*(element(i)-1)+6); % torsion moment around local 3
    
    figure();
    for j = 1:3
        eval(['subplot(2,3,' num2str(j) ');']);
        eval(['ylabel(''F'  num2str(j) '  N '');']);
        xlabel('time s'); hold on; grid on;
        set(gca,'fontsize',8,'fontweight','bold');
        eval(['title([''F' num2str(j) ' vs. time of element '' num2str(element(i))]);']);
        eval(['plot(time(:,1),F' num2str(j) ',''-r'');']);
    end
        for j = 1:3
        eval(['subplot(2,3,' num2str(j+3) ');']);
        eval(['ylabel(''M'  num2str(j) '  Nm '');']);
        xlabel('time s'); hold on; grid on;
        set(gca,'fontsize',8,'fontweight','bold');
        eval(['title([''M' num2str(j) ' vs. time of element '' num2str(element(i))]);']);
        eval(['plot(time(:,1),M' num2str(j) ',''-b'');']);
    end
end
% print('stress_time','-dpng', '-r500');
return