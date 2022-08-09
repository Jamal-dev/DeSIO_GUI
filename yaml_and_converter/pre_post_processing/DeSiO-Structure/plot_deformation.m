% =========================================================================
% Plotting position vector
% =========================================================================
clc;
clear all;
close all

addpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\03_Prepostprocessing\DeSiO-Structure');

model = strucure_readmodel;

% loading of DeSiO result files
q = load([model.strSimName '_q.dres']);
time = load([model.strSimName '_t.dres']);

% function to extract dof-solution for node from solution.m file
node  = [5];
for i = 1:length(node)
    [u] = get_DeSiO_dof_solu(node(i),1,q);
    figure(); hold on; grid on;
    plot(time(:,1),u(:,1),'-r');
    plot(time(:,1),u(:,2),'-b');
    plot(time(:,1),u(:,3),'-k');
    ylabel('displacement m'); xlabel('time s')
    title(['displacement vs. time of node' num2str(node(i))]);
    legend('ux','uy','uz');
end
set(gca,'fontsize',12,'fontweight','bold');
% print('displ_time','-dpng', '-r500');
return