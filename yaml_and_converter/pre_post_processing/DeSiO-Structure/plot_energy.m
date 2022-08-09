% =========================================================================
% Plotting DeSiO Energy-Solution
% =========================================================================
clc;
clear all;
close all;

addpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\03_Prepostprocessing\DeSiO-Structure');

model = strucure_readmodel;

% Loading DeSiO result files
q = load([model.strSimName '_q.dres']);
invariants = load([model.strSimName '_e.dres']);

figure(); hold on; grid on;
plot(time(:,1),invariants(:,1),'-b');
plot(time(:,1),invariants(:,2),'-r');
plot(time(:,1),invariants(:,3),'-g');
leg = legend('kinetic energy', 'potential energy', 'total energy');
ylabel('energy Nm'); xlabel('time s');
set(gca,'fontsize',12,'fontweight','bold');
title(['energy vs. time']);
% print('energy_time','-dpng', '-r500');
return