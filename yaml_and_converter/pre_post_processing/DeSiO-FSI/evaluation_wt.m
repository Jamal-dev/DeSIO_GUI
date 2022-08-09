% script for calculating external moments from director-based formulation
clc; close all; clear all;

addpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\03_Prepostprocessing\DeSiO-FSI');
addpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\03_Prepostprocessing\DeSiO-Aero');
addpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\03_Prepostprocessing\DeSiO-Structure');

model   = uvlm_readmodel;

q = load([model.strSimName '_q.dres']);
t = load([model.strSimName '_t.dres']); t = t(:,1);
E = load([model.strSimName '_e.dres']);

node = 2;

% indizes for node
inz = 12*(node-1)+1:12*(node-1)+12;

% drehwinkel and winkelgeschw
h = [];
alfa = 0;
for i = 1:length(t)-1
    q_node_o = q(i,inz);
    q_node_n = q(i+1,inz);
    dalfa = real(acos(q_node_o(10:12)*q_node_n(10:12)'));
    alfa(i+1) = alfa(i) + dalfa;
end

for i = 1:length(alfa)-1
    omega(i) = (alfa(i+1)-alfa(i))/(t(i+1)-t(i));
end
omega(i+1) = omega(i);

figure(); hold on; grid on;
plot(t,omega/(2*pi)*60);
ylabel('angular velocity U/m');
xlabel('time s');
set(gca,'fontsize',12,'fontweight','bold');

figure(); hold on; grid on;
plot(t(1:size(E,1)),E(:,1),'-g'); 
plot(t(1:size(E,1)),E(:,2),'-k');
plot(t(1:size(E,1)),E(:,3),'-b');
ylabel('energy J')
xlabel('time s')
legend('kinetic energy', 'potential energy','total energy');
set(gca,'fontsize',12,'fontweight','bold');
return

