% script for calculating external moments from director-based formulation
clc; close all; 
clear all;

addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));

model_fsi       = fsi_readmodel;
model_structure = structure_readmodel;
if isempty(model_fsi) && not(isempty(model_structure))
    model = model_structure;
elseif isempty(model_structure) && not(isempty(model_fsi))
    model = model_fsi;
else
    disp(['No model found!']);
    return;
end

q = load([model.strSimName '_q.dres']);
s = load([model.strSimName '_v.dres']);
t = load([model.strSimName '_t.dres']); t = t(:,1);
E = load([model.strSimName '_e.dres']);


% indizes for node
node = 2;          % node at which angular velocity is computed
nd   = [-7.03233176e-01;7.03233176e-01;1.04528463e-01]; % axis in local director cos of node

% drehwinkel and winkelgeschw
inz = 12*(node-1)+1:12*(node-1)+12;
for i = 1:length(t)
    
    d1  = q(i,inz(4:6))'; 
    d2  = q(i,inz(7:9))'; 
    d3  = q(i,inz(10:12))';
    
    w1  = s(i,inz(4:6))'; 
    w2  = s(i,inz(7:9))'; 
    w3  = s(i,inz(10:12))';
    
    n   = nd(1)*d1 + nd(2)*d2 + nd(3)*d3;
    omega(i) = n'*0.5*(cross(d1,w1) + cross(d2,w2) + cross(d3,w3));
end

figure(); hold on; grid on;
plot(t,omega/(2*pi)*60);
ylabel('angular velocity [U/m]');
xlabel('time s');
set(gca,'fontsize',12,'fontweight','bold');

figure(); hold on; grid on;
plot(t(1:size(E,1)),E(:,1),'-b'); 
plot(t(1:size(E,1)),E(:,2),'-r');
plot(t(1:size(E,1)),E(:,3),'-g');
ylabel('energy in Nm')
xlabel('time in s')
legend('kinetic energy', 'potential energy','total energy');
set(gca,'fontsize',12,'fontweight','bold');
return

