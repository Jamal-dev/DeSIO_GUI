% script for calculating external moments from director-based formulation
clc; close all; 
clear all;

addpath(genpath('C:\Users\hente\Desktop\digital-twin-crc-1463\Z01\DeSiO\Pre_post_processing'));

model   = structure_readmodel;

q = load([model.strSimName '_q.dres']);
t = load([model.strSimName '_t.dres']); t = t(:,1);
E = load([model.strSimName '_e.dres']);

node1 = 2;          % hub
node2 = 2+21+1;     % 1st node of blade 1

% indizes for node
inz1 = 12*(node1-1)+1:12*(node1-1)+12;
inz2 = 12*(node2-1)+1:12*(node2-1)+12;

% drehwinkel and winkelgeschw
h = [];
alfa = 0;
for i = 1:length(t)-1
    n_o  = q(i,inz1(1:3)) - q(i,inz2(1:3));
    d1o  = q(i,inz1(4:6)); 
    d2o  = q(i,inz1(7:9)); 
    d3o  = q(i,inz1(10:12));
    
    % nd_o - n_o axis in director coordinate system
    nd_o = [n_o*d1o';n_o*d2o';n_o*d3o'];
    nd_o = nd_o/norm(nd_o);
    
    n_n  = q(i+1,inz1(1:3)) - q(i+1,inz2(1:3));
    d1n  = q(i+1,inz1(4:6)); 
    d2n  = q(i+1,inz1(7:9)); 
    d3n  = q(i+1,inz1(10:12));
    % nd_n - n_o axis in director coordinate system
    nd_n = [n_n*d1n';n_n*d2n';n_n*d3n']; nd_n = nd_n/norm(nd_n);
    nd_o = [n_o*d1n';n_o*d2n';n_o*d3n']; nd_o = nd_o/norm(nd_o);
        
    dalfa = real(acos(nd_o'*nd_n));
    alfa(i+1) = alfa(i) + dalfa;
end

for i = 1:length(alfa)-1
    omega(i+1) = (alfa(i+1)-alfa(i))/(t(i+1)-t(i));
end

figure(); hold on; grid on;
plot(t,alfa*180/pi);
ylabel('angle [grad]');
xlabel('time s');
set(gca,'fontsize',12,'fontweight','bold');

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

