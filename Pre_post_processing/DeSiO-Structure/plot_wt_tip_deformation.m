% =========================================================================
% Plotting DeSiO DOF-Solution
% =========================================================================
clc;
clear all;
close all
format short

addpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\01_Nullversion\Pre_post_processing\DeSiO-Structure');
addpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\01_Nullversion\Pre_post_processing\DeSiO-FSI');
model = fsi_readmodel;

i1  = [1,0,0]';
i2  = [0,1,0]';
i3  = [0,0,1]';

% loading of DeSiO result files
q = load([model.strSimName '_q.dres']);
time = load([model.strSimName '_t.dres']);

% Function to extract dof-solution for node from solution.m file
node    = 53;
noderef = 2;
for i = 1:length(node)
    % relativ Verformung der Blattspitze
    inz    = [12*(node-1)+1:12*(node-1)+12]; 
    inzref = [12*(noderef-1)+1:12*(noderef-1)+12];
    d10    = q(1,inzref(4:6))'/norm(q(1,inzref(4:6)));
    d20    = q(1,inzref(7:9))'/norm(q(1,inzref(7:9)));
    d30    = q(1,inzref(10:12))'/norm(q(1,inzref(10:12)));
    d1n    = d10; d2n = d20; d3n = d30;
    
    xn      = q(1,inz(1:3))';
    xrefn   = q(1,inzref(1:3))';
    inztime = [1:1:size(time,1)];
    for j = 1:length(inztime)
        x    = q(inztime(j),inz(1:3))';
        xref = q(inztime(j),inzref(1:3))';
        d1 = q(inztime(j),inzref(4:6))'/norm(q(j,inzref(4:6)));
        d2 = q(inztime(j),inzref(7:9))'/norm(q(j,inzref(7:9)));
        d3 = q(inztime(j),inzref(10:12))'/norm(q(j,inzref(10:12)));
        R  = d1*d1n' + d2*d2n' + d3*d3n';
        
        % Vector of tip deflection in initial COS
        urel(j,1:3) = x - (R*(xn-xrefn)+xrefn);
        
        % Vector of tip deflection in director cos
        d1no = q(inztime(j),inz(4:6))'/norm(q(j,inz(4:6)));
        d2no = q(inztime(j),inz(7:9))'/norm(q(j,inz(7:9)));
        d3no = q(inztime(j),inz(10:12))'/norm(q(j,inz(10:12)));
        urel_d(j,1) = urel(j,1:3)*d1no;
        urel_d(j,2) = urel(j,1:3)*d2no;
        urel_d(j,3) = urel(j,1:3)*d3no;
        
        d1n = d1; d2n = d2; d3n = d3; xn = x; xrefn = xref;
    end
    figure(); hold on; grid on; 
%     axis([0 600 -0.30 0.30]);
    plot(time(inztime,1),urel_d(:,1),'-r');
    plot(time(inztime,1),urel_d(:,2),'linestyle','-','color',[0.3027    0.8299    0.4687]);
    plot(time(inztime,1),urel_d(:,3),'-b');
    ylabel('tip-deflection m'); xlabel('time s')
    title(['tip-deflection vs. time']);
    legend('edge-wise','flap-wise','longitudinal');
end
saveFigure(gcf,'wt_tip_displ.fig')
return

