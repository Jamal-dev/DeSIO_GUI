% =========================================================================
% Plotting DeSiO DOF-Solution
% =========================================================================
clc;
clear all;
% close all
format short

addpath(genpath('T:\Projekte\Dynamics-2021-06_SFB1463\05_Teilprojekte\Z01\01_Digitaler_Zwilling\01_Nullversion\Pre_post_processing'));

model = structure_readmodel;

i1  = [1,0,0]';
i2  = [0,1,0]';
i3  = [0,0,1]';

% flag for plot output
flag_type = 1; 
% flag_type(1) - relative tip displacement
% flag_type(2) - global tip displacement

% loading of DeSiO result files
q    = load([model.strSimName '_q.dres']);
time = load([model.strSimName '_t.dres']);

% Function to extract dof-solution for node from solution.m file
arrnode = [2+51+101,2+51+201,2+51+301];
% arrnode = [2+51+101];
% reference node to calculate rigid body rotation
n_ref   = 2;
for i = 1:length(arrnode)
    node = arrnode(i);
    
    % initial tip position and directors
    inz = [12*(node-1)+1:12*(node-1)+12]; 
    d10 = q(1,inz(4:6))'/norm(q(1,inz(4:6)));
    d20 = q(1,inz(7:9))'/norm(q(1,inz(7:9)));
    d30 = q(1,inz(10:12))'/norm(q(1,inz(10:12)));
    x0  = q(1,inz(1:3))';

    % initial position and director of reference node
    inz_ref = [12*(n_ref-1)+1:12*(n_ref-1)+12]; 
    d1_ref0 = q(1,inz_ref(4:6))'/norm(q(1,inz_ref(4:6)));
    d2_ref0 = q(1,inz_ref(7:9))'/norm(q(1,inz_ref(7:9)));
    d3_ref0 = q(1,inz_ref(10:12))'/norm(q(1,inz_ref(10:12)));
    x_ref0  = q(1,inz_ref(1:3))';
    
    inztime = [1:1:size(time,1)];
    for j = 1:length(inztime)
        
        % position and director of tip node at t
        x  = q(inztime(j),inz(1:3))';
        d1 = q(inztime(j),inz(4:6))'/norm(q(j,inz(4:6)));
        d2 = q(inztime(j),inz(7:9))'/norm(q(j,inz(7:9)));
        d3 = q(inztime(j),inz(10:12))'/norm(q(j,inz(10:12)));
        
        % position and director of reference node at t
        d1_ref = q(inztime(j),inz_ref(4:6))'/norm(q(1,inz_ref(4:6)));
        d2_ref = q(inztime(j),inz_ref(7:9))'/norm(q(1,inz_ref(7:9)));
        d3_ref = q(inztime(j),inz_ref(10:12))'/norm(q(1,inz_ref(10:12)));
        x_ref  = q(inztime(j),inz_ref(1:3))';
        
        % rotation matrix between reference node at initial and current position
        R      = d1_ref*d1_ref0' + d2_ref*d2_ref0' + d3_ref*d3_ref0';
        
        % vector of tip deflection in director global cos
        u(j,:) = (x-x0)';

        % vector of tip deflection in director cos
        u_d0(j,1) = (x-x0)'*d10;
        u_d0(j,2) = (x-x0)'*d20;
        u_d0(j,3) = (x-x0)'*d30;
        
        % relative tip deflection of tip in blade cos
        u_rel_ref  = x - (x_ref + R*(x0-x_ref0));
        u_rel(j,1) = u_rel_ref'*d1;
        u_rel(j,2) = u_rel_ref'*d2;
        u_rel(j,3) = u_rel_ref'*d3;
    end
    if flag_type == 1
        figure(); hold on; grid on;
        title(['relative tip displacement' 'of blade ' num2str(i) ' in blade cos']);
        plot(time(inztime,1),u_rel(:,1),'-r','linewidth',2);
        plot(time(inztime,1),u_rel(:,2),'linestyle','-','color',[0.3027    0.8299    0.4687],'linewidth',2);
        plot(time(inztime,1),u_rel(:,3),'-b','linewidth',2);
        ylabel('tip-deflection in m'); xlabel('time in s')
        leg = legend('flap-wise','edge-wise','longitudinal');
        set(leg,'box','off'); pbaspect([2 1 1]);
        set(gca,'fontsize',14,'fontweight','bold');
        saveFigure(gcf,['wt_tip_rel_displ_blade_' num2str(i) '.fig']);
        print(['wt_tip_rel_displ_blade_' num2str(i)],'-dpng', '-r500');
    elseif flag_type == 2
        figure(); hold on; grid on; 
        title(['relative tip displacement' 'of blade ' num2str(i) ' in global cos']);
        plot(time(inztime,1),u(:,1),'-r','linewidth',2);
        plot(time(inztime,1),u(:,2),'linestyle','-','color',[0.3027    0.8299    0.4687],'linewidth',2);
        plot(time(inztime,1),u(:,3),'-b','linewidth',2);
        ylabel('tip-deflection in m'); xlabel('time in s')
        leg = legend('u_1','u_2','u_3');
        set(leg,'box','off'); pbaspect([2 1 1]);
        set(gca,'fontsize',14,'fontweight','bold');
        saveFigure(gcf,['wt_tip_gl_displ_blade_' num2str(i) '.fig']);
        print(['wt_tip_gl_displ_blade_' num2str(i)],'-dpng', '-r500');
    end
end
return

