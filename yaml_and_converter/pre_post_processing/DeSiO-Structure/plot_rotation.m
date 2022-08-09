% =========================================================================
% plotting DeSiO solution of rotations
% =========================================================================
clc;
clear all;
close all;

model = strucure_readmodel;

% loading of result files
q    = load([[model.strSimName '_q.dres']);
time = load([[model.strSimName '_t.dres']); time = time(:,1);

% determine rotations of beam nodes
node = [2];
for i = 1:length(node)
    figure(); hold on; grid on;
    d1_0 = q(1,12*(node(i)-1)+4  : 12*(node(i)-1)+6);
    d2_0 = q(1,12*(node(i)-1)+7  : 12*(node(i)-1)+9);
    d3_0 = q(1,12*(node(i)-1)+10 : 12*(node(i)-1)+12);
    phi = zeros(length(time),3);
    for j = 2:length(time)
        d1_n = q(j,12*(node(i)-1)+4  : 12*(node(i)-1)+6);
        d2_n = q(j,12*(node(i)-1)+7  : 12*(node(i)-1)+9);
        d3_n = q(j,12*(node(i)-1)+10 : 12*(node(i)-1)+12);
        delta_d1 = d1_n - d1_0; delta_d2 = d2_n - d2_0; delta_d3 = d3_n - d3_0;
        phi(j,1:3) = phi(j-1,1:3) + 0.5*( cross(d1_n,delta_d1) + cross(d2_n,delta_d2) + cross(d3_n,delta_d3) );
        d1_0 = d1_n; d2_0 = d2_n; d3_0 = d3_n;
    end
    plot(time,phi(:,1),'-','lineWidth',2,'color','r');
    plot(time,phi(:,2),'-','lineWidth',2,'color','g');
    plot(time,phi(:,3),'-','lineWidth',2,'color','b');
    ylabel('rotation [rad]'); xlabel('time [s]')
    title(['rotation vs. time of node' num2str(node(i))]);
    legend('\phi_1','\phi_2','\phi_3');
    set(gca,'fontsize',12,'fontweight','bold');
end

return