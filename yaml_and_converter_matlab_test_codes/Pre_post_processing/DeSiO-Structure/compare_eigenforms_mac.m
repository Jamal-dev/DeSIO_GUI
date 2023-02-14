function compare_eigenforms_mac()
clc
clear all
close all

addpath(genpath('U:\PhD-2021-12_Christian_Hente\Pre_post_processing'));

strfilename = 'test';

scale = 1.0;

% Read simulation 1
curDir = cd;
model1 = structure_readmodel;
[ef1,ev1] = get_eigenforms(model1,1.0);

% Read simulation 2
curDir = cd;
model2 = structure_readmodel;
[ef2,ev2] = get_eigenforms(model2,1.0);

% Modal Assurance Criterion
mac=MAC(ef1,ef2)

% Plot MAC-Data
f1 = figure(); b = bar3(mac); title('MAC');
for k = 1:length(b); zdata = b(k).ZData; b(k).CData = zdata; b(k).FaceColor = 'interp'; end; 
colormap('jet');
colorbar;
axis([0 size(mac,1)+1 0 size(mac,2)+1 0 1]);
savefig(f1,[strfilename '_eigenform_MAC.fig']);
% print([strfilename '_eigenform_MAC'],'-dpng','-r0');
return

function [ef,ev] = get_eigenforms(model,scale)
    step  = fun_load_file([model.strSimName '_steps.dres']);
    inzm  = find(step(:,2)==3);
    t     = fun_load_file([model.strSimName '_t.dres']);
    q     = fun_load_file([model.strSimName '_q.dres']); 
    q0    = q(1,:); dq = q(inzm,:);

    inz   = [(1:12:size(dq,2))',(2:12:size(dq,2))',(3:12:size(dq,2))'];
    inz   = reshape(inz', [], 1);

    ev = t(inzm,1);
    ef = [];
    for j = 1:length(ev)
        ef(j,inz) = q0(inz) + dq(j,inz)*scale';
    end
   return
   
function mac=MAC(phi1,phi2)
    % phi: matrix of the identified mode shapes
    % mac: MAC matrix
    for I=1:size(phi1,1)
        for J=1:size(phi2,1)
            mac(I,J)=Mac(phi1(I,:),phi2(J,:));
        end
    end
return

function mAc=Mac(Phi1,Phi2)
    % This function calculates mac between phi1 and phi2
    mAc= (abs(dot(Phi1,Phi2)))^2/((dot(Phi1,Phi1))*(dot(Phi2,Phi2)));
return

