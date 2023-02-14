% script for calculating external moments from director-based formulation
clc; close all; clear all;

model = strucure_readmodel;

q = load([model.strSimName '_q.dres']);
f = load([model.strSimName '_fext.dres']);
t = load([model.strSimName '_t.dres']); t = t(:,1);

node = 2;

% indizes for node
inz = 12*(node-1)+1:12*(node-1)+12;

for i = 1:length(t)
    q_node = q(i,inz);
    f_node = f(i,inz);
    m_node(i,:) = -( cross(f_node(4:6),q_node(4:6)) + cross(f_node(7:9),q_node(7:9)) + cross(f_node(10:12),q_node(10:12)) );
end

plot(t,m_node(:,1))
return