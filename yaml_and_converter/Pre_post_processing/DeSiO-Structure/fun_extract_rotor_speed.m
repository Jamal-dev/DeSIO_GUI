% script for calculating rotor speed from director-based formulation
function [res] = fun_extract_rotor_speed(q,s,t,node,nd); 

% indizes for node
inz = 12*(node-1)+1:12*(node-1)+12;

for i = 1:length(t)
    
    d1  = q(i,inz(4:6))'; 
    d2  = q(i,inz(7:9))'; 
    d3  = q(i,inz(10:12))';
    
    w1  = s(i,inz(4:6))'; 
    w2  = s(i,inz(7:9))'; 
    w3  = s(i,inz(10:12))';
    
    n   = nd(1)*d1 + nd(2)*d2 + nd(3)*d3;
    omega(i,1) = n'*0.5*(cross(d1,w1) + cross(d2,w2) + cross(d3,w3));
end

omega = omega.*(60/(2*pi));
res.val = omega;
return

