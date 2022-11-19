% Constitutive Laws
clc; clear all
cd 'H:\02_Elementherleitung\03_material_laws'

% Isotropic linear elastic materials
syms Lam nue E G
Lam = G*(E-2*G)/(3*G-E);
nue = G;
% Cijkl = Lam dij dkl + nue (dik*djl + dil*djk)
delta = eye(3);
for i = 1:3
    for j = 1:3
        for k = 1:3
            for l = 1:3
                C(i,j,k,l) = Lam*delta(i,j)*delta(k,l) + nue*(delta(i,k)*delta(j,l) + delta(i,l)*delta(j,k));
            end
        end
    end
end

% Beam with condition that s22 = s33 = 0
inz = {[1,1,1,1],[1,1,2,2],[1,1,3,3],[1,1,1,2],[1,1,1,3];...
       [2,2,1,1],[2,2,2,2],[2,2,3,3],[2,2,1,2],[2,2,1,3];...
       [3,3,1,1],[3,3,2,2],[3,3,3,3],[3,3,1,2],[3,3,1,3];...
       [1,2,1,1],[1,2,2,2],[1,2,3,3],[1,2,1,2],[1,2,1,3];...
       [1,3,1,1],[1,3,2,2],[1,3,3,3],[1,3,1,2],[1,3,1,3]};

% Rearranging to matrix notation   
for i = 1:size(inz,1)
    for j = 1:size(inz,2)
        Cb(i,j) = C(inz{i,j}(1),inz{i,j}(2),inz{i,j}(3),inz{i,j}(4));
    end
end

% Checking, if obtaining correct strains-stress relation
syms s11 s12 s13
b = [s11;0;0;1/2*s12;1/2*s13];
eps = inv(Cb)*b
return
