function [res] = fun_extract_stress_resulants(stress,time,element)
% Function to extract dof-solution for node from solution.m file
for i = 1:length(element)
    res(i).F1 = stress(:,6*(element(i)-1)+1); % shear force in local 1
    res(i).F2 = stress(:,6*(element(i)-1)+2); % shear force in local 2
    res(i).F3 = stress(:,6*(element(i)-1)+3); % normal force in local 3
    res(i).M1 = stress(:,6*(element(i)-1)+4); % bending moment around local 1
    res(i).M2 = stress(:,6*(element(i)-1)+5); % bending moment around local 2
    res(i).M3 = stress(:,6*(element(i)-1)+6); % torsion moment around local 3
end
save('res_stress','res');
return