function res = fun_extract_wt_tip_deformation(q,time,arrnode,n_ref)
% =========================================================================
% Plotting DeSiO DOF-Solution
% =========================================================================
% flag for plot output
% flag_type(1) - relative tip displacement
% flag_type(2) - global tip displacement

% Function to extract dof-solution for node from solution.m file
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
    res(i).node = arrnode(i);
    res(i).u = u_rel;
end
save('res_u_tip','res');
return

