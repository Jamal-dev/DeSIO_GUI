function [M] = uvlm_m_surface(model,t,qs,dp,x0)
% =================================================================================================================
    M = [];
    for k = 1:length(t)
        a = 0;
        mcp = [0,0,0];
        for i = 1:model.nsurfaces
            nodes = model.surfaces(i).nodes;
            for j = 1:model.surfaces(i).nelements
                c   = model.surfaces(i).connectivity(j,:);
                n1  = c(1); n2 = c(2); n3 = c(3); n4 = c(4);
                % Area of surface
                diag1 = qs(k,nodes(n1).indices_q)-qs(k,nodes(n3).indices_q);
                diag2 = qs(k,nodes(n2).indices_q)-qs(k,nodes(n4).indices_q);
                e1 = (diag2-diag1);
                e2 =-(diag1+diag2);
                e1 = e1/norm(e1);
                e2 = e2/norm(e2);
                nj = cross(e1, e2);
                area1 = 0.5d0*norm(cross(qs(k,nodes(n2).indices_q)-qs(k,nodes(n1).indices_q), qs(k,nodes(n4).indices_q)-qs(k,nodes(n1).indices_q)));
                area2 = 0.5d0*norm(cross(qs(k,nodes(n4).indices_q)-qs(k,nodes(n3).indices_q), qs(k,nodes(n2).indices_q)-qs(k,nodes(n3).indices_q)));
                Aj = area1 + area2;
                % coordinates of control point
                xcpj = 0.25*( qs(k,nodes(n1).indices_q) + qs(k,nodes(n2).indices_q) + qs(k,nodes(n3).indices_q) + qs(k,nodes(n4).indices_q) )';
                xaj  = xcpj - x0;
                % resultant force/dynamic pressure due to pressure
                Fj  = dp(k,j+a)*Aj*nj';
%                 mj  = cross(Fj,xaj)';
                mj  = cross(xaj,Fj)';
                mcp = mcp + mj;
            end
            % global lift = projecting Cp_vector in lift direction
            a = a + model.surfaces(i).nelements;
        end
        M(k,1:3) = mcp;
    end
% =================================================================================================================
return