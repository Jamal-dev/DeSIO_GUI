% Matlab-Script zur Extrahierung von DeSiO-Knotenergebnissen
function [var] = get_DeSiO_dof_solu(node,Element,q)

for i = 1:length(node)
    node_i = node(i);
        
    % Filename and variable name to create
    filName = ['n' num2str(node_i) '_de' '.m'];  % zu erzeugende Datei mit Input als der Variablennamen
    varName = ['u_de'];                          % zu erzeugender Variablenname
    
    % Freiheitsgrad aus Knoten ermitteln
    if Element == 1
        ndof = 12*(node_i-1);
    elseif Element == 2
        ndof = 6*(node_i-1);
    end
    
    var = genvarname(varName);
    eval([var '(:,1)= q(:,ndof+1)-q(1,ndof+1);' ]);     % ux
    eval([var '(:,2)= q(:,ndof+2)-q(1,ndof+2);' ]);     % uy
    eval([var '(:,3)= q(:,ndof+3)-q(1,ndof+3);' ]);     % uz
    var = eval(var);
end