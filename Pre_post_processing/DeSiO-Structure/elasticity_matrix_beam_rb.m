function to_call()
    clear all; clc; close all;
    
    mass_matrix_rigidbody();
    elasticity_matrix_beam();
    return
end

function elasticity_matrix_beam()
    syms E A G Ay Az Ip Iy Iz Sy Sz Syz da rho
    syms xhi1 xhi2 xhi3

% Isotropic Elasticity Matrix for beam in plane stress state    
    CBeam = [2*G,  0, 0;
               0,2*G, 0;
               0,  0, E];

    e1 = [1,0,0]'; e2 = [0,1,0]'; e3 = [0,0,1]';
% Parameter vector    
    teta = xhi1*e1+xhi2*e2;
% Skew matrix of parameter vector    
    Theta = skew(teta);
   
    A  = [sym(eye(3,3)),-Theta];
    At = [sym(eye(3,3));-transpose(Theta)];

% Beam elasticity matrix for computation in our notation
    CB  = At*CBeam*A*da
    CTT = CB(1:3,1:3)
    CTK = CB(1:3,4:6)
    CKT = CB(4:6,1:3)
    CKK = CB(4:6,4:6)

% mass matrix
    teta = [1;xhi1*e1+xhi2*e2];
    m = teta(1:3)*rho*transpose(teta(1:3))*da
 return
end

function mass_matrix_rigidbody()
    syms rho dv
    syms xhi1 xhi2 xhi3

    e1 = [1,0,0]'; e2 = [0,1,0]'; e3 = [0,0,1]';
% Parameter vector    
    teta = [1;xhi1*e1+xhi2*e2+xhi3*e3];
% mass matrix
    m = teta(1:4)*rho*transpose(teta(1:4))*dv
 return
end

function A = skew(a)
    A      = sym(zeros(3,3));
    A(1,2) = -a(3);
    A(1,3) =  a(2);
    A(2,1) =  a(3);
    A(2,3) = -a(1);
    A(3,1) = -a(2);
    A(3,2) =  a(1);
end