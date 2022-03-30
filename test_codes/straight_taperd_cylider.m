r1=5;
r2=10;
z1=5;
z2=10;
rad = @(z) r1 + (r2-r1)/(z2-z1) * z;
n_pnts = 25;
z_vec = linspace(z1,z2,n_pnts);
r_vec = rad(z_vec);
X=zeros(n_pnts,n_pnts);
Y=zeros(n_pnts,n_pnts);
Z=zeros(n_pnts,n_pnts);
theta_vec = linspace(0,2*pi,n_pnts);
for i=1:n_pnts
    z = z_vec(i);
    r = rad(z);
    x = r * cos(theta_vec);
    y = r * sin(theta_vec);
    Z(:,i) = z * ones(n_pnts,1);
    X(:,i) = x';
    Y(:,i) = y';
end
figure(1)
surf(X,Y,Z)




