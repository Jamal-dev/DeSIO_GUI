function hplot = plot_triad(dx,dy,dz,d,d_scale)
dx = d_scale*dx; dy = d_scale*dy; dz = d_scale*dz;
hplot(1) = plot3([d(1), d(1)+dx(1)],[d(2), d(2)+dx(2)],[d(3), d(3)+dx(3)],'-r','linewidth',2);
hplot(2) = plot3([d(1), d(1)+dy(1)],[d(2), d(2)+dy(2)],[d(3), d(3)+dy(3)],'-g','linewidth',2);
hplot(3) = plot3([d(1), d(1)+dz(1)],[d(2), d(2)+dz(2)],[d(3), d(3)+dz(3)],'-b','linewidth',2);
