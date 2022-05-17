[i_m, j_m, k_m] = size(jacket_line_end_points);
pnt_left = jacket_line_end_points(:,1,:);
pnt_right = jacket_line_end_points(:,2,:);
pnt_left = squeeze(pnt_left);
pnt_right = squeeze(pnt_right);
close all
k = 1;
scatter3(jacket_corrd(k,:), jacket_corrd(k+1,:), jacket_corrd(k+2,:))
for i = 1:i_m
    hold on 
    
    temp = [pnt_left(i,:); pnt_right(i,:)];
    plot3(temp(:,1),temp(:,2), temp(:,3))
    disp(i);
end
hold off