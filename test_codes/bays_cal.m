close all
clear all
%%
number_legs = 3;
height = 4;
losg = 1;
l_tp = 1;
r_base = 3;
r_top = 2;
number_bays = 2;
number_elements_per_beam = 2;
%%
oblique_angle_plane = atand(height/(r_base-r_top));
angle_pos = 360/number_legs;


r = @(h) r_base + (r_top-r_base)/(height - 0 ) * (h - 0);
L = sqrt((r_base-r_top)^2 + height^2);
%% calculating heights of stands
heights_stands = zeros( number_bays +2,1);
heights_stands(1) = losg;
height_comp_stand = (height - losg - l_tp)/(number_bays);
for i =1: number_bays
    heights_stands(i+1) = height_comp_stand;
end
heights_stands(i+2) = l_tp;
%% calculating lengths of stands
lengths_stands = heights_stands/height * L;
%% absolute height division by number of elements
absolute_height_radius= cell(size(heights_stands,1),1);
position_stands= cell(size(heights_stands,1),1);
h = 0;
for s=1:size(heights_stands,1)
    current_height =  heights_stands(s);
    local_heights = linspace(0,current_height, number_elements_per_beam + 1);
    
    abs_heights = (local_heights + h)';
    % absolute height
    h = h + current_height;
    abs_radius = r(abs_heights);
    absolute_height_radius{s,1} = [abs_heights, abs_radius];
    position_stands{s,1} = convert_hr_position(absolute_height_radius{s,1},number_legs);
end

origin = [ 0, 0, 0];
bays_pnts = calc_bay_pnts(position_stands,number_legs, number_bays, number_elements_per_beam);
plotJacket(position_stands,bays_pnts,number_legs)
[stands_data,bays_data] = findingDirectors(position_stands,bays_pnts,number_legs);

function [stands_data,bays_data] = findingDirectors(position_stands,bays_pnts,number_legs)

    number_stands = size(position_stands,1);
    
    stands_data = cell(number_stands, number_legs);
    for s = 1:number_stands
        for leg =1:number_legs
            beam_pnts = position_stands{s,1}(:,:,leg);
            stands_data{s,leg} = cal_dirs(beam_pnts,true);
        end
        
    end
    number_bays = size(bays_pnts{1,1},1);
    bays_data = cell(1,number_legs);
    for leg =1 :number_legs
        bays_data{1,leg} = cell(number_bays,2);
        for bay = 1:number_bays
            beam_pnts1 = bays_pnts{leg,1}{bay,1};
            beam_pnts2 = bays_pnts{leg,1}{bay,2};
            bays_data{1,leg}{bay,1} = cal_dirs(beam_pnts1);
            bays_data{1,leg}{bay,2} = cal_dirs(beam_pnts2);
        end
    end
    hold off
end

function result = cal_dirs(beam_pnts,if_plotting)
    if nargin ==1
        if_plotting = true;
    end
    number_pnts = size(beam_pnts,1);
    % calculating once because for a line it will remain same
    delta = beam_pnts(2,:) - beam_pnts(1,:);
    x = delta/norm(delta,2);
    y = [-x(2), x(1), 0];
    y = y/norm(y,2);
    z = cross(x,y);
    z = z/norm(z,2);
    x_dirs = repmat(x,number_pnts,1);
    y_dirs = repmat(y,number_pnts,1);
    z_dirs = repmat(z,number_pnts,1);
    result = [beam_pnts, x_dirs, y_dirs, z_dirs];
    if if_plotting
        hold on
        quiver3(beam_pnts(:,1),beam_pnts(:,2),beam_pnts(:,3),x_dirs(:,1),x_dirs(:,2),x_dirs(:,3))
        quiver3(beam_pnts(:,1),beam_pnts(:,2),beam_pnts(:,3),y_dirs(:,1),y_dirs(:,2),y_dirs(:,3))
        quiver3(beam_pnts(:,1),beam_pnts(:,2),beam_pnts(:,3),z_dirs(:,1),z_dirs(:,2),z_dirs(:,3))
    end
end

function plotJacket(position_stands,bays_pnts,number_legs)
    linewidth = 2.5;
    for s = 1: size(position_stands,1)
        stand = position_stands{s,1};
        for leg = 1: number_legs
            x = stand(:,1,leg);
            y = stand(:,2,leg);
            z = stand(:,3,leg);
            plot3(x,y,z,'LineWidth',linewidth);
            hold on
        end
        
    end
    for  k = 1:size(bays_pnts,1)
        leg_junc = bays_pnts{k};
        for bay = 1: size(leg_junc,1)
            line1 = leg_junc{bay,1};
            line2 = leg_junc{bay,2};
            x1 = line1(:,1);
            y1 = line1(:,2);
            z1 = line1(:,3);
            x2 = line2(:,1);
            y2 = line2(:,2);
            z2 = line2(:,3);
            plot3(x1,y1,z1,'LineWidth',linewidth)
            plot3(x2,y2,z2,'LineWidth',linewidth)
        end
    end
    hold off
    axis equal
end
function bays_pnts = calc_bay_pnts(position_stands,number_legs, number_bays, number_elements_per_beam)
    wrapN = @(i) (1 + mod(i-1, number_legs));
    % bays will be between leg1, leg2 ; leg2, leg1; leg3, leg1
    bays_pnts = cell(number_legs,1);
    for leg= 1:number_legs
        bays_pnts{leg} = cell(number_bays,2);
        l1 = leg;
        l2 = leg +1;
        l2 = wrapN(l2);
        for bay =1:number_bays
            temp = position_stands{bay+1,1};
            lower_pnt1 = temp(1,:,l1);
            lower_pnt2 = temp(1,:,l2);
            upper_pnt1 = temp(end,:,l1);
            upper_pnt2 = temp(end,:,l2);
            comp_1 = upper_pnt2 - lower_pnt1;
            comp_2 = upper_pnt1 - lower_pnt2;
            len_1 = norm(comp_1,2);
            len_2 = norm(comp_2,2);

            ratio_1 = linspace(0,len_1,number_elements_per_beam+1);
            ratio_2 = linspace(0,len_2,number_elements_per_beam+1);
            dir_1 = comp_1 /len_1;
            dir_2 = comp_2 /len_2;

            pnts_1 = line3D(lower_pnt1,ratio_1, dir_1);
            pnts_2 = line3D(lower_pnt2,ratio_2, dir_2);
            bays_pnts{leg}{bay,1} = pnts_1;
            bays_pnts{leg}{bay,2} = pnts_2;
        end
    end
end

function pnts = line3D(lower_pnt,mag,dir)
    pnts = zeros(length(mag),3);
    for pnt_num = 1:length(mag)
        pnts(pnt_num,:) = lower_pnt + mag(pnt_num) * dir;
    end
end
function pnts = convert_hr_position(mat,number_legs)
    heights = mat(:,1);
    radius = mat(:,2);
    number_pnts = size(heights,1);
    pnts = zeros(number_pnts,3,number_legs); % 3 is for x y z
    angle_pos = 360/number_legs;
    for k= 1: number_pnts
        h = heights(k,1);
        r = radius(k,1);
        for l=1:number_legs
            pnts(k, 1, l) = r * cosd(l *angle_pos);
            pnts(k, 2, l) = r * sind(l *angle_pos);
            pnts(k, 3, l) = h;
        end
    end
end