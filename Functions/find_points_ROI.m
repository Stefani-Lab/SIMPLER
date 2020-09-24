function [c,slope, origin] = find_points_ROI(xy, angle_roi, positions_roi,x_or_y)

x_inf_positions_roi = positions_roi(1);
y_inf_positions_roi = positions_roi(2);
w_roi = positions_roi(3);
h_roi = positions_roi(4);
x_y_cent = [x_inf_positions_roi+0.5*w_roi y_inf_positions_roi+0.5*h_roi];

if (angle_roi ~= 0) && (angle_roi~= pi)
    y1_roi = x_y_cent(2)+sin(angle_roi)*(h_roi*0.5-w_roi*0.5/tan(angle_roi));
    x1_roi = x_y_cent(1)+((y1_roi-x_y_cent(2))/tan(angle_roi))+0.5*w_roi/sin(angle_roi);
elseif angle_roi == 0
    y1_roi = x_y_cent(2)-w_roi*0.5;
    x1_roi = x_y_cent(1)+h_roi*0.5;
else
    y1_roi = x_y_cent(2)+w_roi*0.5;
    x1_roi = x_y_cent(1)-h_roi*0.5;   
end

x2_roi = x1_roi-cos(angle_roi)*h_roi;
y2_roi = y1_roi-sin(angle_roi)*h_roi;
x3_roi = x2_roi-sin(angle_roi)*w_roi;
y3_roi = y2_roi+cos(angle_roi)*w_roi;
x4_roi = x1_roi-sin(angle_roi)*w_roi;
y4_roi = y1_roi+cos(angle_roi)*w_roi;

xy_roi = [x1_roi y1_roi; x2_roi y2_roi; x3_roi y3_roi; x4_roi y4_roi];

[~,idx] = sort(xy_roi(:,1)); 
sorted_xy_roi = xy_roi(idx,:); 
                                     
slope_1 = abs((sorted_xy_roi(3,2)-sorted_xy_roi(4,2))/(sorted_xy_roi(3,1)-sorted_xy_roi(4,1))); 
if sorted_xy_roi(3,2) > sorted_xy_roi(4,2)
    slope_1 = -slope_1;
end
intercept_1 = sorted_xy_roi(3,2)-slope_1*sorted_xy_roi(3,1);
slope_2 = abs((sorted_xy_roi(3,2)-sorted_xy_roi(1,2))/(sorted_xy_roi(3,1)-sorted_xy_roi(1,1))); 
if sorted_xy_roi(3,2) < sorted_xy_roi(1,2)
    slope_2 = -slope_2;
end
intercept_2 = sorted_xy_roi(3,2)-slope_2*sorted_xy_roi(3,1);
slope_3 = slope_1;
intercept_3 = sorted_xy_roi(2,2)-slope_3*sorted_xy_roi(2,1);
slope_4 = slope_2;
intercept_4 = sorted_xy_roi(2,2)-slope_4*sorted_xy_roi(2,1);

f1 = @(x) slope_1*x+intercept_1;
f2 = @(x) slope_2*x+intercept_2;
f3 = @(x) slope_3*x+intercept_3;
f4 = @(x) slope_4*x+intercept_4;

% x = [sorted_xy_roi(1,1) sorted_xy_roi(4,1)];


if sorted_xy_roi(2,2) > sorted_xy_roi(3,2)
    c = find(xy(:,2) > f1(xy(:,1)) & xy(:,2) > f2(xy(:,1))...
        & xy(:,2) < f3(xy(:,1)) & xy(:,2) < f4(xy(:,1)));
else
    c = find(xy(:,2) > f3(xy(:,1)) & xy(:,2) > f4(xy(:,1))...
        & xy(:,2) < f1(xy(:,1)) & xy(:,2) < f2(xy(:,1)));
end

% distance_1_2 = pdist([x1_roi y1_roi; x2_roi y2_roi]);
% distance_2_3 = pdist([x2_roi y2_roi; x3_roi y3_roi]);

distance_1_2 = pdist([sorted_xy_roi(1,1) sorted_xy_roi(1,2);...
    sorted_xy_roi(2,1) sorted_xy_roi(2,2)]);
distance_2_4 = pdist([sorted_xy_roi(2,1) sorted_xy_roi(2,2);...
    sorted_xy_roi(4,1) sorted_xy_roi(4,2)]);

if x_or_y == 'x'
    if distance_1_2 > distance_2_4
        slope = slope_1;
        P2 = [sorted_xy_roi(3,1) sorted_xy_roi(3,2)];
    else
        slope = slope_2;
        P2 = [sorted_xy_roi(2,1) sorted_xy_roi(2,2)];
    end
else
    if distance_1_2 > distance_2_4
        slope = slope_2;
        P2 = [sorted_xy_roi(2,1) sorted_xy_roi(2,2)];

    else
        slope = slope_1;
        P2 = [sorted_xy_roi(3,1) sorted_xy_roi(3,2)];
    end
end
P1 = [sorted_xy_roi(1,1) sorted_xy_roi(1,2)];
origin = (P1(:) + P2(:)).'/2;

if (1.5*pi-angle_roi == 0) || (1.5*pi-angle_roi == pi) || (1.5*pi-angle_roi == 2*pi)
    y1_roi = x_y_cent(2)-w_roi*0.5;
    x1_roi = x_y_cent(1)-h_roi*0.5;
    y2_roi = x_y_cent(2)+w_roi*0.5;
    x2_roi = x_y_cent(1)+h_roi*0.5;
    
    c = find(xy(:,2) > y1_roi & xy(:,2) < y2_roi & xy(:,1) > x1_roi & xy(:,1) < x2_roi);
    
    if x_or_y == 'x'
        comp_x_y = abs(x2_roi-x1_roi)>(abs(y2_roi-y1_roi));
        slope = 0*comp_x_y + atan(pi/2)*(1-comp_x_y);
        P2 = [x1_roi y1_roi];
        P1 = [x1_roi y2_roi];
        P1 = [x2_roi y2_roi];
        origin = (P1(:) + P2(:)).'/2;
    else
        comp_x_y = abs(x2_roi-x1_roi)>(abs(y2_roi-y1_roi));
        slope = 0*(1-comp_x_y) + atan(pi/2)*comp_x_y;
        P2 = [x1_roi y1_roi];
        P1 = [x2_roi y1_roi]; 
        P1 = [x2_roi y2_roi];
        origin = (P1(:) + P2(:)).'/2;
    end
        
    
elseif (1.5*pi-angle_roi == 0.5*pi) || (1.5*pi-angle_roi == 1.5*pi)
    y1_roi = x_y_cent(2)-h_roi*0.5;
    x1_roi = x_y_cent(1)-w_roi*0.5;
    y2_roi = x_y_cent(2)+h_roi*0.5;
    x2_roi = x_y_cent(1)+w_roi*0.5;
    
    c = find(xy(:,2) > y1_roi & xy(:,2) < y2_roi & xy(:,1) > x1_roi & xy(:,1) < x2_roi);
    if x_or_y == 'x'
        comp_x_y = abs(x2_roi-x1_roi)>(abs(y2_roi-y1_roi));
        slope = 0*comp_x_y + atan(pi/2)*(1-comp_x_y);
        P2 = [x1_roi y1_roi];
        P1 = [x1_roi y2_roi];
        P1 = [x2_roi y2_roi];
        origin = (P1(:) + P2(:)).'/2;
    else
        comp_x_y = abs(x2_roi-x1_roi)>(abs(y2_roi-y1_roi));
        slope = 0*(1-comp_x_y) + atan(pi/2)*comp_x_y;
        P2 = [x1_roi y1_roi];
        P1 = [x2_roi y1_roi]; 
        P1 = [x2_roi y2_roi];
        origin = (P1(:) + P2(:)).'/2;
    end
    
end
    




  

