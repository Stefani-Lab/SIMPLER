function B = render_colored_final(A,z_rend_prec)
% This function creates a z-color-coded matrix

% LUT:
mapname{1,1} = 'magenta';
mapname{2,1} ='mediumorchid';
mapname{3,1} ='darkorchid';
mapname{4,1} ='mediumpurple';
mapname{5,1} ='mediumslateblue';
mapname{6,1} ='slateblue';
mapname{7,1} ='deepskyblue';
mapname{8,1} ='mediumturquoise';
mapname{9,1} ='mediumspringgreen';
mapname{10,1} ='lightgreen';
mapname{11,1} ='yellow';

maprgb = zeros(size(mapname,1),3);
for i = 1:size(mapname,1)
    maprgb(i,:) = rgb(mapname{i,1});
end


interp_range = 300; % 'pseudo-continuum' LUT (300 intervals)
maprgb_interp(:,1) = interp1((1:size(maprgb,1))',...
    maprgb(:,1),(1:(1/interp_range):size(maprgb,1))'); 
maprgb_interp(:,2) = interp1((1:size(maprgb,1))',...
    maprgb(:,2),(1:(1/interp_range):size(maprgb,1))'); 
maprgb_interp(:,3) = interp1((1:size(maprgb,1))',...
    maprgb(:,3),(1:(1/interp_range):size(maprgb,1))'); 

res_px = ceil(z_rend_prec);
num_of_changes = round(size(A,1)/res_px); % Number of z-levels
colormap_chosen = maprgb_interp(1:floor(size(maprgb_interp,1)...
    /(num_of_changes)):size(maprgb_interp,1),:);
colormap_chosen = colormap_chosen(1:num_of_changes,:); % Colors for each z-level

sat = 180;
sat_white = [sat 0];
colormap_chosen_final = cell(size(colormap_chosen,1),1); 
for i = 1:size(colormap_chosen,1)
    % For each color (z-level), a transition from black to bright is built
    colormap_chosen_final{i,1} = [linspace(0,colormap_chosen(i,1),...
        sat_white(1))',linspace(0,colormap_chosen(i,2),sat_white(1))',...
        linspace(0,colormap_chosen(i,3),sat_white(1))';...
        linspace(colormap_chosen(i,1),0.5,sat_white(2))',...
        linspace(colormap_chosen(i,2),0.5,sat_white(2))',...
        linspace(colormap_chosen(i,3),0.5,sat_white(2))';...
        linspace(colormap_chosen(i,1),1,256-sum(sat_white))',...
        linspace(colormap_chosen(i,2),1,256-sum(sat_white))',...
        linspace(colormap_chosen(i,3),1,256-sum(sat_white))']; 
end



A = [A (sat/255)*ones(size(A,1),round(0.05*size(A,2)))]; % A vertical stripe is added to A (on the right side)
x = [1 size(A,2)];
y = 1:res_px:size(A,1);
if y(end)<size(A,1)
    y = [y size(A,1)];
end
grayImage = mat2gray(A); % grayImage is a 'gray' image built from A

grayImage1 = round(double(255*grayImage/(max(max(grayImage)))));

 croppedImage = cell((size(y,2)-1),1);
 rgbCroppedImage = croppedImage;
for i = 1:(size(y,2)-1)
    % Horizontal stripes representing z steps are added to the grayImage.
    % The colors of these stripes are extracted from colormap_chosen_final.
    % As different z fragments give rise to different colors, they end up
    % forming a z-color code.
    croppedImage{i,1} = grayImage1(floor(y(i)):floor(y(i+1)), x(1):x(2));
    rgbCroppedImage{i,1} =double(ind2rgb(croppedImage{i,1},...
        colormap_chosen_final{i,1}));
end
    
rgbImage = cat(3, grayImage, grayImage, grayImage);


for i = 1:(size(y,2)-1)
    rgbImage(floor(y(i)):floor((y(i+1))), x(1):x(2), :) = rgbCroppedImage{i,1};
end

B = rgbImage;
end
    
    