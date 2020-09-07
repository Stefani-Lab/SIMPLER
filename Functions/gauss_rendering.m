function gauss_rendering(Magnif,sigma_width_lateral_nm,sigma_width_axial_nm,...
    camera_px,proyec_r,z,xl,yl,rz_xyz,z_color,auto_contrast)

% This function builds a 2D matrix where each localizations is plotted as a
% normalized Gaussian function.

c_roi = find(proyec_r>xl(1) & proyec_r<xl(2) & z>yl(1) & z<yl(2)); % 'c_roi'
                                            % determines the indices of the
                                            % molecules located within the
                                            % ROI selected with the zoom
                                            % tool.

proyec_r = proyec_r(c_roi); % filtered lateral position

z = z(c_roi); % filtered axial position (or lateral, for (x,y) plot)


% Next, string references to the corresponding axis are generated. These
% strings are used in the title of the Figure showing the rendered image 

axial_ax = 'z';

if rz_xyz == 2 
    lateral_ax = 'r';
elseif rz_xyz == 3
    lateral_ax = 'x';
elseif rz_xyz == 4
    lateral_ax = 'y';
elseif rz_xyz == 5
    lateral_ax = 'x';
    axial_ax = 'y';   
end



px_size_rendering = camera_px/Magnif; % pixel size (nm) in the SR image



proyec_r_reorigin = proyec_r-xl(1)+sigma_width_lateral_nm; % relative lateral positions
z1_reorigin = z-yl(1)+sigma_width_axial_nm; %relative axial (or lateral, for (x,y) plot) positions



max_r = (xl(2)-xl(1))+sigma_width_lateral_nm; % Maximum relative lateral position
max_z = (yl(2)-yl(1))+sigma_width_axial_nm; % Maximum relative axial (or lateral) position


% Convertion of r and z from nm to px

SMr = proyec_r_reorigin/px_size_rendering;
SMz = z1_reorigin/px_size_rendering;
sigma_width_lateral_px = sigma_width_lateral_nm/px_size_rendering;
sigma_width_axial_px = sigma_width_axial_nm/px_size_rendering;

% Definition of pixels affected by the list of SM (+- 5*sigma in both
% dimensions):
A = zeros(ceil(max_r/px_size_rendering),ceil(max_z/px_size_rendering));

for i = 1:length(SMr)
    A(ceil(SMr(i)),ceil(SMz(i))) = 1;
end

sigma_width_nm = max([sigma_width_lateral_nm sigma_width_axial_nm]);
sigma_width_px = sigma_width_nm/px_size_rendering;
sd = strel('square',round(5*sigma_width_px));
A_affected = imdilate(A,sd); % This matrix contains 1 in +- 5 sigma units 
                             % around the SM positions

[r_affected z_affected] = find(A_affected); % r and z positions affected

PSF = @(r,z,SMr,SMz,I) (I/(2*pi*sigma_width_lateral_px...
    *sigma_width_axial_px))*exp(-((((r-SMr)^2)/...
    (2*sigma_width_lateral_px^2))+((z-SMz)^2)/(2*sigma_width_axial_px^2)));
% 'PSF' is a function that calculates the value for a given position (r,z)
% assuming a 2D Gaussian distribution centered at (SMr(k),SMz(k))

B = zeros(ceil(5*sigma_width_px+(max_r/px_size_rendering))...
    ,ceil(5*sigma_width_px+(max_z/px_size_rendering))); % B = empty matrix 
                                                % that will be filled with
                                                % the values given by the
                                                % Gaussian blur of the
                                                % points listed in (SMr,
                                                % SMz)

percent_k_molec = 0;
f = waitbar(percent_k_molec,['Loading gaussian rendering...',...
    num2str(round(100*percent_k_molec)),' %']);

for k = 1:length(SMr) % For each molecule from the list...
    percent_k_molec = k/length(SMr);
    waitbar(percent_k_molec,f,['Loading gaussian rendering...',...
        num2str(round(100*percent_k_molec)),' %']);

    for i = 1:length(r_affected) % For each pixel of the final image with a 
                                 % value different from zero...
        B(r_affected(i),z_affected(i)) = B(r_affected(i),z_affected(i))...
            + PSF(r_affected(i),z_affected(i),SMr(k),SMz(k),1);
        % Each 'affected' pixel (i) will take the value that had at the beggining
        % of the k-loop + the value given by the distance to the
        % k-molecule.
    end

    
end
close(f)

BB = B;

if auto_contrast == 1 % If the 'auto-contrast' option is enabled, the 'imadjust'
                      % function is applied. 
    BB_adj = imadjust(BB);
elseif auto_contrast == 0
    BB_adj = BB;
end


sb_length_candidates = [10 25 50 100 500 1000]; % 'Possible scalebar lengths'
sb_length_min = min(abs(0.25*max([size(B,2) size(B,1)])*...
    px_size_rendering-sb_length_candidates)); % Which 'possible scalebar length'
                                             % is closest to
                                             % 0.25*image_width?
sb_length_idx = find(abs(0.25*max([size(B,2) size(B,1)])*...
    px_size_rendering-sb_length_candidates) == sb_length_min);
sb_length = sb_length_candidates(sb_length_idx(1)); % % Selected scalebar length

figure, imshow(flipud(BB_adj'),'InitialMagnification','fit'); % Display image
daspect([1 1 1]); % Data aspect ratio: 1:1
xlabel(lateral_ax,'FontSize',15)
ylabel(axial_ax,'FontSize',15)
colormap(hot);

if rz_xyz < 5
title({[axial_ax,' vs ',lateral_ax,'; Pixel size (SR image) = ',...
    num2str(px_size_rendering),' nm; ',axial_ax,' range = ',...
    num2str(round(yl(1))),' - ',num2str(round(yl(2))),' nm'],...
    ['Scalebar = ',num2str(sb_length),' nm']});
elseif rz_xyz == 5
    title({[axial_ax,' vs ',lateral_ax,'; Pixel size (SR image) = ',...
    num2str(px_size_rendering),' nm; ',axial_ax,' range = ',...
    num2str(round(yl(1))),' - ',num2str(round(yl(2))),' nm'],...
    ['Scalebar = ',num2str(sb_length),' nm'],['(Full ROI)']});
end



width = round(sb_length/px_size_rendering); 
height = round(0.1*width); 
xCenter = size(BB_adj',2)-round(1.5*0.5*width); 
yCenter = size(BB_adj',1)-round(1.5*0.5*height); 
xLeft = xCenter - width/2;
yBottom = yCenter - height/2;
rectangle('Position', [xLeft, yBottom, width, height],...
    'EdgeColor', 'k', 'FaceColor', 'w', 'LineWidth', 0.5);
h = zoom(); 
h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]); % After performing a zoom-in,
                                                      % the data aspect
                                                      % ratio remains 1:1


if z_color == 1 % If z-color rendering is enabled...
    
    BBB = render_colored_final(BB_adj',1/px_size_rendering);
    figure, imshow(flipdim(BBB,1),'InitialMagnification','fit');
    daspect([1 1 1]); % Data aspect ratio: 1:1
    xlabel(lateral_ax,'FontSize',15)
    ylabel(axial_ax,'FontSize',15)
    title({[axial_ax,' vs ',lateral_ax,'; Pixel size (SR image) = ',...
         num2str(px_size_rendering),' nm; ',axial_ax,' range = ',...
            num2str(round(yl(1))),' - ',num2str(round(yl(2))),' nm'],...
                ['Scalebar = ',num2str(sb_length),' nm']});

    
    width = round(sb_length/px_size_rendering); 
    height = round(0.1*width); 
    xCenter = round(0.95*size(BBB,2))-round(1.5*0.5*width); 
    yCenter = round(size(BBB,1))-round(1.5*0.5*height); 
    xLeft = xCenter - width/2;
    yBottom = yCenter - height/2;
    rectangle('Position', [xLeft, yBottom, width, height],...
        'EdgeColor', 'k', 'FaceColor', 'w', 'LineWidth', 0.5);
    h = zoom();
    h.ActionPostCallback = @(o, e) daspect(e.Axes, [1 1 1]); % After performing a zoom-in,
                                                      % the data aspect
                                                      % ratio remains 1:1

    
    
end



end