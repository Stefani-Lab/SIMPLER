function obtain_bg(handles)
% Clear plots
axes(handles.axes1); 
cla reset;
axes(handles.axes3); 
cla reset;
    axes(handles.axes4); 
cla reset;
    axes(handles.axes5); 
cla reset;
    axes(handles.axes6); 
cla reset;
cla reset;
    axes(handles.axes7); 
cla reset;
    axes(handles.axes8); 
cla reset;
set(handles.axes1,'visible', 'off');
set(handles.open_fig_1_tag,'visible', 'off');
set(handles.axes3,'visible', 'off');
set(handles.open_fig_3_tag,'visible', 'off');
set(handles.export_csv_1_tag,'visible', 'off');
set(handles.export_csv_3_tag,'visible', 'off');
set(handles.export_csv_whole_1_tag,'visible', 'off');
set(handles.export_csv_whole_3_tag,'visible', 'off');
set(handles.panel_ft_input_tag,'visible', 'off');
set(handles.ft_panel_tag,'visible', 'off'); 
set(handles.axes4,'visible', 'off'); 
set(handles.axes5,'visible', 'off'); 
set(handles.axes6,'visible', 'off');
set(handles.lateral_range_1_tag,'visible', 'off'); 
set(handles.lateral_range_2_tag,'visible', 'off'); 
set(handles.lateral_range_3_tag,'visible', 'off'); 
set(handles.auto_ft_panel_tag,'visible', 'off'); 
set(handles.fit_relz_gauss_tag,'visible', 'off'); 
set(handles.z_bin_tag,'visible', 'off'); 
set(handles.z_bin_txt_tag,'visible', 'off'); 
set(handles.yn_circfit_tag,'visible', 'off'); 
set(handles.r_fit_tag,'visible', 'off'); 
set(handles.circ_fit_tag,'visible', 'off'); 
set(handles.hist_settings_panel_tag,'visible', 'off'); 
set(handles.calc_exc_prof_tag,'visible', 'off'); 
set(handles.slider_1_tag,'visible', 'off'); 
set(handles.slider_3_tag,'visible', 'off'); 
set(handles.slider_1_text_tag,'visible', 'off'); 
set(handles.slider_3_text_tag,'visible', 'off'); 
set(handles.axes7,'visible', 'off'); 
set(handles.axes8,'visible', 'off'); 
set(handles.set_small_roi_tag,'visible', 'off'); 
set(handles.roi_angle_text_tag,'visible', 'off'); 
set(handles.angle_ROI_tag,'visible', 'off'); 
set(handles.plot_xz_tag,'visible', 'off'); 
set(handles.plot_yz_tag,'visible', 'off'); 
set(handles.marker_size_roi_text_tag,'visible', 'off'); 
set(handles.slider_markersize_roi_tag,'visible', 'off'); 
set(handles.export_roi_tag,'visible', 'off'); 


% Filename, file format and camera pixel size

filename_wformat = get(handles.filename_tag,'string');
camera_px = str2double(get(handles.px_size_tag,'string'));
file_format = get(handles.fileformat_tag,'value');

%% Importation
if file_format == 1 % Importation procedure for Picasso hdf5 files.
    % load input HDF5 file
    data = h5read(filename_wformat,'/locs');
    bg = data.bg;

    % convert x,y,sd values from 'camera subpixels' to nanometres
    xloc = data.x .* camera_px;
    yloc = data.y .* camera_px;
    
elseif file_format == 2 % Importation procedure for ThunderStorm csv files.
    
   [~, c_x, c_y, ~, c_off] = get_column_idx_ThunderStorm(filename_wformat);
   full_list_ThunderStorm = csvread(filename_wformat,1,0);
   xloc = full_list_ThunderStorm(:,c_x);
   yloc = full_list_ThunderStorm(:,c_y);
   bg = full_list_ThunderStorm(:,c_off);
   
elseif file_format == 3 % Importation procedure for custom csv files.
   full_list = csvread(filename_wformat);
   xloc = full_list(:,2);
   yloc = full_list(:,3);
   if size(full_list,2)>4
       bg = full_list(:,5);
   end
   
end

%% Excitation profile calculation and exportation

Img_bg = zeros(ceil(max(xloc/camera_px)),ceil(max(yloc/camera_px))); % Empty matrix
                                                % to be filled with
                                                % background values
                                                % pixel-wise.
                                                
Count_molec_px = Img_bg; % Empty matrix to be filled with the number of molecules
                         % used to calcualte local background for each pixel.


% If the list contains > 1,000,000 localizations, only 500,000 are used in
% order to speed up the analysis and avoid redundancy.
vector_indices = 1:length(xloc);
c_random = randperm(length(vector_indices));
length_c_random = min([1E6 length(xloc)]);
c = c_random(1:length_c_random);
                         
                         
percent_k_molec = 0;

if length_c_random == 1E6
    f = waitbar(percent_k_molec,['Generating excitation profile from ',...
        'SM local background (1,000,000 localizations used)...',...
    num2str(round(100*percent_k_molec)),' %']); % Waitbar message.
else
    f = waitbar(percent_k_molec,['Generating excitation profile from SM local background...',...
    num2str(round(100*percent_k_molec)),' %']); % Waitbar message.
end


for i = 1:length(xloc(c))    
    percent_k_molec = i/length(c);

    if length_c_random == 1E6
        waitbar(percent_k_molec,f,['Generating excitation',...
            ' profile from SM local background (1,000,000 localizations used)...',...
        num2str(round(100*percent_k_molec)),' %']); 
    else
        waitbar(percent_k_molec,f,['Generating excitation profile from SM local background...',...
        num2str(round(100*percent_k_molec)),' %']);    
    end

    Count_molec_px(ceil(xloc(c(i))/camera_px),ceil(yloc(c(i))/camera_px)) = ...
    Count_molec_px(ceil(xloc(c(i))/camera_px),ceil(yloc(c(i))/camera_px))+1; % The #molecules 
                                                % used to calculate local
                                                % background in current
                                                % pixel is updated.

    Count_i = Count_molec_px(ceil(xloc(c(i))/camera_px),ceil(yloc(c(i))/camera_px));

    if Count_i > 1 % If the current pixel's background has already been estimated, then its 
                 % value is updated through a weighted average.
    Img_bg(ceil(xloc(c(i))/camera_px),ceil(yloc(c(i))/camera_px)) = ...
        ((Count_i-1)/Count_i)*mean(Img_bg(ceil(xloc(c(i))/camera_px),...
        ceil(yloc(c(i))/camera_px)))+(1/Count_i)*bg(c(i));
    else
        Img_bg(ceil(xloc(c(i))/camera_px),ceil(yloc(c(i))/camera_px)) = bg(c(i));
    end
end 
close(f)
if file_format == 1
    N = 5;
else
    N = 4;
end
    
csvwrite(['excitation_profile_',filename_wformat(1:(end-N)),'.csv'],Img_bg);

end