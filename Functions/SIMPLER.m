function simpler_output = SIMPLER(lambda_exc,NA,lambda_em,angle,...
    nS,nI,alpha,N0,filename_wformat,Illum,max_dist,...
    scatter_plot,camera_px,rz_xyz,file_format,hObject,handles,calib_file_name)
%% Calibration (dF & alphaF)

[dF, alphaF] = getParameters_SIMPLER(lambda_exc,NA,lambda_em,angle,nS,nI,alpha);

%% Importation
if file_format == 1 % Importation procedure for Picasso hdf5 files.
    % load input HDF5 file
    info = h5info(filename_wformat);
    data = h5read(filename_wformat,'/locs');
    frame = data.frame;
    photon_raw = data.photons;
    bg = data.bg;

    % convert x,y,sd values from 'camera subpixels' to nanometres
    xloc = data.x .* camera_px;
    yloc = data.y .* camera_px;
    sx = data.lpx .* camera_px;
    sy = data.lpy .* camera_px;
    sd = (sx.*2 + sy.*2) .^0.5;

elseif file_format == 2 % Importation procedure for ThunderStorm csv files.
    
   [c_frame c_x c_y c_ph c_off] = get_column_idx_ThunderStorm(filename_wformat);
   full_list_ThunderStorm = csvread(filename_wformat,1,0);
   frame = full_list_ThunderStorm(:,c_frame);
   xloc = full_list_ThunderStorm(:,c_x);
   yloc = full_list_ThunderStorm(:,c_y);
   photon_raw = full_list_ThunderStorm(:,c_ph);
   bg = full_list_ThunderStorm(:,c_off);
   
elseif file_format == 3 % Importation procedure for custom csv files.
   full_list = csvread(filename_wformat);
   frame = full_list(:,1);
   xloc = full_list(:,2);
   yloc = full_list(:,3);
   photon_raw = full_list(:,4);
   if size(full_list,2)>4
       bg = full_list(:,5);
   end
   
end
% Correct photon counts

% To perform this correction, the linear dependency between local laser
% power intensity and  local background photons is used. Photons are
% converted to the value they would have if the whole image was illuminated
% with the maximum laser power intensity. This section is executed if the
% user has chosen to perform correction due to non-flat illumination.


if Illum == 1
    profile = csvread(calib_file_name);
    phot=double(photon_raw);
    phot_corr=zeros(size(photon_raw,1),1);
    max_bg = prctile(profile(:),97.5);
    for i=1:size(phot,1)
        if (ceil(xloc(i)/camera_px) < size(profile,1)) &&...
             (ceil(yloc(i)/camera_px)< size(profile,2))
        phot_corr(i,1)=phot(i,1)*(double(max_bg)/(profile(ceil(xloc(i)...
            /camera_px),ceil(yloc(i)/camera_px))));
        elseif (ceil(xloc(i)/camera_px) > size(profile,1)) &&...
             (ceil(yloc(i)/camera_px)< size(profile,2))
         phot_corr(i,1)=phot(i,1)*(double(max_bg)/(profile(floor(xloc(i)...
            /camera_px),ceil(yloc(i)/camera_px))));
        elseif (ceil(xloc(i)/camera_px) < size(profile,1)) &&...
             (ceil(yloc(i)/camera_px)> size(profile,2))
         phot_corr(i,1)=phot(i,1)*(double(max_bg)/(profile(ceil(xloc(i)...
            /camera_px),floor(yloc(i)/camera_px))));
        elseif (ceil(xloc(i)/camera_px) > size(profile,1)) &&...
             (ceil(yloc(i)/camera_px)> size(profile,2))
         phot_corr(i,1)=phot(i,1)*(double(max_bg)/(profile(floor(xloc(i)...
            /camera_px),floor(yloc(i)/camera_px))));
        end
        
    end
else
    phot_corr = double(photon_raw);
end
    
% build the output array
listLocalizations = horzcat(xloc, yloc, double(frame),phot_corr);

[~,idx] = sort(listLocalizations(:,3)); % sort just the third column
sortedmat = listLocalizations(idx,:);   % sort the whole matrix using the 
                                        % 'sorted' indices

listLocalizations = unique(listLocalizations,'rows','stable');


%% Removing locs with no locs in i-1 and i+1 frame

% We keep a molecule if there is another one in the previous and next frame,
% located at a distance < max_dist, where max_dist is introduced by user
% (20 nm by default).

min_div=2.5*1e3; % We divide the list into sub-lists of 'min_div' locs, 
% to minimize memory usage;

Ntimes_min_div=ceil(size(listLocalizations,1)/min_div); 
truefalse_sum_roi_acum = 0;
listLocalizations_filtered=zeros(min_div*Ntimes_min_div,size(listLocalizations,2));

if (rz_xyz == 4) || (rz_xyz == 5)
    percent_k_molec = 0;
    f = waitbar(percent_k_molec,['Running frame filtering step...',...
        num2str(round(100*percent_k_molec)),' %']);
end

for N=1:Ntimes_min_div
    
    if (rz_xyz == 4) || (rz_xyz == 5)
        percent_k_molec = N/Ntimes_min_div;
    waitbar(percent_k_molec,f,['Running frame filtering step...',...
        num2str(round(100*percent_k_molec)),' %']);
    end
    
    min_div_N = min_div;
    
    if N == Ntimes_min_div
        min_div_N = size(listLocalizations,1)-min_div*(Ntimes_min_div-1);
    end
    
    truefalse=zeros(min_div_N,min_div_N);%This matrix relates each localization
                                       % with the rest of localizations. It
                                       % will take a value of 1 if the
                                       % molecules i and j (i =
                                       % row, j = column) are located at
                                       % distances < max_dist, and are
                                       % detected in frames N and (N+1) or 
                                       % (N-1). In any other case, it will
                                       % take a value of 0.
    daa=zeros(size(listLocalizations,1),1);
    frame_dif=zeros(size(listLocalizations,1),1);
    
    for i=(1+min_div_N*(N-1)):(min_div_N*N)
        for j=(1+min_div_N*(N-1)):(min_div_N*N)
            daa(j-(min_div_N*(N-1)),1)=((xloc(i,1)-xloc(j,1)).^2+...
                (yloc(i,1)-yloc(j,1)).^2).^(1/2);
            frame_dif(j-(min_div_N*(N-1)),1)=((listLocalizations(i,3)...
                -listLocalizations(j,3)).^2).^(1/2);
            if daa(j-(min_div_N*(N-1)),1)<max_dist && frame_dif(j-...
                    (min_div_N*(N-1)),1)==1
                truefalse(i-(min_div_N*(N-1)),j-(min_div_N*(N-1)))=1;
            end
        end
    end
    truefalse_sum=sum(truefalse')'; % For each row (i.e. each molecule)
                                    % we calculate the sum of every
                                    % column from the 'truefalse'
                                    % matrix.

    truefalse_sum_roi=roicolor(truefalse_sum,2,2); % We label with a '1'
                                                   % every molecule (row)
                                                   % with at least 2 other
                                                   % molecules located at
                                                   % d<max_dist and
                                                   % detected at
                                                   % frame+1 or frame-1
                                                   % (i.e. those rows
                                                   % whose columnwise
                                                   % sum is equal to 2)
    truefalse_sum_roi=double(truefalse_sum_roi);

    truefalse_sum_roi_acum = [truefalse_sum_roi_acum;...
        truefalse_sum_roi]; % We save the molecules' labels of each loop
                            % in a the '_acum' vector.                             
end

if (rz_xyz == 4) || (rz_xyz == 5)
    close(f)
end

idx_filtered = find(truefalse_sum_roi_acum(2:end)); % idx_filtered saves the
                                                    % indices of the kept
                                                    % molecules.
       
 %% Applying filter to original list
 
frame1 = listLocalizations(idx_filtered,3);
x1 = listLocalizations(idx_filtered,1);
y1 = listLocalizations(idx_filtered,2);
photons1 = listLocalizations(idx_filtered,4);

   
   
%% Clear plots and turn off their visibility.

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



%% Z calculation and scatter plots

function mycallback(o,e) % This function is used to set the marker
                         % size of the scatter plots when the zoom is
                         % re-adjusted,in order to set a default value 
                         % of ~ 4 nm.

        xlimits = xlim;
        h11 = findobj(gcf,'Type','scatter');
        if sum(size(h11))>0
            h11.SizeData = 50*250/(xlimits(2)-xlimits(1));
        end
end

% Small ROI and Large ROI cases

if rz_xyz == 2 | rz_xyz == 3 | rz_xyz == 4
    z1 = (log(alphaF*N0)-log(photons1-(1-alphaF)*N0))/(1/dF);
    z1=real(z1);


    % Scatter plots / Projection r vs. z, x vs. z , y vs. z

    if rz_xyz == 2 % If the selected operation is "(r,z) Small ROI", we 
                   % perform a first step where the main axis is obtained
                   % from a linear fit of the (x,y) data.

        x=x1; y=y1; z=z1;
        M=[x y z];
        P=polyfit(x,y,1); 
        yfit = P(1)*x+P(2);
        Poly_fun = @(x) P(1)*x+P(2);


        Origin_X = min(x);
        Origin_Y = Poly_fun(Origin_X);

        tita=atan(P(1));
        tita1=atan((y-Origin_Y)./(x-Origin_X));
        r=((x-Origin_X).^2+(y-Origin_Y).^2).^(1/2);
        tita2=tita-tita1;
        proyec_r=cos(tita2).*r;

        simpler_output = [x1 y1 proyec_r z1 photons1 frame1]; % The generated output
                                                       % includes the
                                                       % (x,y,r,z,photons,frame)
                                                       % data.


        if scatter_plot == 1 % The following commands are executed if the
                             % scatter plot option is enabled.
            axes(handles.axes1); 
            set(handles.open_fig_1_tag,'visible', 'on'); 
            set(handles.export_csv_1_tag,'visible', 'on');
            set(handles.export_csv_whole_1_tag,'visible', 'on');
            set(handles.slider_1_tag,'visible', 'on');  
            set(handles.slider_1_text_tag,'visible', 'on'); 
            set(handles.slider_1_tag,'Value',0.66);



            h1 = scatter(proyec_r,z,40*ones(size(z,1),1),z,'o','filled');
            xlabel('r (nm)');
            ylabel('z (nm)');
            daspect([1 1 1])
            pbaspect([1 1 1])
            set(handles.axes3,'visible', 'off');
            set(handles.open_fig_3_tag,'visible', 'off');
            set(handles.export_csv_3_tag,'visible', 'off');
            set(handles.export_csv_whole_3_tag,'visible', 'off');
            xlimits_inicial = xlim;
            h1.SizeData = 50*250/(xlimits_inicial(2)-xlimits_inicial(1));
            h = zoom;
            h.ActionPostCallback = @mycallback;

        end



    elseif rz_xyz == 3 % If the selected operation is '(x,y,z)Small roi',
                       % the scatter plots correspond to x vs. z and y vs.z
                       % and the output matrix will include (x,y,z,photons)
                       % information.
        z = z1;
        if scatter_plot== 1
            axes(handles.axes1);
            set(handles.open_fig_1_tag,'visible', 'on'); 
            set(handles.export_csv_1_tag,'visible', 'on');
            set(handles.export_csv_whole_1_tag,'visible', 'off');
            set(handles.slider_1_tag,'visible', 'on'); 
            set(handles.slider_3_tag,'visible', 'on'); 
            set(handles.slider_1_text_tag,'visible', 'on'); 
            set(handles.slider_3_text_tag,'visible', 'on'); 
            set(handles.slider_1_tag,'Value',0.66);
            set(handles.slider_3_tag,'Value',0.66);


            h1 = scatter(x1,z,50*ones(size(z,1),1),z,'o','filled');
            xlabel('x (nm)');
            ylabel('z (nm)');
            daspect([1 1 1])
            pbaspect([1 1 1])
            xlimits_inicial = xlim;
            h1.SizeData = 50*250/(xlimits_inicial(2)-xlimits_inicial(1));
            h = zoom;
            h.ActionPostCallback = @mycallback;

            axes(handles.axes3);
            set(handles.open_fig_3_tag,'visible', 'on'); 
            set(handles.export_csv_3_tag,'visible', 'on');
            set(handles.export_csv_whole_3_tag,'visible', 'on');

            h2 = scatter(y1,z,50*ones(size(z,1),1),z,'o','filled');
            xlabel('y (nm)');
            ylabel('z (nm)');
            daspect([1 1 1])
            pbaspect([1 1 1])
            xlimits_inicial = xlim;
            h2.SizeData = 50*250/(xlimits_inicial(2)-xlimits_inicial(1));
            h = zoom;
            h.ActionPostCallback = @mycallback;

        end
        simpler_output = [x1 y1 z1 photons1 ones(length(x1),1) frame1]; % The output
                                                                 % for this 
                                                                 % operation includes
                                                                 % (x,y,z,photons,frame)
                                                                 % information.
    elseif rz_xyz == 4 % Large ROI case
        simpler_output = [x1 y1 z1 photons1 ones(length(x1),1) frame1];
        % The output matrix includes (x,y,z,photons,frame) information.
        if scatter_plot== 1
            axes(handles.axes7);
            set(handles.axes7,'visible','on');
            set(handles.slider_1_tag,'visible', 'on'); 
            set(handles.slider_1_text_tag,'visible', 'on'); 
            set(handles.slider_1_tag,'Value',0.66);
            set(handles.plot_xz_tag,'visible', 'on'); 
            set(handles.plot_yz_tag,'visible', 'on'); 
            set(handles.set_small_roi_tag,'visible', 'on'); 
            set(handles.roi_angle_text_tag,'visible', 'on'); 
            set(handles.angle_ROI_tag,'visible', 'on'); 
            set(handles.angle_ROI_tag,'Value',0.5);
            set(handles.angle_ROI_tag,'enable','off');
            set(handles.render_tag,'enable','off');
            h1 = scatter(x1,y1,50*ones(size(z1,1),1),z1,'o','filled');
            xlabel('x (nm)');
            ylabel('y (nm)');
            daspect([1 1 1])
            pbaspect([1 1 1])
            xlimits_inicial = xlim;
            h1.SizeData = 50*250/(xlimits_inicial(2)-xlimits_inicial(1));
            h = zoom;
            h.ActionPostCallback = @mycallback;
            

            
        end
        
    end
end
% N0 Calibration case

if rz_xyz == 5 % For the "N0 Calibration" operation, there is no "Z 
               % calculation", because the aim of this procedure is to
               % obtain N0 from a sample which is supposed to contain
               % molecules located at z ~ 0.
               
    xl = [min(x1) max(x1)]; yl = [min(y1) max(y1)];
    c = 1:length(x1);
    axes(handles.axes3);
    set(handles.axes3,'visible', 'on'); 
    set(handles.open_fig_3_tag,'visible', 'on');  
    set(handles.axes1,'visible', 'of'); 
    set(handles.open_fig_1_tag,'visible', 'of'); 
   
   if exist('histogram') >0
       histogram(photons1(c)); % N0 histogram plot
       h = findobj(gcf, 'Type', 'histogram');
       b = h.Values;
       L = h.BinLimits;
       a = L(1):(L(2)-L(1))/(length(b)-1):L(2);

   else
       hist(photons1(c));
       [b a] = hist(photons1(c));
       L = [a(1) a(end)];
   end
    xlimit = L;
    xl_fit = xlimit;   
    f = fit(a',b','gauss1'); % Gaussian fit of the N0 distribution
    
    MyCoeffs = coeffvalues(f);
    a1 = MyCoeffs(1); b1 = MyCoeffs(2); c1 = MyCoeffs(3);    
    sigma = c1/sqrt(2);
    x_plot = xl_fit(1):xl_fit(2);
    y_plot = f(x_plot);
    hold on, plot(x_plot,y_plot,'LineWidth',3);
    xlabel('N_0 (photons)');
    ylabel('Counts');
    title({'N_0 histogram',['N_0 = ',num2str(round(b1)), ...
        ' photons; \sigma_N_0 = ',num2str(round(sigma)),' photons'],...
        ['Area = ', num2str(round(abs(1E-6*(xl(2)-xl(1))*(yl(2)-yl(1))))),...
        ' \mum^2']});    

    N0 = b1;
    z1 = (log(alphaF*N0)-log(photons1-(1-alphaF)*N0))/(1/dF);
    z1=real(z1);

    simpler_output = [x1 y1 z1 photons1 zeros(length(x1),1) frame1];

end
num_locs_orig = size(listLocalizations,1);
num_locs_final = length(x1);


if rz_xyz < 4   
    f = msgbox(['Number of localizations in raw data = ', ...
        num2str(num_locs_orig),'; Number of localizations post-frame filtering step = ',...
        num2str(num_locs_final)]);
elseif rz_xyz == 4
    f = msgbox(['Number of localizations in raw data = ', ...
        num2str(num_locs_orig),'; Number of localizations post-frame filtering step = ',...
        num2str(num_locs_final),'. Please wait until (x,y,z) data is exported to .csv']);
else
    f = msgbox(['Number of localizations in raw data = ', ...
        num2str(num_locs_orig),'; Number of localizations post-frame filtering step = ',...
        num2str(num_locs_final),'. Please wait until N0 calibration results',...
        ' are exported to .csv']);
end


end



