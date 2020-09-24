function output_ft_input = fine_tuning(filename_wformat,handles)

    % clear plots

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
    set(handles.render_tag,'enable','on');
    set(handles.axes5,'visible', 'off'); 
    set(handles.axes6,'visible', 'off'); 
    set(handles.auto_ft_panel_tag,'visible', 'off'); 
    set(handles.calc_exc_prof_tag,'visible', 'off'); 
    set(handles.fit_relz_gauss_tag,'visible', 'on'); 
    set(handles.z_bin_tag,'visible', 'on'); 
    set(handles.z_bin_txt_tag,'visible', 'on'); 
    set(handles.hist_settings_panel_tag,'visible', 'on'); 
    set(handles.slider_1_tag,'visible', 'off'); 
    set(handles.slider_3_tag,'visible', 'off'); 
    set(handles.slider_1_text_tag,'visible', 'off'); 
    set(handles.slider_3_text_tag,'visible', 'off'); 
    set(handles.slider_1_tag,'visible', 'off'); 
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



    % Read csv file containing (lateral,axial) positions of known
    % structures.
    data = csvread(filename_wformat);
    lateral_matrix = data(:,1:2:end); % Lateral positions are obtained from
                                      % odd columns.
    axial_matrix = data(:,2:2:end);   % Axial positions are obtained from
                                      % even columns.
    
    
    lateral = lateral_matrix(:); % Lateral positions are unified in one vector.
    axial = axial_matrix(:); % The same is done for axial positions.
    

    % The values used for the axial positions calculation of the known
    % structures through SIMPLER are obtained from the 'run_SIMPLER'
    % interface. It is important to set them correctly before running the
    % operation.
    
    set(handles.N0_input_ft_tag,'string',get(handles.N0_tag,'string'));
    set(handles.alpha_input_ft_tag,'string',get(handles.alpha_tag,'string'));
    set(handles.angle_input_ft_tag,'string',get(handles.angle_tag,'string'));
    set(handles.excwave_input_ft_tag,'string',get(handles.exc_lambda_tag,'string'));
    set(handles.emwave_input_ft_tag,'string',get(handles.em_wavelength_tag,'string'));
    set(handles.ni_input_ft_tag,'string',get(handles.ni_tag,'string'));
    set(handles.ns_input_ft_tag,'string',get(handles.ns_tag,'string'));
    set(handles.na_input_ft_tag,'Value',get(handles.pop_up_NA_tag,'Value'));

    % Once the parameters are read from the graphical user interface, the
    % dF and alphaF values are calculated.


    N0_input = str2num((get(handles.N0_input_ft_tag,'string')));
    alpha_input = str2num((get(handles.alpha_input_ft_tag,'string')));
    angle_input = str2num((get(handles.angle_input_ft_tag,'string')));
    lambda_exc_input = str2num((get(handles.excwave_input_ft_tag,'string')));
    lambda_em_input =str2num((get(handles.emwave_input_ft_tag,'string')));
    nI_input = str2num((get(handles.ni_input_ft_tag,'string')));
    nS_input = str2num((get(handles.ns_input_ft_tag,'string')));
    y8 = (get(handles.na_input_ft_tag,'Value'));
    NA_pop = y8;
    NA_input = 1.42*(NA_pop==1)+1.45*(NA_pop==2)+1.49*(NA_pop==3);
    [dF_input, alphaF_input] = getParameters_SIMPLER(lambda_exc_input,NA_input,lambda_em_input,...
        angle_input,nS_input,nI_input,alpha_input);

    % The number of photons for each localization are retrieved from the 
    % axial position and the dF, alphaF and N0 values obtained in the
    % above step.
    
    photons_matrix = zeros(size(axial_matrix,1),size(axial_matrix,2));
    photons_median_matrix = zeros(size(axial_matrix,1),size(axial_matrix,2));
    lateral_median_matrix = zeros(size(axial_matrix,1),size(axial_matrix,2));

    % The next function allows to obtain a custom 'median' value, which is
    % calculated as the mean value between the p-centile 10% and p-centile
    % 90% from a given distribution. We use this function in order to
    % re-center the localizations from the known structures around [0,0]
    
    median_perc90_10_center = @(A) mean([prctile(A,90) prctile(A,10)]);

    for i = 1:size(lateral_matrix,2)
        c = find(axial_matrix(:,i)); % To calculate the 'median' values, we
                                     % use valid localizations, i.e. those
                                     % with axial positions different from 0;
                                     % there are elements filled with z = 0 
                                     % and lateral = 0 in the 'data' matrix,
                                     % because not every known structures
                                     % have the same number of
                                     % localizations.
        photons_matrix(:,i) = alphaF_input*N0_input*exp(-axial_matrix(:,i)...
            /dF_input)+(1-alphaF_input)*N0_input;
        photons_median_matrix(:,i) = ones(size(photons_matrix,1),1)*...
            median_perc90_10_center(photons_matrix(c,i));
        lateral_median_matrix(:,i) = ones(size(lateral_matrix,1),1)*...
            median_perc90_10_center(lateral_matrix(c,i));
    end

    photons = photons_matrix(:); % Number of photons for each localization
    
    photons_median = photons_median_matrix(:); % Median value of the 
                                               % number of photons for the
                                               % structure to which each
                                               % localization belongs
                                               
    lateral_median = lateral_median_matrix(:); % Lateral positions median 
                                               % value for the structure to
                                               % wich each localization
                                               % belongs
                                               
    axial_median = (log(alphaF_input*N0_input)-log(photons_median-(1-...
        alphaF_input)*N0_input))/(1/dF_input); % Median value of the 
                                               % axial position for the
                                               % structure to which each
                                               % localization belongs

    % Some elements from the axial vector are zero because the 'data' matrix
    % contains structures with different number of localizations and thus
    % there are columns filled with zeros. Now, we remove those elements. 
    
    c = find(axial == 0);

    photons(c) = [];
    photons_median(c) = [];

    lateral_median(c) = [];
    lateral(c) = [];
    axial_median(c) = [];
    axial(c) = [];
    axial = (log(alphaF_input*N0_input)-log(photons-(1-alphaF_input)*N0_input))/(1/dF_input);

    % The output for this function will contain the information of the
    % parameters used to retrieve the number of photons from the axial
    % positions; the number of photons and lateral position of every
    % localization; and the median values for both the number of photons and
    % the lateral positions.
    % These outputs will be used by the 'Update Scatter' buttom as
    % input information, in order to recalculate every axial position when
    % either the angle or alpha (or both) are changed.
    
    output_ft_input(:,1) = [lambda_exc_input; NA_input; lambda_em_input;...
        nS_input; nI_input; photons];
    output_ft_input(:,2) = [ones(5,1); lateral];
    output_ft_input(:,3) = [ones(5,1); photons_median];
    output_ft_input(:,4) = [ones(5,1); lateral_median];

    % Scatter plot of the relative (lateral,axial) positions (centered in
    % [0,0])

    axes(handles.axes5); 
    set(handles.axes5,'visible', 'on'); 
    set(handles.panel_ft_input_tag,'visible', 'on');
    set(handles.ft_panel_tag,'visible', 'on'); 
    set(handles.auto_ft_panel_tag,'visible', 'on'); 


    scatter((lateral-lateral_median),(axial-axial_median),5)
    xlabel('Lateral (nm)');
    ylabel('Axial (nm)'); 
    title('Combined structures');

    box on
    daspect([1 1 1])
    zlimits = [min((axial-axial_median))-0.1*(max(axial-axial_median)-...
        min(axial-axial_median)) max((axial-axial_median))+0.1*(max(axial...
        -axial_median)-min(axial-axial_median))]; %'zlimits' establishes an axial 
                                           % range centered at the axial
                                           % 'median' position, and spanning
                                           % over an axial length 20% greater
                                           % than the range given by
                                           % (max-min)axial positions.
    zlimits_axis = ylim;
    set(handles.axes5, 'Units', 'pixels');
    pos_box = plotboxpos(handles.axes5); % This function will let us match
                                         % the axis of the relative-z
                                         % histogram with the axial axis of
                                         % the 'Combined structures'
                                         % scatter plot.
    
    
    
    lat_lim = xlim;
    set(handles.lateral_range_3_tag,'string',num2str(round(lat_lim(2))));
    set(handles.lateral_range_1_tag,'string',num2str(round(lat_lim(1))));

    % Plot as vertical lines the lateral limits selected for the relative-z
    % histogram.
        if str2num(get(handles.lateral_range_3_tag,'string'))...
                >str2num(get(handles.lateral_range_1_tag,'string'))
            hold on, plot([str2num(get(handles.lateral_range_1_tag,...
                'string')) str2num(get(handles.lateral_range_1_tag,'string'))],...
            [zlimits(1)+5 zlimits(2)-5],'LineStyle',':','LineWidth',1.5,'Color','k');
            hold on, plot([str2num(get(handles.lateral_range_3_tag,...
                'string')) str2num(get(handles.lateral_range_3_tag,'string'))],...
            [zlimits(1)+5 zlimits(2)-5],'LineStyle',':','LineWidth',1.5,'Color','k');      
        end


    % Histogram of axial positions
    axes(handles.axes4);
    set(handles.axes4,'visible', 'on'); 
    if exist('histogram') >0
        histogram(axial);
    else
        hist(axial);
    end
    xlabel('z-position (nm)');
    title('z histogram');
    axes(handles.axes6);
    set(handles.axes6,'visible', 'on'); 

    set(handles.lateral_range_1_tag,'visible', 'on'); 
    set(handles.lateral_range_2_tag,'visible', 'on'); 
    set(handles.lateral_range_3_tag,'visible', 'on');

    set(handles.lateral_range_1_tag,'string',num2str(round(min(lateral-lateral_median))));
    set(handles.lateral_range_3_tag,'string',num2str(round(max(lateral-lateral_median))));

    % Histogram of relative-axial positions
    nbins = zlimits(1):str2num(get(handles.z_bin_tag,'string')):zlimits(2);

    if exist('histogram') >0
        histogram((axial-axial_median),nbins);
    else
       hist((axial-axial_median),nbins);
    end
        set(gca,'xlim',zlimits_axis); 
        set(gca,'view',[90 -90]) % Rotate histogram
        set(handles.axes6, 'Units', 'pixels');
            pos_box1 = get(handles.axes6,'Position'); 
            pos_box1(4) = pos_box(4);% Match of relative-z histogram z-axis with 
                                     % the corresponding axial axis from the scatter
                                     % plot.
            pos_box1(2) = pos_box(2);
            set(handles.axes6,'Position',pos_box1);
            title('relative-z histogram');


    set(handles.export_csv_whole_1_tag,'visible', 'on'); 
    set(handles.yn_circfit_tag,'visible', 'on'); 
    set(handles.r_fit_tag,'visible', 'on'); 
    set(handles.circ_fit_tag,'visible', 'on'); 

end
