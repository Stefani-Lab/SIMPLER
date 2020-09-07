function update_axial(lateral,photons,dF_ft,alphaF_ft, N0_input,...
        handles,photons_median,lateral_median)
    
    axial_ft = (log(alphaF_ft*N0_input)-log(photons-(1-alphaF_ft)...
    *N0_input))/(1/dF_ft); % The updated axial positions ('axial_ft') are calculated
                           % from the input parameters introduced in the
                           % 'alpha-angle adjustment' panel
    axial_ft = real(axial_ft);
    axial_ft_median = (log(alphaF_ft*N0_input)-log(photons_median-...
        (1-alphaF_ft)*N0_input))/(1/dF_ft); % The axial 'median' position
                                            % is also updated.
    axial_ft_median = real(axial_ft_median);

    % Scatter plot of updated (lateral, axial) positions
    axes(handles.axes6); 
        cla reset;
    axes(handles.axes5); 
        cla reset;
    scatter((lateral-lateral_median),(axial_ft-axial_ft_median),5)
    box on
    xlabel('Lateral (nm)');
    ylabel('Axial (nm)');
    title('Combined structures');
    zlimits = [min((axial_ft-axial_ft_median))-0.1*(max(axial_ft-...
            axial_ft_median)-min(axial_ft-axial_ft_median))...
         max((axial_ft-axial_ft_median))+0.1*(max(axial_ft-axial_ft_median)...
         -min(axial_ft-axial_ft_median))]; %'zlimits' establishes an axial 
                                           % range centered at the axial
                                           % 'median' position, and spanning
                                           % over an axial length 20% greater
                                           % than the range given by
                                           % (max-min)axial positions.
                            
      
    zlimits_axis = ylim;

    daspect([1 1 1])
    pos_box = plotboxpos(handles.axes5); % This function will let us match
                                         % the axis of the relative-z
                                         % histogram with the axial axis of
                                         % the 'Combined structures'
                                         % scatter plot.

    
    % Histogram of axial positions
    
    axes(handles.axes4);  
    
    if exist('histogram') >0
        histogram(axial_ft); % Histogram of updated axial positions
    else
        hist(axial_ft);
    end
    xlabel('z-position (nm)');
    title('z histogram');

    axes(handles.axes6);
    set(handles.axes6,'visible', 'on'); 
    set(handles.lateral_range_1_tag,'visible', 'on'); 
    set(handles.lateral_range_2_tag,'visible', 'on'); 
    set(handles.lateral_range_3_tag,'visible', 'on');

    lat_ax = [(lateral-lateral_median) (axial_ft-axial_ft_median)]; 
    % lat_ax represents the (lateral,axial) relative positions, centered
    % at [0,0]
    
    % Plot as vertical lines the lateral limits selected for the relative-z
    % histogram.
    if str2num(get(handles.lateral_range_3_tag,'string'))>...
            str2num(get(handles.lateral_range_1_tag,'string'))
        c1 = find(lat_ax(:,1)>str2num(get(handles.lateral_range_1_tag,'string')) &...
            lat_ax(:,1)<str2num(get(handles.lateral_range_3_tag,'string')));
        axes(handles.axes5); 
        hold on, plot([str2num(get(handles.lateral_range_1_tag,'string'))...
            str2num(get(handles.lateral_range_1_tag,'string'))],...
            [zlimits(1)+5 zlimits(2)-5],'LineStyle',':','LineWidth',1.5,'Color','k');
         hold on, plot([str2num(get(handles.lateral_range_3_tag,'string'))...
             str2num(get(handles.lateral_range_3_tag,'string'))],...
            [zlimits(1)+5 zlimits(2)-5],'LineStyle',':','LineWidth',1.5,'Color','k');       
    else
        c1 = 1:size(lat_ax,1);
    end
    axes(handles.axes6); 
    nbins = zlimits(1):str2num(get(handles.z_bin_tag,'string')):zlimits(2);
    if exist('histogram') >0
        histogram(lat_ax(c1,2),nbins); % Histogram of relative-z positions
    else
       hist(lat_ax(c1,2),nbins);
    end
    set(gca,'xlim',zlimits_axis);
    set(gca,'view',[90 -90]) % Rotate histogram
    set(handles.axes6, 'Units', 'pixels');
    pos_box1 = get(handles.axes6,'Position'); 
    pos_box1(4) = pos_box(4); % Match of relative-z histogram z-axis with 
                              % the corresponding axial axis from the scatter
                              % plot.
    pos_box1(2) = pos_box(2);
    set(handles.axes6,'Position',pos_box1);
    title('relative-z histogram');


    

    
