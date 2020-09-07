function export_function(roi,simpler_output,handles)
% First, the function checks which was the last run operation by analyzing 
% the fifth column of 'simpler_output' array.
if sum(simpler_output(:,5)) == size(simpler_output,1) 
    rz_xyz = 3; % Last operation run: Small ROI (x,y,z) or Large ROI
elseif sum(simpler_output(:,5)) > size(simpler_output,1) 
    rz_xyz = 2; % Last operation run: Small ROI (r,z)
elseif sum(simpler_output(:,5)) == 0
    rz_xyz = 5; % Last operation run: N0 Calibration
else
    rz_xyz = 6; % Last operation run: Incidence angle and alpha adjustment
end 

filename_wformat = (get(handles.filename_tag,'string')); % Filename

if rz_xyz == 2 % Last operation run: Small ROI (r,z)
    
    if roi == 0 % 'Export whole dataset' option
        c = 1:size(simpler_output,1); 
    else % 'Export current view' option
        axes(handles.axes1); 
        xl = xlim; % 'x' limits are obtained from the scatter plot
        yl = ylim; % idem for 'y' limits
        c = find(simpler_output(:,3)>xl(1) & simpler_output(:,3)<xl(2) &...
            simpler_output(:,4)>yl(1) & simpler_output(:,4)<yl(2)); 
        % 'c' represents the indices of the molecules shown in the current
        % view of the scatter plot
    end
    x1 = simpler_output(c,1); % x positions of the filtered molecules
    y1 = simpler_output(c,2); % y positions of the filtered molecules
    proyec_r = simpler_output(c,3); % r positions of the filtered molecules
    frame = simpler_output(c,6); % frames of the filtered molecules
    photons = simpler_output(c,5); % photons of the filtered molecules
    z = simpler_output(c,4); % z positions of the filtered molecules
    exportation_file_rz(x1,y1,proyec_r,z,frame,photons,filename_wformat,handles); % This function
                                                % performs the exportation

elseif rz_xyz == 3 % Last operation run: Small ROI (x,y,z) or Large ROI
    
     if roi == 0 % 'Export whole dataset' option
        c1 = 1:size(simpler_output,1);
     elseif roi == 1 % 'Export current view' option (upper figure)
    axes(handles.axes1); 
    xl = xlim; % 'x' limits are obtained from the scatter plot
    yl = ylim; % idem for 'y' limits
    c1 = find(simpler_output(:,1)>xl(1) & simpler_output(:,1)<xl(2) &...
        simpler_output(:,3)>yl(1) & simpler_output(:,3)<yl(2));
    % 'c1' represents the indices of the molecules shown in the current
        % view of the scatter plot
     elseif roi == 3 % 'Export current view' option (lower figure)
    axes(handles.axes3); 
    xl = xlim; % 'x' limits are obtained from the scatter plot
    yl = ylim; % idem for 'y' limits
    c1 = find(simpler_output(:,2)>xl(1) & simpler_output(:,2)<xl(2) &...
        simpler_output(:,3)>yl(1) & simpler_output(:,3)<yl(2));
    % 'c1' represents the indices of the molecules shown in the current
        % view of the scatter plot
     end    
    x1 = simpler_output(c1,1); % x positions of the filtered molecules
    y1 = simpler_output(c1,2); % y positions of the filtered molecules
    z = simpler_output(c1,3); % z positions of the filtered molecules
    frame = simpler_output(c1,6); % frames of the filtered molecules
    photons = simpler_output(c1,4); % photons of the filtered molecules
    exportation_file_xyz(x1,y1,z,frame,photons,filename_wformat,handles); % This function
                                                % performs the exportation
                                                

elseif rz_xyz == 5 % Last operation run: N0 calibration
    
    z = simpler_output(:,3); % z positions
    x1 = simpler_output(:,1); % x positions
    y1 = simpler_output(:,2); % y positions
    photons = simpler_output(:,4); % photons
    frame = simpler_output(:,6); % frames of the filtered molecules
    exportation_file_xyz(x1,y1,z,frame,photons,filename_wformat,handles); % Exportation
    csvwrite('N0_list.csv',photons); % N0 list exportation  
    
elseif rz_xyz == 6 % Last operation run: Incidence angle and alpha adjustment
    % Data is obtained from the scatter plot of combined structures
    axes(handles.axes5);
    h = findobj(gcf,'Type','scatter');
    z = h.YData'; % relative axial position (z)
    x1 = ones(length(z),1); % As the lateral information is 'relative-r position' 
                            % from different structures, there is no x
                            % position information available. For that
                            % reason, a vector of 'ones' is used. 
    y1 = ones(length(z),1);% As the lateral information is 'relative-r position' 
                            % from different structures, there is no y
                            % position information available. For that
                            % reason, a vector of 'ones' is used. 
    proyec_r = h.XData'; % relative lateral position (r)
    frame = y1; % Since there are different structures from different measurements
                % involved, the frame information will not be exported and
                % for that reason a 'ones' vector is assigned to this
                % variable.
    photons = y1; % Idem
    exportation_file_rz(x1,y1,proyec_r,z,frame,photons,filename_wformat,handles); % Exportation   
end
end
