function gauss_rendering_function(Magnif,sigma_width_nm,camera_px,simpler_output,z_color,...
        handles,axref,auto_contrast)
    
    % This function first discerns which was the last operation run, by checking
    % the values of the fifth column of the 'simpler_output' column.
    % Depending on the operation 
    
       x1 = simpler_output(:,1); % x position
       y1 = simpler_output(:,2); % y position
            
        if sum(simpler_output(:,5)) > size(simpler_output,1) % If this is true,
                                                         % the last
                                                         % operation run
                                                         % could be 'Small
                                                         % ROI (r,z)' or 
                                                         % 'Incidence angle
                                                         % and alpha
                                                         % adjustment'.
            rz_xyz_bis = 2;
            if axref == 1 % Last operation run = 'Small ROI (r,z)'
            axes(handles.axes1); 
            elseif axref == 5 % Last operation run = 'Incidence angle and 
                              % alpha adjustment'
                axes(handles.axes5);
            end

            g = msgbox(['Select from the scatter',...
                ' plot the ROI to be renderized (w/zoom tool).',...
                ' Then, click "ok" button from this message box']);
            waitfor(g);
            xl = xlim; % lateral limits (r) are extracted from the scatter plot after zoom in
            yl = ylim; % axial limits (z) are extracted from the scatter plot after zoom in
            
            proyec_r = simpler_output(:,3); % lateral positions (r)
            z = simpler_output(:,4); % axial positions (z)
            gauss_rendering(Magnif,sigma_width_nm(1),sigma_width_nm(2),...
                camera_px,proyec_r,z,xl,yl,rz_xyz_bis,z_color,auto_contrast);
            % 'gauss_rendering.m' builds the matrix-image (Gaussian
            % rendering)


        elseif size(axref,1) > 1 % Last operation run = 'Large ROI'
                axes(handles.axes8); 
                rz_xyz_bis = 2;          
            proyec_r = axref(:,2); % lateral positions (r)
            xl = [min(proyec_r)-1 max(proyec_r)+1];
            z = simpler_output(axref(:,1),3); % axial positions (z)
            yl  = ylim;
            gauss_rendering(Magnif,sigma_width_nm(1),sigma_width_nm(2),...
                camera_px,proyec_r,z,xl,yl,rz_xyz_bis,z_color,auto_contrast);
            % 'gauss_rendering.m' builds the matrix-image (Gaussian
            % rendering)

        else     
            rz_xyz_bis = 3; % Last operation run = 'Small ROI(x,y,z)'
            z = simpler_output(:,3); % axial positions (z)
            axes(handles.axes1); 
            xl_0 = xlim; % lateral limits (r) are extracted from the scatter plot
            g = msgbox(['Select from the (x,z) scatter plot the ROI to be',...
                ' renderized (w/zoom tool). Then, click "ok" button',...
                ' from this message box']);
            waitfor(g);
            xl_1 = xlim; % lateral limits (x) are extracted from the scatter plot after zoom in
            yl_1 = ylim; % axial limits (z) are extracted from the scatter plot after zoom in
            
            axes(handles.axes3);
            yl_0 = xlim; % lateral limits (y) are extracted from the scatter plot
            g1 = msgbox(['Select from the (y,z) scatter plot the',...
                ' ROI to be renderized (w/zoom tool). ',...
                'Then, click "ok" button from this message box']);
            waitfor(g1);
            xl_2 = xlim; % lateral limits (y) are extracted from the scatter plot after zoom in
            yl_2 = ylim; % axial limits (z) are extracted from the scatter plot after zoom in
      
            % Gaussian rendering (x,y positions)
            gauss_rendering(Magnif,sigma_width_nm(1),sigma_width_nm(1),...
                camera_px,x1,y1,xl_0,yl_0,rz_xyz_bis+2,0,auto_contrast);            
            % Gaussian rendering (x,z positions)
            gauss_rendering(Magnif,sigma_width_nm(1),sigma_width_nm(2),...
                camera_px,x1,z,xl_1,yl_1,rz_xyz_bis,z_color,auto_contrast);
            % Gaussian rendering (y,z positions)
            gauss_rendering(Magnif,sigma_width_nm(1),sigma_width_nm(2),...
                camera_px,y1,z,xl_2,yl_2,(rz_xyz_bis+1),z_color,auto_contrast);
        end

