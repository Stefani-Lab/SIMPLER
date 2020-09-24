function exportation_file_rz(x1,y1,proyec_r,z,frame,photons,hdf5_filename,handles)

list_xyz = [frame x1 y1 z photons];
list_rz = [frame proyec_r z photons];

% Save list as .csv
out_file_extn = '.csv';     % the extension for the output file.
out_file_delimiter = ',';   % the data-delimiter for the output file
add_headers = true;         % add a headers line to the output file; 


if sum(get(handles.axes5,'visible')) == 315 % This line checks if
                                         % the 'r vs. z' plot
                                         % corresponds to
                                         % 'Combined
                                         % structures' (from
                                         % the
                                         % 'Incidence angle and
                                         % alpha adjustment
                                         % operation', or to
                                         % 'Small ROI (r,z)').
                                         % If the boolean is
                                         % True, then the case
                                         % is 'Small ROI
                                         % (r,z)' and a .csv
                                         % file including (x,y,z)
                                         % data is also
                                         % exported
                                         % (additionally to the
                                         % (r,z).csv file)

    [~, hdf5_name, ~] = fileparts(hdf5_filename); % get the name of
                   % the input file without the file-extension
    [file,path]= uiputfile([hdf5_name,'_SIMPLER_xyz.csv']);
    csv_out_folder = path;
    csv_out_fname = fullfile(csv_out_folder,...
    file); % create output name                                     
                                                                                  
    
    if add_headers
        fid = fopen(csv_out_fname, 'wt'); % create and open the output
                                          % file for writing
        fprintf(fid, '%s,%s,%s,%s,%s\n','"frame"', '"x [nm]"',...
            '"y [nm]"', '"z [nm]"','"intensity [photon]"');  
        % write the header line first
        fclose(fid); % close the file so that dlmwrite can use it 
                     % in the next step
    end

    dlmwrite(csv_out_fname, list_xyz, 'delimiter',...
    out_file_delimiter, 'precision', '%10.4f', '-append'); 
    % append the rest of the data to the file
    
    
    [~, hdf5_name, ~] = fileparts(hdf5_filename); % get the name of
                   % the input file without the file-extension
    [file,path]= uiputfile([hdf5_name,'_SIMPLER_rz.csv']);
    csv_out_folder = path;
    csv_out_fname = fullfile(csv_out_folder,...
    file); % create output name     

    if add_headers
        fid = fopen(csv_out_fname, 'wt'); % create and open the output 
                                          % file for writing
        fprintf(fid, '%s,%s,%s,%s\n', '"frame"', '"x [nm]"',...
            '"y [nm]"', '"intensity [photon]"');  % write the 
                                                    % header line first
        fclose(fid); % close the file so that dlmwrite can use it in the next step
    end
    dlmwrite(csv_out_fname, list_rz, 'delimiter',...
    out_file_delimiter, 'precision', '%10.4f', '-append'); % append 
                                % the rest of the data to the file


else % If the last operation run is 'Incidence angle and alpha adjustment':

        [~, hdf5_name, ~] = fileparts(hdf5_filename); % get the name of
                   % the input file without the file-extension
    [file,path]= uiputfile([hdf5_name,'_SIMPLER_rz.csv']);
    csv_out_folder = path;
        csv_out_fname = fullfile(csv_out_folder,...
    file); % create output name     
    if add_headers
        fid = fopen(csv_out_fname, 'wt'); % create and open the output 
                                          % file for writing
        fprintf(fid, '%s,%s\n', '"x [nm]"', '"y [nm]"');  % write the 
                                                    % header line first
        fclose(fid); % close the file so that dlmwrite can use it in the next step
    end
    dlmwrite(csv_out_fname, list_rz(:,2:3), 'delimiter',...
    out_file_delimiter, 'precision', '%10.4f', '-append'); % append 
                                % the rest of the data to the file
end
end
        
      
       
  