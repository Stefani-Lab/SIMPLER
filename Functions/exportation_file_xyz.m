function exportation_file_xyz(x1,y1,z,frame,photons,hdf5_filename,handles)

list_xyz = [frame x1 y1 z photons];
list_xz = [frame x1 z photons];
list_yz = [frame y1 z photons];

% Save list as .csv
csv_out_folder = pwd;       % the folder to save the output (converted) file.
out_file_extn = '.csv';     % the extension for the output file.
out_file_delimiter = ',';   % the data-delimiter for the output file
add_headers = true;         % add a headers line to the output file; 

       
        [~, hdf5_name, ~] = fileparts(hdf5_filename); % get the name of the
                                    % input file without the file-extension
        csv_out_fname = fullfile(csv_out_folder,...
            [hdf5_name, '_SIMPLER_xyz',out_file_extn]); % create output name
        if add_headers
            fid = fopen(csv_out_fname, 'wt'); % create and open the output
                                              % file for writing
            fprintf(fid, '%s,%s,%s,%s,%s\n',...
                '"frame"','"x [nm]"','"y [nm]"', '"z [nm]"', '"intensity [photon]"');  
            % write the header line first
            fclose(fid); % close the file so that dlmwrite can use it in the next step
        end
        dlmwrite(csv_out_fname, list_xyz,...
            'delimiter', out_file_delimiter,...
            'precision', '%10.4f', '-append'); % append the rest of the data to the file
        
        
        if sum(get(handles.axes1,'visible'))== 221 % This line checks if there
            % is a scatter plot displayed. If it does, the last operation
            % run is 'Small ROI(x,y,z)' and then two extra .csv files are
            % exported: (x,z) and (y,z) (they represent side views).
            
            [~, hdf5_name, ~] = fileparts(hdf5_filename); % get the name of 
                            % the input file without the file-extension
            csv_out_fname = fullfile(csv_out_folder,...
            [hdf5_name, '_SIMPLER_xz',out_file_extn]); % create output name
            if add_headers
                fid = fopen(csv_out_fname, 'wt'); % create and open the output file for writing
                fprintf(fid, '%s,%s,%s,%s\n', '"frame"','"x [nm]"','"y [nm]"', '"intensity [photon]"');
                % write the header line first
                fclose(fid); % close the file so that dlmwrite can use it in the next step
            end
            dlmwrite(csv_out_fname, list_xz, 'delimiter',...
                out_file_delimiter,...
                'precision', '%10.4f', '-append'); % append the rest of the data to the file

            [~, hdf5_name, ~] = fileparts(hdf5_filename); % get the name of 
                            % the input file without the file-extension
            csv_out_fname = fullfile(csv_out_folder,...
                [hdf5_name, '_SIMPLER_yz',out_file_extn]); % create output name
            if add_headers
                fid = fopen(csv_out_fname, 'wt'); % create and open the output file for writing
                fprintf(fid, '%s,%s,%s,%s\n', '"frame"','"x [nm]"','"y [nm]"', '"intensity [photon]"');
                % write the header line first
                fclose(fid); % close the file so that dlmwrite can use it in the next step
            end
            dlmwrite(csv_out_fname, list_yz,...
                'delimiter', out_file_delimiter,...
                'precision', '%10.4f', '-append'); % append the rest of the data to the file
        end
end  
        
      
       
  