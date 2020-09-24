function varargout = run_SIMPLER(varargin)
% run_SIMPLER is a graphical user interface that runs in Matlab. This app 
% allows users of SIMPLER to perform all necessary operations to decode 
% the axial positions of single molecules directly from 2D-SMLM-TIRF data.
% The software also includes modules to perform the following operations: 
%
% - Determination of N0 from 2D-SMLM-TIRF data of emitters bound 
% or adsorbed to the coverslip 
%
% - Calculation of the excitation intensity profile from 2D-SMLM-TIRF data
% of molecules spread all over the field of view; the list must include
% the background or offset information for each emitter. 
%
% - Adjustment of the calibration parameters using the SIMPLER 3D
% reconstructions of standard structures of well-defined,
% known geometry as feedback. 
%
%For further details, please refer to the Software Documentation, where we
% include a detailed explanation of how to load files, set the different
% parameters, calculate z-coordinates of single molecules, and run the
% different operations available.


% Begin initialization code - DO NOT EDIT
gui_Singleton = 1;
gui_State = struct('gui_Name',       mfilename, ...
                   'gui_Singleton',  gui_Singleton, ...
                   'gui_OpeningFcn', @run_SIMPLER_OpeningFcn, ...
                   'gui_OutputFcn',  @run_SIMPLER_OutputFcn, ...
                   'gui_LayoutFcn',  [] , ...
                   'gui_Callback',   []);
if nargin && ischar(varargin{1})
    gui_State.gui_Callback = str2func(varargin{1});
end

if nargout
    [varargout{1:nargout}] = gui_mainfcn(gui_State, varargin{:});
else
    gui_mainfcn(gui_State, varargin{:});
end




% End initialization code - DO NOT EDIT


% --- Executes just before run_SIMPLER is made visible.
function run_SIMPLER_OpeningFcn(hObject, eventdata, handles, varargin)
% This function has no output args, see OutputFcn.
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% varargin   command line arguments to run_SIMPLER (see VARARGIN)

% First, plots are cleared and previous SIMPLER otuputs are deleted.
clear simpler_output output_ft_input
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
handles.output = hObject;
set(handles.axes1,'visible', 'off'); 
set(handles.axes3,'visible', 'off');
set(handles.open_fig_1_tag,'visible', 'off'); 
set(handles.open_fig_3_tag,'visible', 'off'); 
set(handles.ft_panel_tag,'visible', 'off'); 
set(handles.export_csv_1_tag,'visible', 'off');
set(handles.export_csv_3_tag,'visible', 'off');
set(handles.export_csv_whole_1_tag,'visible', 'off');
set(handles.export_csv_whole_3_tag,'visible', 'off');
set(handles.panel_ft_input_tag,'visible', 'off');
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





addpath('Functions')
addpath('Example data')





% Update handles structure
guidata(hObject, handles);

% UIWAIT makes run_SIMPLER wait for user response (see UIRESUME)
% uiwait(handles.figure1);


% --- Outputs from this function are returned to the command line.
function varargout = run_SIMPLER_OutputFcn(hObject, eventdata, handles) 
% varargout  cell array for returning output args (see VARARGOUT);
% hObject    handle to figure
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Get default command line output from handles structure
varargout{1} = handles.output;



function N0_tag_Callback(hObject, eventdata, handles)
% hObject    handle to N0_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% Hints: get(hObject,'String') returns contents of N0_tag as text
%        str2double(get(hObject,'String')) returns contents of N0_tag as a double


% --- Executes during object creation, after setting all properties.
function N0_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N0_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called
% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_tag_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_tag as text
%        str2double(get(hObject,'String')) returns contents of alpha_tag as a double


% --- Executes during object creation, after setting all properties.
function alpha_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_tag (see GCBO)

% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle_tag_Callback(hObject, eventdata, handles)
% hObject    handle to angle_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_tag as text
%        str2double(get(hObject,'String')) returns contents of angle_tag as a double


% --- Executes during object creation, after setting all properties.
function angle_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in pushbutton1.
function pushbutton1_Callback(hObject, eventdata, handles)
% This funcion is executed when the 'run_SIMPLER' button is pushed. First,
% it reads the calibration input parameters. Then, depending on the
% selected operation, it runs other functions (SIMPLER.m -which contains
% different commands paths depending on the required operation-
% or 'fine_tuning.m'. Also, exportation functions may be executed.

N0 = str2double(get(handles.N0_tag,'string')); % N0
alpha = str2double(get(handles.alpha_tag,'string')); % alpha
angle = str2double(get(handles.angle_tag,'string')); % Incidence angle
filename_wformat = get(handles.filename_tag,'string'); % Filename (including format)
yn_filename_wformat = exist(filename_wformat,'file'); % This line checks if 
                                                      % the file containing the
                                                      % SMLM data exists.
if yn_filename_wformat == 0 % If file is not found, an error message is shown.
    w = msgbox(['Filename not found. Please check. If it is not located',...
        ' in the current directory, include the directory path',...
        ' in the filename: "C:\Dir\filename.hdf5"']);
    waitfor(w);
end
lambda_exc = str2double(get(handles.exc_lambda_tag,'string')); % Excitation 
                                                               % wavelength.
lambda_em = str2double(get(handles.em_wavelength_tag,'string'));% Emission 
                                                               % wavelength.
nI = str2double(get(handles.ni_tag,'string')); % Refractive index of the
                                               % incidence medium.
nS = str2double(get(handles.ns_tag,'string')); % Refractive index of the
                                               % sample medium.
NA_pop = get(handles.pop_up_NA_tag,'Value');
NA = 1.42*(NA_pop==1)+1.45*(NA_pop==2)+1.49*(NA_pop==3); % Numerical aperture
Illum = get(handles.corr_nonflat_tag,'Value'); % This line checks if the 'Correction
                                               % (non-flat illumination')
                                               % box is checked.
max_dist = str2double(get(handles.max_dist_tag,'string')); % Maximum distance
                                                         % tolerable to
                                                         % link
                                                         % localizations
                                                         % from consecutive
                                                         % frames.
scatter_plot = get(handles.scatter_full_tag,'Value'); % 'scatter_plot' is 1
                                                      % if the user has
                                                      % checked the
                                                      % 'Scatter plot'
                                                      % visualization
                                                      % option, or 0
                                                      % otherwise.
camera_px = str2double(get(handles.px_size_tag,'string')); % Camera pixel size (nm)
rz_xyz = get(handles.pop_up_tag,'Value'); % This variable contains the information
                                          % of the selected operation.
file_format = get(handles.fileformat_tag,'Value'); % File format.

if Illum ==1 % If the correction for uneven illumination option is selected,
             % the function checks if the filename of the excitation
             % profile file exists. If it doesn't, an error message is
             % shown.
    calib_file_name = get(handles.excprof_filename_tag,'string');
    yn_calib_file_name = exist(calib_file_name,'file');
    if yn_calib_file_name == 0 && sum(get(handles.excprof_filename_tag,...
        'enable'))== 221
        ww = msgbox(['Excitation profile filename not found.',...
            ' Please check. If it is not located in the current',...
            ' directory, include the directory path',...
            ' in the filename: "C:\Dir\filename.hdf5"']);
    waitfor(ww);
    end
else % If no correction is needed, the variables related to the excitation
     % profile file are set to 1.
    calib_file_name = 1;
    yn_calib_file_name = 1;
end
% Before running other functions ('SIMPLER.m' or 'fine_tuning.m'), the
% function checks if there is a mismatch between the file format introduced
% and the one deduced from the filename.
if (file_format == 1 && (sum(filename_wformat(end-3:end)=='hdf5')==4))...
        || ((file_format == 2 || file_format == 3)...
        && (sum(filename_wformat(end-3:end)=='.csv')==4))
    
    % If theres is no mismatch, then the code continues this way:
    
    if file_format == 3 && rz_xyz < 6 % If the file format is 'Other' (csv)
                                      % and the operation is not 'Incidence
                                      % angle and alpha adjustment',
                                      % a message box will remind the user
                                      % how the columns from the file must
                                      % be sorted.
        r = msgbox(['You have chosen "Other" csv file format. Please ',...
            'make sure that the order of the columns in your list is',...
            ' as follows: 1) frame, 2) x (nm), 3) y (nm), 4) intensity',...
            ' (photons) and 5) offset or background (photons; only',...
            ' needed if using the correction for non-flat illumination.',...
            ' Also, remove headers. Once verified, click Ok to continue.']);
        waitfor(r);
    end
    if rz_xyz == 5 % If the selected operation is N0 Calibration, a message
                   % box describing the output will be shown.
        u = msgbox(['You have chosen "N0 calibration". After running',...
            ' SIMPLER, you will get a histogram of the number of ',...
            'photons together with the results of a 1-peak gaussian fit.',...
            ' Click Ok to confirm and run calibration.']);
        waitfor(u);
    end
    if rz_xyz == 6 && yn_filename_wformat>0 % If the selected operation is 
                                            % 'Incidence angle & alpha adjustment',
                                            % a message box describing 
                                            % the output will be shown.
       uu = msgbox(['You have chosen "Incidence angle & alpha ',...
           'adjustment". After running SIMPLER, you will get a scatter',...
           ' plot of all the known structures plotted together. Click ',...
           'Ok to confirm and run this operation.']);
        waitfor(uu);
        % When the 'Incidence angle & alpha adjustment' operation is
        % chosen, the 'fine_tuning.m' function is run, and the results form
        % that execution are saved in the 'output_ft_input' array.
         output_ft_input = fine_tuning(filename_wformat,handles);
         assignin('base','output_ft_input', output_ft_input) % The output of 
                                                             % the
                                                             % 'fine_tuning.m'
                                                             % function is
                                                             % kept in the
                                                             % workspace.
                                                             
        if evalin('base', 'exist(''simpler_output'',''var'')')
            simpler_output = evalin('base','simpler_output'); % Here, the script 
                                                          %checks if the
                                                          %'SIMPLER.m'
                                                          %function has
                                                          %been executed
                                                          %previously. If it
                                                          %has, then the
                                                          %fifth column
                                                          %from the
                                                          %'simpler_output'
                                                          %array is
                                                          %overwritten with
                                                          %0.5 values. This
                                                          %information is
                                                          %used if the
                                                          %user decides to
                                                          %export the data
                                                          %of the scatter
                                                          %plot of the
                                                          %combined
                                                          %structures. If
                                                          %the array
                                                          %doesn't exist,
                                                          %then
                                                          %'simpler_output'
                                                          %is created and
                                                          %the fifth column
                                                          %is set to 0.5. 
            simpler_output(:,5) = 0.5*ones(size(simpler_output,1),1);
        else
            simpler_output(:,5) = 0.5;
        end
    assignin('base','simpler_output',simpler_output) % 'simpler_output' is 
                                                       % kept in the
                                                       % workspace.    
    end
     % If the chosen operation is Small ROI, Large ROI or N0 Calibration,
     % the function 'SIMPLER.m' is executed. The array 'simpler_output'
     % contains the results obtained after running SIMPLER.
     if rz_xyz >1 && rz_xyz<6 && yn_filename_wformat>0 && yn_calib_file_name>0 
      simpler_output = SIMPLER(lambda_exc,NA,lambda_em,angle,nS,nI,...
          alpha,N0,filename_wformat,Illum,max_dist,scatter_plot,camera_px,...
          rz_xyz,file_format,hObject,handles,calib_file_name);       
      
      assignin('base','simpler_output',simpler_output) % 'simpler_output' is 
                                                       % kept in the
                                                       % workspace.
      set(handles.export_csv_whole_3_tag,'enable','on');
         if (rz_xyz == 4) || (rz_xyz == 5) % If the selected operation is 
                                           % Large ROI or N0 Calibration,
                                           % the results are exported to
                                           % .csv files.
            export_function(0,simpler_output,handles);
           if rz_xyz == 4
                ff = msgbox('Done! Large ROI data has been exported to .csv');
           else
                ff = msgbox(['Done! (x,y,z) data and a list with N0 values',...
                    ' have been exported to .csv']);
           end
         end
        if sum(simpler_output(:,5)) > 0 && scatter_plot == 1 % If the Operation
                                                             % selected are
                                                             % visualized with
                                                             % scatter plots,
                                                             % the button that
                                                             % allows to
                                                             % perform a
                                                             % Gaussian
                                                             % rendering is
                                                             % enabled.
                                                             % Otherwise, it is
                                                             % disabled.
            set(handles.render_tag,'enable','on');
        end
        if scatter_plot == 0 || sum(get(handles.axes7,'visible')) == 221 % In any other case,
                                                        % or if the
                                                        % selected
                                                        % operation was
                                                        % "Large ROI", the
                                                        % Gaussian
                                                        % rendering option
                                                        % is disabled.
            set(handles.render_tag,'enable','off');
        end     
    elseif rz_xyz==1 % If no Operation is chosen, an error message is shown.
        w = msgbox(['Please choose the operation to ',...
            'be performed in the "Select Operation" menu (bottom)']);
    end
else
    q = msgbox(['Mismatch between file extension from',...
            ' filename field and from file format. Please check and retry.']);
     waitfor(q);
end




function filename_tag_Callback(hObject, eventdata, handles)
% hObject    handle to filename_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of filename_tag as text
%        str2double(get(hObject,'String')) returns contents of filename_tag as a double


% --- Executes during object creation, after setting all properties.
function filename_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to filename_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function em_wavelength_tag_Callback(hObject, eventdata, handles)
% hObject    handle to em_wavelength_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of em_wavelength_tag as text
%        str2double(get(hObject,'String')) returns contents of em_wavelength_tag as a double


% --- Executes during object creation, after setting all properties.
function em_wavelength_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to em_wavelength_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function exc_lambda_tag_Callback(hObject, eventdata, handles)
% hObject    handle to exc_lambda_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exc_lambda_tag as text
%        str2double(get(hObject,'String')) returns contents of exc_lambda_tag as a double


% --- Executes during object creation, after setting all properties.
function exc_lambda_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exc_lambda_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ni_tag_Callback(hObject, eventdata, handles)
% hObject    handle to ni_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ni_tag as text
%        str2double(get(hObject,'String')) returns contents of ni_tag as a double


% --- Executes during object creation, after setting all properties.
function ni_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ni_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ns_tag_Callback(hObject, eventdata, handles)
% hObject    handle to ns_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ns_tag as text
%        str2double(get(hObject,'String')) returns contents of ns_tag as a double


% --- Executes during object creation, after setting all properties.
function ns_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ns_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function NA_tag_Callback(hObject, eventdata, handles)
% hObject    handle to NA_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NA_tag as text
%        str2double(get(hObject,'String')) returns contents of NA_tag as a double


% --- Executes during object creation, after setting all properties.
function NA_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in corr_nonflat_tag.
function corr_nonflat_tag_Callback(hObject, eventdata, handles)
% This function acts like a 'switch', enabling or disabling the editable
% textbox where the excitation profile filename should be entered (it is
% enabled only if the 'Correction (non-flat illumination)' box is checked).
yn = get(hObject,'Value');
if yn == 1
    set(handles.excprof_filename_tag,'enable','on');
else
    set(handles.excprof_filename_tag,'enable','off');
end

    


% Hint: get(hObject,'Value') returns toggle state of corr_nonflat_tag



function edit23_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_tag as text
%        str2double(get(hObject,'String')) returns contents of alpha_tag as a double


% --- Executes during object creation, after setting all properties.
function edit23_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit24_Callback(hObject, eventdata, handles)
% hObject    handle to NA_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of NA_tag as text
%        str2double(get(hObject,'String')) returns contents of NA_tag as a double


% --- Executes during object creation, after setting all properties.
function edit24_CreateFcn(hObject, eventdata, handles)
% hObject    handle to NA_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit20_Callback(hObject, eventdata, handles)
% hObject    handle to exc_lambda_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of exc_lambda_tag as text
%        str2double(get(hObject,'String')) returns contents of exc_lambda_tag as a double


% --- Executes during object creation, after setting all properties.
function edit20_CreateFcn(hObject, eventdata, handles)
% hObject    handle to exc_lambda_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_dist_tag_Callback(hObject, eventdata, handles)
% hObject    handle to max_dist_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_dist_tag as text
%        str2double(get(hObject,'String')) returns contents of max_dist_tag as a double


% --- Executes during object creation, after setting all properties.
function max_dist_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_dist_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function px_size_tag_Callback(hObject, eventdata, handles)
% hObject    handle to px_size_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of px_size_tag as text
%        str2double(get(hObject,'String')) returns contents of px_size_tag as a double


% --- Executes during object creation, after setting all properties.
function px_size_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to px_size_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in scatter_plot_tirf_tag.
function scatter_plot_tirf_tag_Callback(hObject, eventdata, handles)
% hObject    handle to scatter_plot_tirf_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scatter_plot_tirf_tag


% --- Executes on button press in scatter_plot_tag.
function scatter_plot_tag_Callback(hObject, eventdata, handles)
% hObject    handle to scatter_plot_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scatter_plot_tag


% --- Executes on button press in gaussian_tag.
function gaussian_tag_Callback(hObject, eventdata, handles)
% hObject    handle to gaussian_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gaussian_tag



function gauss_width_tag_Callback(hObject, eventdata, handles)
% hObject    handle to gauss_width_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gauss_width_tag as text
%        str2double(get(hObject,'String')) returns contents of gauss_width_tag as a double


% --- Executes during object creation, after setting all properties.
function gauss_width_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gauss_width_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function magnification_tag_Callback(hObject, eventdata, handles)
% hObject    handle to magnification_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of magnification_tag as text
%        str2double(get(hObject,'String')) returns contents of magnification_tag as a double


% --- Executes during object creation, after setting all properties.
function magnification_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to magnification_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in scatter_full_tag.
function scatter_full_tag_Callback(hObject, eventdata, handles)
% hObject    handle to scatter_full_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scatter_full_tag


% --- Executes on button press in scatter_tag.
function scatter_tag_Callback(hObject, eventdata, handles)
% hObject    handle to scatter_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of scatter_tag


% --- Executes on button press in gauss_rend_tag.
function gauss_rend_tag_Callback(hObject, eventdata, handles)
% hObject    handle to gauss_rend_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of gauss_rend_tag


% --- Executes on button press in export_tag.
function export_tag_Callback(hObject, eventdata, handles)

% hObject    handle to export_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of export_tag


% --- Executes on selection change in pop_up_tag.
function pop_up_tag_Callback(hObject, eventdata, handles)
% This function changes the visibility of certain buttons, as well as the
% 'enable-disable' state of some textboxes, depending on which operation
% the user selects. For example, if the user's choice is 'N0 Calibration',
% then the N0 editable textbox gets disabled, because this operation
% assumes that the fluorescent molecules are located at z = 0, and
% estimates a N0 value from the distribution of #photons through the whole
% field of view.
if get(handles.pop_up_tag,'Value') == 7
   set(handles.calc_exc_prof_tag,'visible','on');
   set(handles.corr_nonflat_tag,'enable','off');
   set(handles.excprof_filename_tag,'enable','off');
   set(handles.pushbutton1,'enable','off');
else
    set(handles.pushbutton1,'enable','on');
   set(handles.calc_exc_prof_tag,'visible','off');   
    if get(handles.pop_up_tag,'Value') == 6
        set(handles.corr_nonflat_tag,'enable','off');
        set(handles.max_dist_tag,'enable','off');
        set(handles.excprof_filename_tag,'enable','off');
    else
         set(handles.corr_nonflat_tag,'enable','on');
        set(handles.max_dist_tag,'enable','on'); 
        if get(handles.corr_nonflat_tag,'Value')==1
        set(handles.excprof_filename_tag,'enable','on');
        end
    end
    if get(handles.pop_up_tag,'Value') == 5  
        set(handles.N0_tag,'enable','off');
    else
        set(handles.N0_tag,'enable','on');
    end
end


% --- Executes during object creation, after setting all properties.
function pop_up_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_up_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in listbox1.
function listbox1_Callback(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns listbox1 contents as cell array
%        contents{get(hObject,'Value')} returns selected item from listbox1


% --- Executes during object creation, after setting all properties.
function listbox1_CreateFcn(hObject, eventdata, handles)
% hObject    handle to listbox1 (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: listbox controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gauss_width_axial_tag_Callback(hObject, eventdata, handles)
% hObject    handle to gauss_width_axial_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gauss_width_axial_tag as text
%        str2double(get(hObject,'String')) returns contents of gauss_width_axial_tag as a double


% --- Executes during object creation, after setting all properties.
function gauss_width_axial_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gauss_width_axial_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in pop_up_NA_tag.
function pop_up_NA_tag_Callback(hObject, eventdata, handles)
% hObject    handle to pop_up_NA_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns pop_up_NA_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from pop_up_NA_tag


% --- Executes during object creation, after setting all properties.
function pop_up_NA_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to pop_up_NA_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in fileformat_tag.
function fileformat_tag_Callback(hObject, eventdata, handles)
% hObject    handle to fileformat_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns fileformat_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from fileformat_tag


% --- Executes during object creation, after setting all properties.
function fileformat_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to fileformat_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in zcolor_tag.
function zcolor_tag_Callback(hObject, eventdata, handles)
% hObject    handle to zcolor_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of zcolor_tag


% --- Executes on button press in showcalib_tag.
function showcalib_tag_Callback(hObject, eventdata, handles)

% hObject    handle to showcalib_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of showcalib_tag



function excprof_filename_tag_Callback(hObject, eventdata, handles)
% hObject    handle to excprof_filename_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of excprof_filename_tag as text
%        str2double(get(hObject,'String')) returns contents of excprof_filename_tag as a double


% --- Executes during object creation, after setting all properties.
function excprof_filename_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to excprof_filename_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in export_csv_whole_3_tag.
function export_csv_whole_3_tag_Callback(hObject, eventdata, handles)
% This function is executed when the 'Export CSV (whole dataset)' button is
% pushed (for y vs. z). It calls 'export_function.m', where the exportation
% commands can be found.
simpler_output = evalin('base','simpler_output');
export_function(0,simpler_output,handles);
% hObject    handle to export_csv_whole_3_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in open_fig_3_tag.
function open_fig_3_tag_Callback(hObject, eventdata, handles)
% This function open the y vs. z (or N0 histogram) figure in a separate window.
simpler_output = evalin('base','simpler_output');
  if sum(simpler_output(:,5)) == size(simpler_output,1) % If the current plot
                                                        % at the bottom of
                                                        % the
                                                        % 'Visualization'
                                                        % panel is 'y vs
                                                        % z', then...
        y1 = simpler_output(:,2); % y values are obtained from simpler_output
        z = simpler_output(:,3); % z values are obtained from simpler_output
        figure, scatter(y1,z,10) % y vs. z is plotted in a new figure
        xlabel('y (nm)');
        ylabel('z (nm)');
        daspect([1 1 1])
        pbaspect([1 1 1])
  elseif sum(simpler_output(:,5)) == 0  % If the current plot at the bottom of
                                                        % the 'Visualization'
                                                        % panel is the N0
                                                        % histogram,
                                                        % then...
      x1 = simpler_output(:,1);
      y1 = simpler_output(:,2);
      photons1 = simpler_output(:,4); % #photons for each SM 
      
      xl = [min(x1) max(x1)]; yl = [min(y1) max(y1)]; % x and y limits
      c = 1:length(x1);
      figure
      if exist('histogram') >0
          histogram(photons1(c)); % plot N0 histogram in new figure
          h = findobj(gcf, 'Type', 'histogram');
          b = h.Values;
          L = h.BinLimits;
          a = L(1):(L(2)-L(1))/(length(b)-1):L(2);
      else
          hist(photons1(c)); 
          [b a] = hist(photons1(c));
          L = [a(1) a(end)];
      end
      xlimit =L;
      xl_fit = xlimit;   
      f = fit(a',b','gauss1'); % Perform 1-peak Gaussian fit on the N0 distribution
      
      MyCoeffs = coeffvalues(f);
      a1 = MyCoeffs(1); b1 = MyCoeffs(2); c1 = MyCoeffs(3);
      sigma = c1/sqrt(2);
      x_plot = xl_fit(1):xl_fit(2);
      y_plot = f(x_plot);
      hold on, plot(x_plot,y_plot,'LineWidth',3);
      xlabel('N_0 (photons)');
      ylabel('Counts');
      title({'N_0 histogram',['N_0 = ',num2str(round(b1)), ' photons; \sigma_N_0 = ',...
        num2str(round(sigma)),' photons'],['Area = ',num2str(...
        round(abs(1E-6*(xl(2)-xl(1))*(yl(2)-yl(1))))),' \mum^2']});    
  end


% --- Executes on button press in open_fig_1_tag.
function open_fig_1_tag_Callback(hObject, eventdata, handles)
simpler_output = evalin('base','simpler_output');
  if sum(simpler_output(:,5)) > size(simpler_output,1)
        proyec_r = simpler_output(:,3);
        z = simpler_output(:,4);
        figure, scatter(proyec_r,z,10)
        xlabel('r (nm)');
        ylabel('z (nm)');
        daspect([1 1 1])
        pbaspect([1 1 1])
  elseif sum(simpler_output(:,5)) == size(simpler_output,1)
        x1 = simpler_output(:,1);
        z = simpler_output(:,3);
        figure, scatter(x1,z,10)
        xlabel('x (nm)');
        ylabel('z (nm)');
        daspect([1 1 1])
        pbaspect([1 1 1])
  end 

% hObject    handle to open_fig_1_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on button press in render_tag.
function render_tag_Callback(hObject, eventdata, handles)
% This function performs the Gaussian rendering of the data displayed in
% the scatter plots. 
check_ax1 = get(handles.axes1,'visible');
check_ax8 = get(handles.axes8,'visible');
%First, it checks if the last operation run was 'Small ROI' or 'Incidence
%angle and alpha adjustment'. If check_ax1 is visible (check_ax1 = 'on'),
%it means that the last performed operation was 'Small ROI', and axref is
%set to 1. Additionally, the 'simpler_output' array is loaded from
%workspace. From that array, the [lateral,axial] positions are read and
%used to perform the Gaussian rendering.
if evalin('base', 'exist(''simpler_output'',''var'')') && check_ax1(2) == 'n'
    simpler_output = evalin('base','simpler_output');
    axref = 1;
elseif evalin('base', 'exist(''simpler_output'',''var'')') && check_ax8(2) == 'n'
    simpler_output = evalin('base','simpler_output');
    roi_index_proyec_r = evalin('base','roi_index_proyec_r');
    axref = roi_index_proyec_r; 
    % In this case, the last operation run was 'Large ROI'.
    
else % In this case, the last run operation corresponds to 'Incidence angle
    % and alpha adjustment'. The results from that operation are loaded as
    % 'output_ft_intput' and axref takes a value of 5. 
    % The relative-lateral positions are obtained from 'output_ft_input',
    % as well as the calibration parameters used in the SIMPLER
    % measurements, except for 'angle' and 'alpha', which are taken from
    % the 'Tune' panel. After reading all the parameters, dF and alphaF are
    % obtained through 'getParameters_SIMPLER, and with those values in hand
    % the axial positions are calculated. A new array called
    % 'simpler_output' is generated, which contains the [axial,lateral]
    % positions information. This array is used within this function but it
    % is not exported to the workspace (thus, it doesn't replace the last 
    % generated and exported to workspace 'simpler_output' matrix). 
    
    output_ft_input = evalin('base','output_ft_input');
    angle_ft =  str2double((get(handles.angle_ft_tag,'string')));
    alpha_ft =  str2double((get(handles.alpha_ft_tag,'string')));
    photons = output_ft_input(6:end,1);
    photons_median = output_ft_input(6:end,3);
    lambda_exc_input = output_ft_input(1,1);
    NA_input = output_ft_input(2,1);
    lambda_em_input = output_ft_input(3,1);
    nS_input = output_ft_input(4,1);
    nI_input = output_ft_input(5,1);
    lateral = output_ft_input(6:end,2);
    lateral_median = output_ft_input(6:end,4);
    N0_input = str2double((get(handles.N0_input_ft_tag,'string')));
    [dF_ft, alphaF_ft] = getParameters_SIMPLER(lambda_exc_input,NA_input,lambda_em_input,...
        angle_ft,nS_input,nI_input,alpha_ft);
    axial_ft = (log(alphaF_ft*N0_input)-log(photons-(1-alphaF_ft)*N0_input))/(1/dF_ft);
    axial_ft_median = (log(alphaF_ft*N0_input)-log(photons_median-(1-alphaF_ft)*N0_input))/(1/dF_ft);
    simpler_output(:,1) = ones(length(photons),1);
    simpler_output(:,2) = ones(length(photons),1);
    simpler_output(:,3) = lateral-lateral_median;
    simpler_output(:,4) = axial_ft-axial_ft_median;
    simpler_output(:,5) = photons;
    axref = 5;
end
% Once the information of the [lateral,axial] positions is obtained, the function
% reads the rendering options set by the user:
auto_contrast = get(handles.autoadjust_contrast_tag,'Value');
Magnif = str2double(get(handles.magnification_tag,'string'));
camera_px = str2double(get(handles.px_size_tag,'string'));
sigma_width_nm_lateral = str2double(get(handles.gauss_width_tag,'string'));
sigma_width_nm_axial = str2double(get(handles.gauss_width_axial_tag,'string'));
sigma_width_nm = [sigma_width_nm_lateral sigma_width_nm_axial];
z_color = get(handles.zcolor_tag,'Value');
gauss_rendering_function(Magnif,sigma_width_nm,camera_px,simpler_output,z_color,...
        handles,axref,auto_contrast); % This function continues the
                                              % Gaussian rendering
                                              % execution





function calib_output_tag_Callback(hObject, eventdata, handles)
% hObject    handle to calib_output_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of calib_output_tag as text
%        str2double(get(hObject,'String')) returns contents of calib_output_tag as a double


% --- Executes during object creation, after setting all properties.
function calib_output_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to calib_output_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in refresh_calib_tag.
function refresh_calib_tag_Callback(hObject, eventdata, handles)
% This function calculates alphaF and dF with the calibration parameters
% currently set. It runs after the 'Update calibration(dF and alphaF)'
% button is pushed. The output dF and alphaF values are shown just below
% this button.
N0 = str2double(get(handles.N0_tag,'string')); % N0
alpha = str2double(get(handles.alpha_tag,'string')); % alpha
angle = str2double(get(handles.angle_tag,'string')); % Incidence angle
lambda_exc = str2double(get(handles.exc_lambda_tag,'string')); % Excitation wavelength
lambda_em = str2double(get(handles.em_wavelength_tag,'string')); % Emission wavelength
nI = str2double(get(handles.ni_tag,'string')); % Refractive index (incidence medium)
nS = str2double(get(handles.ns_tag,'string')); % Refractive index (sample)
NA_pop = get(handles.pop_up_NA_tag,'Value'); 
NA = 1.42*(NA_pop==1)+1.45*(NA_pop==2)+1.49*(NA_pop==3); % Numerical aperture
[dF, alphaF] = getParameters_SIMPLER(lambda_exc,NA,lambda_em,angle,nS,nI,alpha);
calib_output_str = ['dF = ',sprintf('%1.1f',dF),'; alphaF = ', sprintf('%1.2f',alphaF)];
 set(handles.calib_output_tag,'string',calib_output_str); 



function tita_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to tita_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of tita_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of tita_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function tita_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to tita_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of alpha_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function alpha_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function edit44_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of alpha_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function edit44_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to angle_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of angle_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function angle_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in export_csv_1_tag.
function export_csv_1_tag_Callback(hObject, eventdata, handles)
% This function executes 'export_function.m' when 'Export CSV (current
% view)' is pushed (for x vs. z plot)
simpler_output = evalin('base','simpler_output');
export_function(1,simpler_output,handles);



% --- Executes on button press in export_csv_3_tag.
function export_csv_3_tag_Callback(hObject, eventdata, handles)
% This function executes 'export_function.m' when 'Export CSV (current
% view)' is pushed (for y vs. z plot)
simpler_output = evalin('base','simpler_output');
export_function(3,simpler_output,handles);



% --- Executes on button press in export_csv_whole_1_tag.
function export_csv_whole_1_tag_Callback(hObject, eventdata, handles)
% This function executes 'export_function.m' when 'Export CSV (whole
% dataset)' is pushed (for x or r vs. z plot). The information of wether
% the current scatter plot represents (x vs. z) (Small ROI), (r vs. z)
% (Small ROI) or (r vs. z) (Combined structures) is deduced from the values
% of the fifth column of 'simpler_output' array.

simpler_output = evalin('base','simpler_output');

export_function(0,simpler_output,handles);
% hObject    handle to export_csv_whole_1_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)



function angle_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to angle_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_input_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of angle_input_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function angle_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function N0_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to N0_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of N0_input_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of N0_input_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function N0_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to N0_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_input_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of alpha_input_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function alpha_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ni_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to ni_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ni_input_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of ni_input_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function ni_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ni_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function emwave_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to emwave_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of emwave_input_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of emwave_input_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function emwave_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to emwave_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function ns_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to ns_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of ns_input_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of ns_input_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function ns_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ns_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in na_input_ft_tag.
function na_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to na_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns na_input_ft_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from na_input_ft_tag


% --- Executes during object creation, after setting all properties.
function na_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to na_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function excwave_input_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to excwave_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of excwave_input_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of excwave_input_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function excwave_input_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to excwave_input_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in update_ft_tag.
function update_ft_tag_Callback(hObject, eventdata, handles)
% This function is executed when the 'Update scatter' button is pushed,
% within the 'Incidence angle and alpha adjustment' operation. 
% It reads the calibration parameters from the 'Parameters used' panel,
% except for the angle and alpha values, which are obtained from the 'Tune'
% panel. The #photons for each SM is taken from the 'otuput_ft_input' array,
% as well as the lateral positions.
% Then, dF and alphaF are calculated and the axial positions are retrieved 
% from the new calibration output parameters. 

output_ft_input = evalin('base','output_ft_input');
angle_ft =  str2double((get(handles.angle_ft_tag,'string')));
alpha_ft =  str2double((get(handles.alpha_ft_tag,'string')));
photons = output_ft_input(6:end,1);
photons_median = output_ft_input(6:end,3);
lambda_exc_input = output_ft_input(1,1);
NA_input = output_ft_input(2,1);
lambda_em_input = output_ft_input(3,1);
nS_input = output_ft_input(4,1);
nI_input = output_ft_input(5,1);
lateral = output_ft_input(6:end,2);
lateral_median = output_ft_input(6:end,4);
N0_input = str2double((get(handles.N0_input_ft_tag,'string')));
[dF_ft, alphaF_ft] = getParameters_SIMPLER(lambda_exc_input,NA_input,lambda_em_input,...
    angle_ft,nS_input,nI_input,alpha_ft);
update_axial(lateral,photons,dF_ft,alphaF_ft,...
    N0_input, handles,photons_median,lateral_median); % This function carries
% out the update of all the plots involved in the 'Incidence angle and
% alpha adjustment' operation. For further details see the comments on
% 'update_axial.m'



function axial_range_tag_Callback(hObject, eventdata, handles)
% hObject    handle to axial_range_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of axial_range_tag as text
%        str2double(get(hObject,'String')) returns contents of axial_range_tag as a double


% --- Executes during object creation, after setting all properties.
function axial_range_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to axial_range_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function max_ax_dist_tag_Callback(hObject, eventdata, handles)
% hObject    handle to max_ax_dist_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of max_ax_dist_tag as text
%        str2double(get(hObject,'String')) returns contents of max_ax_dist_tag as a double


% --- Executes during object creation, after setting all properties.
function max_ax_dist_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to max_ax_dist_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lateral_range_1_tag_Callback(hObject, eventdata, handles)
update_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to lateral_range_1_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lateral_range_1_tag as text
%        str2double(get(hObject,'String')) returns contents of lateral_range_1_tag as a double


% --- Executes during object creation, after setting all properties.
function lateral_range_1_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lateral_range_1_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function lateral_range_3_tag_Callback(hObject, eventdata, handles)
update_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to lateral_range_3_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of lateral_range_3_tag as text
%        str2double(get(hObject,'String')) returns contents of lateral_range_3_tag as a double


% --- Executes during object creation, after setting all properties.
function lateral_range_3_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to lateral_range_3_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in autoadjust_contrast_tag.
function autoadjust_contrast_tag_Callback(hObject, eventdata, handles)
% hObject    handle to autoadjust_contrast_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of autoadjust_contrast_tag



function alpha_step_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_step_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_step_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of alpha_step_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function alpha_step_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_step_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function alpha_range_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to alpha_range_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of alpha_range_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of alpha_range_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function alpha_range_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to alpha_range_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle_step_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to angle_step_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_step_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of angle_step_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function angle_step_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_step_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function angle_range_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to angle_range_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of angle_range_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of angle_range_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function angle_range_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_range_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z2_auto_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to z2_auto_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z2_auto_ft_tag as text
%        str2double(get(hObject,'String')) returns contents of z2_auto_ft_tag as a double


% --- Executes during object creation, after setting all properties.
function z2_auto_ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z2_auto_ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function z1_auto_Ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to z1_auto_Ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z1_auto_Ft_tag as text
%        str2double(get(hObject,'String')) returns contents of z1_auto_Ft_tag as a double


% --- Executes during object creation, after setting all properties.
function z1_auto_Ft_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z1_auto_Ft_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



% --- Executes on button press in run_auto_ft_tag.
function run_auto_ft_tag_Callback(hObject, eventdata, handles)
% This function retrieves the relative-axial positions of the structures
% with known geometry used for the angle and alpha adjustment at different
% [alpha, angle] values. For both parameters, the range and the steps are
% set in the 'Fit for different angles and alpha'. In each condition,
% either a circular or a 2-peaks Gaussian fit on the [r,z] data is
% performed. The radius (circular fit) or the axial position of the 2
% peaks (Gaussian fit) are saved and in the end are plotted as follows:
% for each alpha, either the diameters (2*r) or the differences between the 
% two peaks (dz = z2-z1) are plotted against the angle of incidence. 
angle = str2double(get(handles.angle_ft_tag,'string')); % Central value of the 'angles' set
alpha = str2double(get(handles.alpha_ft_tag,'string')); % Central value of the 'alphas' set
angle_range = str2double(get(handles.angle_range_ft_tag,'string')); % Angle range
steps_angle = str2double(get(handles.angle_step_ft_tag,'string')); % Steps (angle)
alpha_range = str2double(get(handles.alpha_range_ft_tag,'string')); % Alpha range
steps_alpha = str2double(get(handles.alpha_step_ft_tag,'string')); % Steps (alpha)
angles = (angle-0.5*angle_range):steps_angle:(angle+0.5*angle_range); % 'Angle's (vector) 
alphas = (alpha-0.5*alpha_range):steps_alpha:(alpha+0.5*alpha_range); % 'Alpha's (vector) 

% The following loops span the whole range of [angle, alpha] values set by
% the user. The fit results are kept in 'coef2' and 'coef4' variables.
for i = 1:length(angles)
    for k = 1:length(alphas)
        set(handles.angle_ft_tag,'string',num2str(angles(i)));
        set(handles.alpha_ft_tag,'string',num2str(alphas(k)));
        if get(handles.ft_auto_method_tag,'Value') == 2 % For 2-peaks Gaussian
            update_ft_tag_Callback(hObject, eventdata, handles); % The scatter plot is updated        
            fit_zhist_tag_Callback(hObject, eventdata, handles); % This function performs the
                                                                 % 2-peaks
                                                                 % Gaussian
                                                                 % fit
            MyCoeffs_output = evalin('base','MyCoeffs_output');
            coef2(i,k) = MyCoeffs_output(2); % z1
            coef4(i,k) =MyCoeffs_output(5); % z2
        else % For circular fit
            update_ft_tag_Callback(hObject, eventdata, handles); % The scatter plot is updated
            axes(handles.axes5); 
            h = findobj(gcf, 'Type', 'scatter');
            r = h.XData'; % r data (read from plot)
            z = h.YData'; % z data (idem)           
            XY = [r z];
            circle_output = circle(XY); % circular fit
            rad = 2*circle_output(3); % diameter
            coef2(i,k) = rad;
            coef4(i,k) = circle_output(4);
            th = 0:pi/50:2*pi;
            xunit = 0.5*rad * cos(th)+circle_output(1);
            yunit = 0.5*rad * sin(th)+circle_output(2);
            hold on, plot(xunit,yunit,'LineWidth',1.5,'Color','g') % plot circle
            pause(0.25)
            clear r z
        end            
    end
end

figure(1)

% Plot (results)
for k = 1:length(alphas)
    if get(handles.ft_auto_method_tag,'Value') == 2 % 2-peaks Gaussian fit
        dif_z_k = abs(coef2(:,k)-coef4(:,k)); % dz = z1-z2
    else % circular fit
        dif_z_k = coef2(:,k);  % Diameter
    end
        hold on, plot(angles,dif_z_k);
        text(angles(end),dif_z_k(end),['Alpha = ',num2str(alphas(k))])
    end

xlimits = [angles(1)-2 angles(end)+2];
set(gca,'xlim',xlimits);
ylimits = get(gca,'ylim');
hold on, plot([angle angle],[ylimits(1) ylimits(2)]);
xlabel('Angle of incidence');
if get(handles.ft_auto_method_tag,'Value') == 2
    ylabel('\Deltaz / nm');
else
    ylabel('Diameter / nm');
end

    



% --- Executes on selection change in num_peaks_fitgs_tag.
function num_peaks_fitgs_tag_Callback(hObject, eventdata, handles)
% hObject    handle to num_peaks_fitgs_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)
% Hints: contents = cellstr(get(hObject,'String')) returns num_peaks_fitgs_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from num_peaks_fitgs_tag


% --- Executes during object creation, after setting all properties.
function num_peaks_fitgs_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to num_peaks_fitgs_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in fit_zhist_tag.
function MyCoeffs_output = fit_zhist_tag_Callback(hObject, eventdata, handles)
% This function performs the 1 or 2 peak(s) Gaussian fit to the relative-z
% histogram.
axes(handles.axes6)
   if exist('histogram') >0
    h = findobj(gcf, 'Type', 'histogram'); % Values and bin centres are read
                                           % from the histogram.
    b = h.Values;
    L = h.BinLimits;
    a = L(1):(L(2)-L(1))/(length(b)-1):L(2);

   else
    h = findobj(gcf, 'Type', 'patch');
    a = mean([h.XData(2,:);h.XData(3,:)]);
    b = mean([h.YData(2,:);h.YData(3,:)]);
    L = [a(1) a(end)];
   end
    xlimit =L;
    xl_fit = xlimit;
    N_pks = get(handles.num_peaks_fitgs_tag,'Value'); % Number of peaks (1 or 2)
    f = fit(a',b',['gauss',num2str(N_pks)],'Lower',[0 -inf 2 0 0 2],...
        'Upper',[inf 0 inf inf inf inf]); % Fit 
    xlim6 = xlim;
        hold on, plot(xlim6(1):0.5:xlim6(2),f(xlim6(1):0.5:xlim6(2)),...
            'LineWidth',1.5);
        pause(0.25)
    
    MyCoeffs_output = coeffvalues(f); % Coefficients obtained from fit
    assignin('base','MyCoeffs_output',MyCoeffs_output); % Fit results are
                                                    % exported to the workspace.
    set(handles.gauss_ft_z_out_tag,'string',['z1 = ',...
        num2str(0.1*round(10*MyCoeffs_output(2))),' nm; z2 = ',...
        num2str(0.1*round(10*MyCoeffs_output(5))),' nm']); % Fit results are shown
                                                           % in the GUI.
    



function z_bin_tag_Callback(hObject, eventdata, handles)
update_ft_tag_Callback(hObject, eventdata, handles)
% hObject    handle to z_bin_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of z_bin_tag as text
%        str2double(get(hObject,'String')) returns contents of z_bin_tag as a double


% --- Executes during object creation, after setting all properties.
function z_bin_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to z_bin_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in yn_circfit_tag.
function yn_circfit_tag_Callback(hObject, eventdata, handles)
% This function executes a circular fit of the 'Combined structures'
% scatter plot if the user checks the 'Fit circle' box. The result
% (diameter) is shown in the 'd/nm' box.
if get(handles.yn_circfit_tag,'Value') == 1
    axes(handles.axes5);
    h = findobj(gcf, 'Type', 'scatter');
    r = h.XData';
    z = h.YData';
    XY = [r z];
    circle_output = circle(XY);
    rad = circle_output(3);
    set(handles.circ_fit_tag,'string',num2str(2*0.1*round(10*rad)));
    th = 0:pi/50:2*pi;
    xunit = rad * cos(th)+circle_output(1);
    yunit = rad * sin(th)+circle_output(2);
    hold on, plot(xunit,yunit,'LineWidth',1.5,'Color','g');
else
update_ft_tag_Callback(hObject, eventdata, handles)
set(handles.circ_fit_tag,'string','-');
end

    
% hObject    handle to yn_circfit_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hint: get(hObject,'Value') returns toggle state of yn_circfit_tag



function circ_fit_tag_Callback(hObject, eventdata, handles)
% hObject    handle to circ_fit_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of circ_fit_tag as text
%        str2double(get(hObject,'String')) returns contents of circ_fit_tag as a double


% --- Executes during object creation, after setting all properties.
function circ_fit_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to circ_fit_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on selection change in ft_auto_method_tag.
function ft_auto_method_tag_Callback(hObject, eventdata, handles)
% hObject    handle to ft_auto_method_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: contents = cellstr(get(hObject,'String')) returns ft_auto_method_tag contents as cell array
%        contents{get(hObject,'Value')} returns selected item from ft_auto_method_tag


% --- Executes during object creation, after setting all properties.
function ft_auto_method_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to ft_auto_method_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: popupmenu controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end



function gauss_ft_z_out_tag_Callback(hObject, eventdata, handles)
% hObject    handle to gauss_ft_z_out_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'String') returns contents of gauss_ft_z_out_tag as text
%        str2double(get(hObject,'String')) returns contents of gauss_ft_z_out_tag as a double


% --- Executes during object creation, after setting all properties.
function gauss_ft_z_out_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to gauss_ft_z_out_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: edit controls usually have a white background on Windows.
%       See ISPC and COMPUTER.
if ispc && isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor','white');
end


% --- Executes on button press in calc_exc_prof_tag.
function calc_exc_prof_tag_Callback(hObject, eventdata, handles)
% This function is executed when the selected operation is 'Calculate
% excitation profile from background'. First, it verifies that the filename
% of the SMLM file exists, and then it runs 'obtain_bg.m' function. The
% resulting matrix (which contains the spatial information of the excitation
% profile pixelwise) is exported to a .csv file.
if exist(get(handles.filename_tag,'string'))
    obtain_bg(handles);
    ff = msgbox('Done! Excitation profile has been exported to .csv');
    set(handles.render_tag,'enable','off');
else
    w = msgbox(['Filename not found. Please check. If it is not ',...
        'located in the current directory, include the directory',...
        ' path in the filename: "C:\Dir\filename.hdf5"']);
end
    

% hObject    handle to calc_exc_prof_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)


% --- Executes on slider movement.
function slider_1_tag_Callback(hObject, eventdata, handles)
% This function controls the 'Marker size' slider.

slider_value = get(handles.slider_1_tag,'Value');

ax1_vis = get(handles.axes1,'visible');
ax2_vis = get(handles.axes7,'visible');

if sum(ax1_vis) == 221
    axes(handles.axes1);
    xlimits = xlim;
    h11 = findobj(gcf,'Type','scatter');
else
    axes(handles.axes7);
    xlimits = xlim;
    h11 = findobj(gcf,'Type','scatter');
end

    if sum(size(h11))>0
    h11(1).SizeData = (10+60*slider_value)*250/(xlimits(2)-xlimits(1));
    end
% hObject    handle to slider_1_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_1_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_1_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on slider movement.
function slider_3_tag_Callback(hObject, eventdata, handles)
% This function controls the 'Marker size' slider.
slider_value = get(handles.slider_3_tag,'Value');

axes(handles.axes3);
    xlimits = xlim;
    h33 = findobj(gcf,'Type','scatter');

    if sum(size(h33))>0
    h33(1).SizeData = (10+60*slider_value)*250/(xlimits(2)-xlimits(1));
    end
% hObject    handle to slider_3_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    structure with handles and user data (see GUIDATA)

% Hints: get(hObject,'Value') returns position of slider
%        get(hObject,'Min') and get(hObject,'Max') to determine range of slider


% --- Executes during object creation, after setting all properties.
function slider_3_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_3_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in set_small_roi_tag.
function set_small_roi_tag_Callback(hObject, eventdata, handles)
roi = drawrectangle('StripeColor','y');
assignin('base','roi',roi)
set(handles.angle_ROI_tag,'enable','on');
set(handles.angle_ROI_tag,'SliderStep',[0.001 0.1]);


% --- Executes on slider movement.
function angle_ROI_tag_Callback(hObject, eventdata, handles)
roi = evalin('base','roi');
slider_value = get(handles.angle_ROI_tag,'Value');

roi.RotationAngle = 360*slider_value;
 



% --- Executes during object creation, after setting all properties.
function angle_ROI_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to angle_ROI_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in plot_yz_tag.
function plot_yz_tag_Callback(hObject, eventdata, handles)
roi = evalin('base','roi');
simpler_output = evalin('base','simpler_output');
angle_roi = (270-(roi.RotationAngle))*pi/180;
positions_roi = roi.Position;
x1 = simpler_output(:,1);
y1 = simpler_output(:,2);
z1 = simpler_output(:,3);
[c,slope, origin] = find_points_ROI([x1 y1],...
    angle_roi, positions_roi,'y');

tita = atan(slope);

if (slope>0) && (slope ~= atan(pi/2))
    tita1= atan(abs((y1(c)-origin(2))./(x1(c)-origin(1))));
    r = ((x1(c)-origin(1)).^2+(y1(c)-origin(2)).^2).^(1/2);

    for i = 1:length(c)
        if (x1(c(i)) > origin(1)) && (y1(c(i)) > origin(2))
            tita2(i,1) = abs(tita-tita1(i));
        elseif (x1(c(i)) > origin(1)) && (y1(c(i)) < origin(2))
            tita2(i,1) = (tita+tita1(i));
        elseif (x1(c(i)) < origin(1)) && (y1(c(i)) > origin(2))
            tita2(i,1) = pi-(tita+tita1(i));
        end
    end
    proyec_r = cos(tita2).*r;
elseif (slope<0) && (slope ~= atan(-pi/2))
    tita = atan(-slope);
        tita1= atan(abs((y1(c)-origin(2))./(x1(c)-origin(1))));
    r = ((x1(c)-origin(1)).^2+(y1(c)-origin(2)).^2).^(1/2);

    for i = 1:length(c)
        if (x1(c(i)) > origin(1)) && (y1(c(i)) < origin(2))
            tita2(i,1) = abs(tita-tita1(i));
        elseif (x1(c(i)) > origin(1)) && (y1(c(i)) > origin(2))
            tita2(i,1) = (tita+tita1(i));
        elseif (x1(c(i)) < origin(1)) && (y1(c(i)) < origin(2))
            tita2(i,1) = pi-(tita+tita1(i));
        end
    end
    proyec_r = cos(tita2).*r;
end
   
if slope == 0
    proyec_r = x1(c)-min(x1(c));
end
if slope == atan(pi/2)
    proyec_r = y1(c)-min(y1(c));
end

assignin('base','roi_index_proyec_r', [c proyec_r]);


axes(handles.axes8);
set(handles.axes8,'visible','on');
scatter(simpler_output(c,1), simpler_output(c,2));
h1 = scatter(proyec_r,z1(c),50*ones(size(z1(c),1),1),z1(c),'o','filled');
xlabel('y (nm)');
ylabel('z (nm)');
daspect([1 1 1])
pbaspect([1 1 1])
xlimits_inicial = xlim;
h1.SizeData = 50*250/(xlimits_inicial(2)-xlimits_inicial(1));
daspect([1 1 1])
set(handles.marker_size_roi_text_tag,'visible', 'on'); 
set(handles.slider_markersize_roi_tag,'visible', 'on'); 
set(handles.export_roi_tag,'visible', 'on'); 
set(handles.render_tag,'enable','on');




% --- Executes on button press in plot_xz_tag.
function plot_xz_tag_Callback(hObject, eventdata, handles)

roi = evalin('base','roi');
simpler_output = evalin('base','simpler_output');
angle_roi = (270-(roi.RotationAngle))*pi/180;
positions_roi = roi.Position;
x1 = simpler_output(:,1);
y1 = simpler_output(:,2);
z1 = simpler_output(:,3);
[c,slope, origin] = find_points_ROI([x1 y1],...
    angle_roi, positions_roi,'x');

 

tita = atan(slope);

if (slope>0) && (slope ~= atan(pi/2))
    tita1= atan(abs((y1(c)-origin(2))./(x1(c)-origin(1))));
    r = ((x1(c)-origin(1)).^2+(y1(c)-origin(2)).^2).^(1/2);
    for i = 1:length(c)
        if (x1(c(i)) > origin(1)) && (y1(c(i)) > origin(2))
            tita2(i,1) = abs(tita-tita1(i));
        elseif (x1(c(i)) > origin(1)) && (y1(c(i)) < origin(2))
            tita2(i,1) = (tita+tita1(i));
        elseif (x1(c(i)) < origin(1)) && (y1(c(i)) > origin(2))
            tita2(i,1) = pi-(tita+tita1(i));
        end
    end
    proyec_r = cos(tita2).*r;
elseif (slope<0) && (slope ~= atan(-pi/2))
    tita = atan(-slope);  
    tita1= atan(abs((y1(c)-origin(2))./(x1(c)-origin(1))));
    r = ((x1(c)-origin(1)).^2+(y1(c)-origin(2)).^2).^(1/2);

    for i = 1:length(c)
        if (x1(c(i)) > origin(1)) && (y1(c(i)) < origin(2))
            tita2(i,1) = abs(tita-tita1(i));
        elseif (x1(c(i)) > origin(1)) && (y1(c(i)) > origin(2))
            tita2(i,1) = (tita+tita1(i));
        elseif (x1(c(i)) < origin(1)) && (y1(c(i)) < origin(2))
            tita2(i,1) = pi-(tita+tita1(i));
        end
    end
    max(tita1)
    max(tita2)
    proyec_r = cos(tita2).*r;
end
   
if slope == 0
    proyec_r = x1(c)-min(x1(c));
end
if slope == atan(pi/2)
    proyec_r = y1(c)-min(y1(c));
end

assignin('base','roi_index_proyec_r', [c proyec_r]);

axes(handles.axes8);
set(handles.axes8,'visible','on');
scatter(simpler_output(c,1), simpler_output(c,2));
h1 = scatter(proyec_r,z1(c),50*ones(size(z1(c),1),1),z1(c),'o','filled');
xlabel('x (nm)');
ylabel('z (nm)');
daspect([1 1 1])
pbaspect([1 1 1])
xlimits_inicial = xlim;
h1.SizeData = 50*250/(xlimits_inicial(2)-xlimits_inicial(1));
daspect([1 1 1])
set(handles.marker_size_roi_text_tag,'visible', 'on'); 
set(handles.slider_markersize_roi_tag,'visible', 'on'); 
set(handles.export_roi_tag,'visible', 'on'); 
set(handles.render_tag,'enable','on');



% --- Executes on slider movement.
function slider_markersize_roi_tag_Callback(hObject, eventdata, handles)
slider_value = get(handles.slider_markersize_roi_tag,'Value');

axes(handles.axes8);
xlimits = xlim;
h11 = findobj(gcf,'Type','scatter');
if sum(size(h11))>0
    h11(1).SizeData = (10+60*slider_value)*250/(xlimits(2)-xlimits(1));
end


% --- Executes during object creation, after setting all properties.
function slider_markersize_roi_tag_CreateFcn(hObject, eventdata, handles)
% hObject    handle to slider_markersize_roi_tag (see GCBO)
% eventdata  reserved - to be defined in a future version of MATLAB
% handles    empty - handles not created until after all CreateFcns called

% Hint: slider controls usually have a light gray background.
if isequal(get(hObject,'BackgroundColor'), get(0,'defaultUicontrolBackgroundColor'))
    set(hObject,'BackgroundColor',[.9 .9 .9]);
end


% --- Executes on button press in export_csv_roi_tag.
function export_csv_roi_tag_Callback(hObject, eventdata, handles)



% --- Executes on button press in export_roi_tag.
function export_roi_tag_Callback(hObject, eventdata, handles)
simpler_output_in = evalin('base','simpler_output');
roi_index_proyec_r = evalin('base','roi_index_proyec_r');
simpler_output = simpler_output_in(roi_index_proyec_r(:,1),:);
simpler_output = [simpler_output roi_index_proyec_r(:,2)];
export_function(4,simpler_output,handles);


% --- Executes on button press in browse_tag.
function browse_tag_Callback(hObject, eventdata, handles)
[file,path]= uigetfile({'*.hdf5;*.csv','.csv or .hdf5 files'},'Select a File');
set(handles.filename_tag,'string',fullfile(path,file))


% --- Executes on button press in browse_1_tag.
function browse_1_tag_Callback(hObject, eventdata, handles)
[file,path]= uigetfile('*.csv');
set(handles.excprof_filename_tag,'string',fullfile(path,file))
