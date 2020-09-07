function [c_frame, c_x, c_y, c_ph, c_off] = get_column_idx_ThunderStorm(filename_wformat)

data = readtable(filename_wformat);
headers_thunder_storm = data.Properties.VariableNames;

c_head_x = zeros(size(headers_thunder_storm,2),1);
c_head_y = c_head_x;
c_head_ph = c_head_x;
c_frame = c_head_x;
c_offset = c_head_x;

for i = 1:size(headers_thunder_storm,2)
    length_x_str = length('x_nm_');
    length_str_i = length(num2str(headers_thunder_storm{1,i}));    
    if length_x_str == length_str_i
        c_head_x(i,1) = sum(num2str(headers_thunder_storm{1,i}) == 'x_nm_') == length_x_str;
    end
end
for i = 1:size(headers_thunder_storm,2)
    length_y_str = length('y_nm_');
    length_str_i = length(num2str(headers_thunder_storm{1,i}));    
    if length_y_str == length_str_i
        c_head_y(i,1) = sum(num2str(headers_thunder_storm{1,i}) == 'y_nm_') == length_y_str;
    end
end
for i = 1:size(headers_thunder_storm,2)
    length_ph_str = length('intensity_photon_');
    length_str_i = length(num2str(headers_thunder_storm{1,i}));    
    if length_ph_str == length_str_i
        c_head_ph(i,1) = sum(num2str(headers_thunder_storm{1,i}) == 'intensity_photon_') == length_ph_str;
    end
end
for i = 1:size(headers_thunder_storm,2)
    length_frame_str = length('frame');
    length_str_i = length(num2str(headers_thunder_storm{1,i}));    
    if length_frame_str == length_str_i
        c_head_frame(i,1) = sum(num2str(headers_thunder_storm{1,i}) == 'frame') == length_frame_str;
    end
end
for i = 1:size(headers_thunder_storm,2)
    length_offset_str = length('offset_photon_');
    length_str_i = length(num2str(headers_thunder_storm{1,i}));    
    if length_offset_str == length_str_i
        c_head_offset(i,1) = sum(num2str(headers_thunder_storm{1,i}) == 'offset_photon_') == length_offset_str;
    end
end

c_x = find(c_head_x);
c_y = find(c_head_y);
c_ph = find(c_head_ph);
c_frame = find(c_head_frame);
if exist('c_head_offset')
c_off = find(c_head_offset);
else
    c_off = ones(length(c_x),1);
end

