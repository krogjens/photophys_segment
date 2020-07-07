function [NA,pixelSize,waveLength] = get_optic_params_direct(file)
fid = fopen(file);
data = textscan(fid,'%s%s%s%f%s%f%s%f');
fclose(fid);
foundWL = 0;
foundNA = 0; 
for j = 1:length(data{1})
    if strcmp(data{1}(j),'wavelength')
        waveLength = str2double(data{2}{j});
        foundWL = 1;        
    end
    if strcmp(data{1}(j),'pixelsize')
        pixelSize = str2double(data{2}{j});
        foundPS = 1;
    end
    if strcmp(data{1}(j),'NA')
        NA = str2double(data{2}{j});
        foundNA = 1;
    end
end
if ~foundWL
    fprintf('Did not find "wavelength" row in %s file \n',file);
end
if ~foundPS
    fprintf('Did not find "pixelsize" row in %s file \n',file);
end
if ~foundNA
    fprintf('Did not find "NA" row in %s file \n',file);
end
