function [NA,pixelSize,waveLength] = get_optic_params(folder)

content = dir(folder);
if length(content) < 1
	fprintf('No folder found at location "%s"',folder);
end
i = 1;
while i < length(content) + 1
    isoptics = strcmp(content(i).name,'optics.txt');
    if isoptics
        break;
    end
    i = i+1;
end
if ~isoptics
    fprintf('Error: optic parameter file not found in %s folder',folder);
end
fid = fopen([folder,content(i).name]);
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
    fprintf('Did not find "wavelength" row in %s file \n',content(i).name);
end
if ~foundPS
    fprintf('Did not find "pixelsize" row in %s file \n',content(i).name);
end
if ~foundNA
    fprintf('Did not find "NA" row in %s file \n',content(i).name);
end
