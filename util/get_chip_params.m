function [chipPars] = get_chip_params(chipParsFile)

fid = fopen(chipParsFile);
data = textscan(fid,'%s%s%s%f%s%f%s%f');
fclose(fid);
foundWL = 0;
foundNA = 0;
foundCO = 0; 
for j = 1:length(data{1})
    if strcmp(data{1}(j),'gain')
        chipPars.gain = str2double(data{2}{j});
        foundGa = 1;        
    end
    if strcmp(data{1}(j),'adfactor')
        chipPars.adFactor = str2double(data{2}{j});
        foundAF = 1;
    end
    if strcmp(data{1}(j),'countoffset')
        chipPars.countOffset = str2double(data{2}{j});
        foundCO = 1;
    end
	if strcmp(data{1}{j},'ronoise')
		chipPars.roNoise = str2double(data{2}{j});
	end
end
if ~foundGa
    fprintf('Did not find "gain" row in %s file \n',chipParsFile);
end
if ~foundAF
    fprintf('Did not find "adfactor" row in %s file \n',chipParsFile);
end
if ~foundCO
    fprintf('Did not find "countoffset" row in %s file \n',chipParsFile);
end
