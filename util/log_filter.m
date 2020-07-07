function [logim,optics] = log_filter(images,experiment,actions)

imAverage = images.imAverage;
if isfield(experiment,'opticsFile')
	[optics.NA,optics.pixelSize,optics.waveLength] = get_optic_params_direct(experiment.opticsFile);
else
	folder = experiment.targetFolder;
	[optics.NA,optics.pixelSize,optics.waveLength] = get_optic_params(folder);
end
optics.sigma = 1.22*optics.waveLength/(2*optics.NA) / optics.pixelSize;
if images.imageNumber == 1
	fprintf('Width of point spread function estimated to be %.2f pixels.\n',optics.sigma);
end

n = ceil(6*optics.sigma);
n = n + 1 - mod(n,2);
filt = fspecial('log',n,optics.sigma);
if ~isa(imAverage,'double')
	imAverage = double(imAverage);
	fprintf('Averaged image has been converted to type "double".\n');
end
logim = imfilter(imAverage,filt,'replicate');

if actions.showMolecules
	figure(1+(images.imageNumber-1)*5)
	imshow(imbinarize(logim,0),'InitialMagnification','fit')
end
