function [lambda,chipPars,bgRat,m,threshInt,nOut,theta] = background_est(images,chipPars)

gain = chipPars.gain;
adFactor = chipPars.adFactor;
countOffset = chipPars.countOffset;
roNoise = chipPars.roNoise;
r = gain/adFactor;

% Sort data
sortI = sort(images.imAverage(:));
m = numel(sortI);

% Specify FDR
alpha = .5;

% Specify likelihood function
invpdf = @(data,lambda) SAD_logl_inv(data,[lambda 0 r adFactor roNoise countOffset]); % Specify likelihood function

% Intialize
hasOutliers = 1;
nOutliers = 0;
runs = 0;
while hasOutliers
    runs = runs + 1;
	% Remove outliers
	m = m - nOutliers;
	truncI = sortI(1:m);
	% Fit lambda
	lamGuess = abs((truncI(round(m/2)) - countOffset)/r);
	lambda = mle(truncI,'logpdf',invpdf,'start',lamGuess,'lowerBound',0);
	% Calculate p-values
	theta = [lambda 0 r adFactor roNoise countOffset];
	pVals = SAD_pval_inv(truncI,theta); 
	pVals = fliplr(pVals);
	% Find outliers
	threshold = (1:m)/m*alpha;
	outliers = find(pVals < threshold);
	if ~isempty(outliers)
		nOutliers = outliers(end);
		hasOutliers = 1;
	else
		hasOutliers = 0;
	end
end
bgRat = 1-pVals(1);
threshInt = sortI(m);
nOut = numel(sortI) - m;
fprintf('Estimated (certain) background (n_ic <= %i) area fraction: %.4f.\n',threshInt,1-nOut/numel(sortI))
fprintf('Estimated value of lambda = %.2f\n',lambda)

end
