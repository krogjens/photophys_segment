function [movies,optics,newLabelIm,sigProbs,fdrEst,frrEst] = trunc_seg(images,experiment,actions,threshCoeff)
    
    if nargin < 4
        threshCoeff = 1;
    end
% Filter with LoG 
     [logim,optics] = log_filter(images,experiment,actions);

% Binarize to get regions
     thedges = imbinarize(logim,0);
	 images.logim = thedges;
     [boundaries,labelIm] = bwboundaries(thedges,'holes'); 
     %We allow holes to be regions, since the middle of a bright region can have opposite curvature than the edge
    
% Obtain photophysics parameters (load or infer)
    % Moment calibration
    chipPars = get_chip_params(experiment.chipParsFile);
    
	% MLE for background strength and readout noise (uses logl)
    [lambda,chipPars,bgRat,noBgPx,threshInt,~,theta] = background_est(images,chipPars);
	
% Calculate p values for each region
    [newLabelIm,boundaries,sigProbs,fdrEst,frrEst] = filter_regions_bayes(images,labelIm,boundaries,actions,bgRat,threshInt,theta,threshCoeff);
   
% Store (and plot) each molecule
    movies = store_movies(images,newLabelIm,boundaries,sigProbs,experiment,actions);
    
% Paper plot with inset distribution
	figure(5)
	clf
	sortI = sort(images.imAverage(:));
	hbg = histogram(sortI(1:noBgPx));
	hbg.BinWidth = .999;
	hold on
	hol = histogram(sortI(noBgPx+1:end));
	hol.BinWidth = hbg.BinWidth;
	theta = [lambda 0 chipPars.gain/chipPars.adFactor chipPars.adFactor chipPars.roNoise chipPars.countOffset];
	maxInt = sortI(noBgPx);
	[cdfEnd,scale,cdf,offset] = SAD_cdf_inv(maxInt,1,theta);
	pdf = diff(cdf/cdfEnd)/scale;
	intVals = maxInt - offset(1):scale:maxInt+offset(2);
	pdf = pdf/sum(pdf(1:offset(1)));
	plot(intVals(2:end),hol.BinWidth*pdf*noBgPx,'--','Color','black','LineWidth',1)
	hold off
    xlim([min(images.imAverage(:))*0.8 1.7*median(images.imAverage(:))])
    xlabel('$n_{ic}$','Interpreter','latex','FontSize',12)
    ylabel('Counts','Interpreter','latex','FontSize',12)
    
    fig = gcf;
	fig.PaperUnits = 'inches';
	fig.PaperPosition = [0 0 4 3];
        
end
