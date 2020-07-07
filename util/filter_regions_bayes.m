function [accLabelIm,accBoundaries,sigProbs,fdrEstimate,frrEstimate] = filter_regions_bayes(images,labelIm,boundaries,actions,bgRat,threshInt,theta,threshCoeff)
   
    nReg = length(boundaries);
    % The naive Bayesian classifier approach 

    iMax = max(images.imAverage(:));
    tt = graythresh(images.imAverage/iMax);
    segThresh = threshCoeff*tt*iMax;
    
    nPx =  numel(images.imAverage);
    pMin = 1/ (nPx + 2);
    segCDF = 1-SAD_pval_inv(segThresh+1,theta);

    p_bg_black =  segCDF;
    if p_bg_black > 1
        fprintf('p_bg_black found > 1: %f\n',p_bg_black);
        p_bg_black = 1-1/(nPx + 2);
    elseif p_bg_black < pMin
        fprintf('p_bg_black found < pMin: %f\n',p_bg_black)
        p_bg_black = 1/(nPx + 2);
    end
    p_bg_white = 1-p_bg_black;
    
    p_s_black = 0.5;
    p_s_white = 1-p_s_black;
    
    [significant,sigProbsTemp] = classify_regions(images,labelIm,nReg,segThresh,p_bg_black,p_bg_white,p_s_black,p_s_white);
   
    edgeDist = calc_edge_distance(images,labelIm,boundaries);
	edgeTresh = 1; % Minimum distance to edge for regions
    
	nonEdge = edgeDist > edgeTresh;
	accepted = significant & nonEdge;
    frrEstimate = mean(sigProbsTemp(~accepted));
    
    newLabelIm = zeros(size(labelIm));
    regIdx = 0;
    for k = 1:nReg
        if accepted(k)
            regIdx = regIdx + 1;
            newLabelIm(labelIm == k) = regIdx;
        end
    end
    
    acceptIm = imbinarize(newLabelIm,0.5);
    [accBoundaries,accLabelIm] = bwboundaries(acceptIm,'noholes'); % Since regions have been considered, we dont allow holes to be new regions.
    nRegReduced = size(accBoundaries,1);
    
    [~,sigProbs] = classify_regions(images,accLabelIm,nRegReduced,segThresh,p_bg_black,p_bg_white,p_s_black,p_s_white);
    fdrEstimate = mean(1-sigProbs);
   
    
	if actions.showScores
		%figure(2+(images.imageNumber-1)*5)
		%h1 = histogram(pValAll,0:0.1:1.1);
		%title([images.imageName,' regional p values'])
		%xlabel('P')
%	    line([pTresh pTresh],[0 max(h1.Values)],'LineStyle','--','Color','black','LineWidth',2)
%	    text(1.1*pTresh,2/3*max(h1.Values),'User set threshold','FontSize',14)
% 		figure(99)
% 		clf
% 		histogram(pDist(abs(pDist) < pDistThresh))
% 		vv = var(pDist(abs(pDist) < pDistThresh))
% 		xlim([-15 40])
% 		xlabel('$d$','Interpreter','latex')
% 		ylabel('Count','Interpreter','latex')
% 		axes('Position',[0.6 0.6 0.25 0.25])
% 		box on
%		ylim([0 5])
%		figure(100)
		%hold on
		%histogram(pDist(pDist>pDistThresh))
		%hold off
		%xlabel('d','Interpreter','latex')
%		xlim([p200 2000])
		%fig = gcf;
		%fig.PaperUnits = 'inches';
		%fig.PaperPosition = [0 0 4 3];
%		print('~/Dropbox/dnabarcoding/ifm/draft/gain100beads_regdists','-depsc')
	end
end

function [significant,probs] = classify_regions(images,labelIm,nReg,segThresh,p_bg_black,p_bg_white,p_s_black,p_s_white)
    significant = zeros(1,nReg);
    probs = zeros(1,nReg);
    for i = 1:nReg % for each region
        regIdx = labelIm == i;
        counts = images.imAverage(regIdx);
        regSize = numel(counts);
        n_black = sum(counts <= segThresh); %number of pixels less than I_th
        n_white = regSize - n_black;
        log_p_bg = n_black * log(p_bg_black)  + n_white * log(p_bg_white); % Up to a constant that is equal for log_p_s
        % and VV for p_s;
        log_p_s = n_white * log(p_s_white) + n_black * log(p_s_black);
        probSig = exp(log_p_s - logsumexp2(log_p_s,log_p_bg));
        probs(i) = probSig; 
        if log_p_s > log_p_bg %, set signal
            significant(i) = 1;
        else
            significant(i) = 0;
        end
    end
end
