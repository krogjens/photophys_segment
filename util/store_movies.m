function movies = store_movies(images,newLabelIm,boundaries,sigProbs,experiment,actions)

hasDots = isfield(images,'dotIm');
if hasDots
	dotIm = images.dotIm;
end
edgePx = 1;
% Generate molecule images for barcode extraction
if hasDots
	[molM,bwM,dotM,pos] = generate_molecule_images_fast(boundaries,newLabelIm,images.registeredIm,dotIm,experiment.targetFolder,images.runNo,images.imageName,edgePx,actions);
else
	[molM,bwM,dotM,pos] = generate_molecule_images_fast(boundaries,newLabelIm,images.registeredIm,[],experiment.targetFolder,images.runNo,images.imageName,edgePx,actions);
end
movies.imageName = images.imageName;
movies.molM = molM;
movies.bwM = bwM;
movies.pos = pos;

% Set up figures and plots of molecules and pvalues (Needed for
% compatibility with mark_bars etc)
if actions.showMolecules
    molFigNum = 30+(images.imageNumber-1)*5;
    movies.molFigNum = molFigNum;
    figure(molFigNum)
    clf
    ax1 = axes;
    imshow(mat2gray(images.imAverage),'InitialMagnification','fit');
    parMap = parula;
    ax1.Visible = 'off';
    ax1.XTick = [];
    ax1.YTick = [];
	hold on
    
    for k = 1:length(boundaries)
        lineColor = parMap(ceil((sigProbs(k)-0.5)*2*256),:);
		plot(boundaries{k}(:,2),boundaries{k}(:,1),'LineWidth',1,'Color',lineColor);
    end
    hold off
    ax2 = axes;
    linkaxes([ax1,ax2])
    ax2.Visible = 'off';
    ax2.XTick = [];
    ax2.YTick = [];
   
    cb = colorbar(ax2,'Ticks',[0,.2,0.4,0.6,0.8,1],...
         'TickLabels',{0.5,0.6,0.7,0.8,0.9,1},'Position',[.85 .11 .055 .815]); %.815
    
    clabeltext = text(575,320,'$p(\mathbf{s}|n_d,n_b)$','Color','black','Interpreter','latex');
    set(clabeltext,'Rotation',-90)
end

if hasDots 
	movies.dotM = dotM;
    if actions.showMolecules
        dotFigNum = 4+(images.imageNumber-1)*5;
        movies.dotFigNum = dotFigNum;
    end
end
