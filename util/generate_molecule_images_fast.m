function [molM,bwM,dotM,pos] = generate_molecule_images_fast(D,L,registeredIm,dotIm,folder,runNo,imageName,edgePx,actions)
cellLength = length(D);
pos = cell(1,cellLength);
molM = cell(1,cellLength);
bwM = cell(1,cellLength);
dotM = cell(1,cellLength);

for k = 1:cellLength
    yMin = max(1,min(D{k}(:,1))-edgePx);
    yMax = min(edgePx+max(D{k}(:,1)),size(registeredIm{1},1));
    xMin = max(1,min(D{k}(:,2))-edgePx);
    xMax = min(edgePx+max(D{k}(:,2)),size(registeredIm{1},2));
	dY = yMax-yMin+1;
	dX = xMax-xMin+1;
    molMov = nan(dY,dX,length(registeredIm));
    bwMov = nan(dY,dX);

	 labelCut = L(yMin:yMax,xMin:xMax);
     mask = double(labelCut == k);
	 for l = 1:length(registeredIm)
		 regCut = registeredIm{l}(yMin:yMax,xMin:xMax);
	 	 molMov(:,:,l) = double(regCut).*mask;
		 molMov(molMov == 0) = NaN;
	 end
	 if isempty(dotIm)
		 dotMov = [];
	 else
		 dotCut = dotIm(yMin:yMax,xMin:xMax);
		 dotMov = double(dotCut).*mask;
         dotMov(dotMov == 0) = NaN;
	 end
	 bwMov(mask == 1) = 1;
 
	 if actions.saveMolecules
         outputFolderName = [folder,'/molecules_run',num2str(runNo),'/'];
         if ~exist(outputFolderName,'dir')
             mkdir(outputFolderName);
         end
		 molSaveName = [folder,'/molecules_run',num2str(runNo),'/',imageName,'_mol',num2str(k),'mov.tif'];
		 imwrite(uint16(molMov(:,:,1)),molSaveName);
		 for i = 2:length(registeredIm)
   	     	imwrite(uint16(molMov(:,:,i)),molSaveName,'writemode','append');
		 end
		 if ~isempty(dotMov)
			 dotSaveName = [folder,'molecules_run',num2str(runNo),'/',imageName,'_dot',num2str(k),'mov.tif'];
			 imwrite(uint16(dotMov),dotSaveName);
		 end
		 bwSaveName = [folder,'molecules_run',num2str(runNo),'/',imageName,'_bwmol',num2str(k),'mov.tif'];	 
		 imwrite(uint16(bwMov(:,:,1)),bwSaveName);
	 end
%	 pos{k} = D{k}(1,:);
	 pos{k} = [yMin xMin];
  	 molM{k} = molMov;
	 bwM{k} = bwMov;
	 dotM{k} = dotMov;
end

molM = molM(~cellfun('isempty',molM));
bwM = bwM(~cellfun('isempty',bwM));
if isempty(dotIm)
	dotM = [];
else
	dotM = dotM(~cellfun('isempty',dotM));
end
