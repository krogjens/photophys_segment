function segmentation_skeleton()

im = imread('lung_sample.tif'); % Specify target image
experiment.targetFolder = 'testfolder/'; % Specify folder to store molecule images
experiment.opticsFile = 'lund_20x_optics.txt'; % Optionally specified optics file location (full path needed). If this is specified, the software assumes the optical setup to be identical for all images
experiment.chipParsFile = 'lund_gain300_chippars.txt'; % Specify file with 

actions.showScores = 1; 			% Show a histogram of the scores for the regions as well as the dots, which helps setting tweak parameters
actions.showMolecules = 1;			% Show plots of detected molecules and images
actions.saveMolecules = 0;			% Save individual molecules

images.imAverage = double(im);
images.registeredIm{1} = double(im);
images.imageName = 'LungSample';
images.imageNumber = 1;
images.runNo = 1;

addpath('util')
t = cputime;
[segmentImgs,~,newLabelIm,sigProbs,fdrEst,~] = trunc_seg(images,experiment,actions);
e = cputime - t;
output.segments = segmentImgs;
output.sigProbs = sigProbs;
output.fdrEst = fdrEst;
save('segmentationoutput','output')

fprintf('Found %i regions .\n',max(newLabelIm(:)));
fprintf('Total segmentation time: %.2f.s\n', e);
end
