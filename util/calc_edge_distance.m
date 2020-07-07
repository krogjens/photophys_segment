function edgeDist = calc_edge_distance(images,labelIm,boundaries)

% Calculate regional image counts and store region sizes
N = length(boundaries);
edgeDist = zeros(1,N);
dims = size(labelIm);

for k = 1:N
	yDist = min(dims(1) - max(boundaries{k}(:,1)),min(boundaries{k}(:,1)));
	xDist = min(dims(2) - max(boundaries{k}(:,2)),min(boundaries{k}(:,2)));
	edgeDist(k) = min(xDist,yDist);
end
end
