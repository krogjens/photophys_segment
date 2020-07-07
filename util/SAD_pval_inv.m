function [pvals] = SAD_pval_inv(sortI,theta)

m = numel(sortI);
pvals = zeros(1,m);
[cdf,offset] = inv_cf(sortI,theta);
oldI = sortI(1);
cdfIdx = 1 + offset(1);
pvals(1) = 1 - cdf(cdfIdx);
for idx = 2:m
	newI = sortI(idx);
	diff = newI - oldI;
	if diff > 0
		cdfIdx = cdfIdx + diff;
		oldI = newI;
	end
	pvals(idx) = 1 - cdf(cdfIdx-1);
end

end

function [cdf_I,offset] = inv_cf(sortI,theta)

lambda = theta(1);
r = theta(3);
f = theta(4);
sig = theta(5);
delta = theta(6);
m = numel(sortI);

meanI = lambda*r+delta;
varI = (sig/f)^2 + 2*lambda*r^2;

ymin = max(0,min(sortI(1),meanI-6*sqrt(varI)));
ymin = round(ymin);
offset(1) = round(sortI(1) - ymin);
ymax = max(sortI(m),meanI+6*sqrt(varI));
ymax = round(ymax);
offset(2) = ymax - sortI(m);
range = ymax - ymin;
y = round(linspace(ymin,ymax,range+1));
cfY_I = @(t) exp(-t.^2*(sig/f)^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*delta);
N = 2^12;
dt = 2*pi/range; 
t = (1:N)' * dt;
cft = cfY_I(t);
cft(N) = cft(N)/2;
[n,m] = size(y);
y = y(:);
E = exp(-1i*y*t');

cdf_I = (meanI - y)/2 + imag(E * (cft ./ t));
cdf_I = 0.5 - (cdf_I * dt)/pi;
cdf_I = reshape(max(0,cdf_I),n,m);

%cdf_I = cdf_I/cdf_I(end);
%cdf = cdf_I(1+offset(1)); 
end
