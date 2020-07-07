function logl = SAD_logl_inv(truncI,theta)

nI = numel(truncI); 
[pdf,cdf,offset] = inv_cf(truncI,theta);
oldI = truncI(1);
pdfIdx = 1 + offset(1);
loglcont = log(pdf(pdfIdx));
logl = loglcont;
for idx = 2:nI
	newI = truncI(idx);
	if newI > oldI
		pdfIdx = pdfIdx + newI - oldI;
		loglcont = log(pdf(pdfIdx));
		oldI = newI;
	end
	logl = logl + loglcont;
end
logl = logl - nI*log(cdf(end-offset(2)));

end

function [pdf,cdf_I,offset] = inv_cf(truncI,theta)

lambda = theta(1);
r = theta(3);
f = theta(4);
sig = theta(5);
delta = theta(6);
offset = zeros(1,2);

meanI = lambda*r+delta;
varI = (sig/f)^2 + 2*lambda*r^2;
ymin = max(0,min(truncI(1),meanI-6*sqrt(varI)));
ymin = round(ymin);
offset(1) = truncI(1) - ymin;
ymax = max(truncI(end),meanI+6*sqrt(varI));
ymax = round(ymax);
offset(2) = ymax - truncI(end);
range = ymax - ymin;
y = round(linspace(ymin,ymax,range+1));
cfY_I = @(t) exp(-t.^2*(sig/f)^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*delta);
N = 2^8;
dt = 2*pi/range;
t = (1:N)' * dt;
cft = cfY_I(t);
cft(N) = cft(N)/2;
[n,m] = size(y);
y = y(:);
E = exp(-1i*y*t');
pdf = 0.5 + real(E * cft);
pdf = (pdf*dt)/pi;
%pdf = reshape(max(0,pdf),n,m);
pdf = reshape(max(1e-14,pdf),n,m);

cdf_I = (meanI - y)/2 + imag(E * (cft ./ t));
cdf_I = 0.5 - (cdf_I * dt)/pi;
cdf_I = reshape(max(0,cdf_I),n,m);
if cdf_I(end-offset(2)) < 1e-10
%	fprintf('Warning: Found cdf = %e at truncation point\n',cdf_I(end-offset(2)));
%	cdf_I(end-offset(2)) = 1; % This is wrong but helps punish bad parameter values
end
end
