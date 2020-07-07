function [cdf,scale,cdf_I,offset] = SAD_cdf_inv(I,nPx,theta)

lambda = theta(1)*nPx;
r = theta(3);
f = theta(4);
if nPx > 1
	sig = theta(5)*sqrt(nPx);
else
	sig = theta(5);
end
delta = theta(6)*nPx;
offset = zeros(1,2);

meanI = lambda*r+delta;
varI = (sig/f)^2 + 2*lambda*r^2;

ymin = max(0,min(I,meanI-6*sqrt(varI)));
ymin = round(ymin);
offset(1) = round(I - ymin);
ymax = max(I,meanI+6*sqrt(varI));
ymax = round(ymax);
offset(2) = ymax - I;
range = ymax - ymin;
rmax = 5000;
if range > rmax
    scale = range/rmax;
    range = rmax;
    offset(1) = round(offset(1)/scale);
%    fprintf('CDF has been scaled to reduce computation.\n')
else
	scale = 1;
end
y = round(linspace(ymin,ymax,range+1));
cfY_I = @(t) exp(-t.^2*(sig/f)^2/2 + lambda./(1-1i*r*t) - lambda + 1i*t*delta);
N = 2^6;
dt = 2*pi/(ymax-ymin);%2*pi/range; 
t = (1:N)' * dt;
cft = cfY_I(t);
cft(N) = cft(N)/2;
[n,m] = size(y);
y = y(:);
E = exp(-1i*y*t');

cdf_I = (meanI - y)/2 + imag(E * (cft ./ t));
cdf_I = 0.5 - (cdf_I * dt)/pi;
cdf_I = reshape(max(0,cdf_I),n,m);
cdf_I = cdf_I/cdf_I(end);
cdf = cdf_I(1+offset(1)); 
end
