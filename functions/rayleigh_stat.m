function [r, pval] = rayleigh_stat(phi)
%calculate rayleigh statistic for circular distribution of phases

% phi should be a vector
n = numel(phi);
r = (sqrt(sum(cos(phi)).^2) + sum(sin(phi)).^2)/n;

% compute p value using approxation in Zar, p. 617
pval = exp(sqrt(1+4*n+4*(n^2-r^2))-(1+2*n));

end

