function s = nansum(x, d)
x(isnan(x)) = 0; if nargin < 2, s = sum(x); else, s = sum(x, d); end
end
