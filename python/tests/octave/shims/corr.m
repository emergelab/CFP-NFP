function [r, p] = corr(X, Y, varargin)
% Shim for MATLAB Statistics Toolbox corr(X,Y,'rows','complete'[,'Type','Spearman']) -- Pearson/Spearman r and
% two-sided t-test p for every column of X against every column of Y. Used only to run the original CFP code in Octave.
if nargin < 2 || ischar(Y)
  % corr(X): plain Pearson matrix (what funaNRS uses) -- vectorised
  Xc = X - mean(X); s = sqrt(sum(Xc.^2)); r = (Xc' * Xc) ./ (s' * s); p = []; return
end
typ = 'Pearson';
for i = 1:2:numel(varargin), if strcmpi(varargin{i}, 'Type'), typ = varargin{i+1}; end, end
r = nan(size(X,2), size(Y,2)); p = r;
for a = 1:size(X,2)
  for b = 1:size(Y,2)
    ok = ~isnan(X(:,a)) & ~isnan(Y(:,b)); x = X(ok,a); y = Y(ok,b);
    if strcmpi(typ, 'Spearman'), x = tiedrank(x); y = tiedrank(y); end
    n = numel(x); xc = x - mean(x); yc = y - mean(y);
    rr = sum(xc.*yc) / sqrt(sum(xc.^2) * sum(yc.^2));
    r(a,b) = rr;
    t = rr * sqrt((n-2) / (1 - rr^2));
    p(a,b) = 2 * tcdf(-abs(t), n-2);
    if abs(rr) == 1, p(a,b) = 0; end
  end
end
end
