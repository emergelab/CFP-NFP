function m = nanmean(x, d)
if nargin < 2, d = find(size(x) ~= 1, 1); if isempty(d), d = 1; end, end
ok = ~isnan(x); x(~ok) = 0; m = sum(x, d) ./ sum(ok, d);
end
