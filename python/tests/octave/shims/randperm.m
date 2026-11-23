function v = randperm(n, varargin)
% Shim: pre-recorded permutations (global CFP_PERMS) for funpPM's null.
global CFP_PERMS CFP_PI
CFP_PI = CFP_PI + 1; v = CFP_PERMS{CFP_PI}(:)'; assert(numel(v) == n);
end
