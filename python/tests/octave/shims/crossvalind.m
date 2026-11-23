function idx = crossvalind(~, n, k)
% Shim: hands out pre-recorded fold labels (global CFP_FOLDS, a cell queue) so Octave and Python use identical folds.
global CFP_FOLDS CFP_FI
CFP_FI = CFP_FI + 1; idx = CFP_FOLDS{CFP_FI}(:); assert(numel(idx) == n);
end
