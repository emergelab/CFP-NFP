function [score, r, pval, pred] = applyFP(features, posW, negW, behav, k)
% Apply a published fingerprint (e.g. FPs/AD_CFP.mat, FPs/AD_NFP.mat) to new subjects.
%
% features : new subjects' features with subjects in the LAST dimension, laid out like the
%            fingerprint (e.g. m x m x n NRS matrices for a CFP level, or 424 x n for an NFP column)
% posW     : positive fingerprint weights for that level/measure (same shape as one subject)
% negW     : negative fingerprint weights (same shape)
% behav    : (optional) n x 1 outcome for the new subjects
% k        : (optional) folds for calibrated predictions in the new data (default 10)
%
% score    : n x 1 fingerprint score = sum(features .* posW) - sum(features .* negW)
% r, pval  : Pearson correlation of score with behav (scale-free; no fitting needed)
% pred     : n x 1 cross-validated predictions on the scale of behav: the score -> behav line is
%            fitted on the training folds of the NEW data and applied to the held-out fold
%
% The score itself has arbitrary units; use r/pval, or pred, rather than reading the score as an outcome.
%
% Inputs must match how the fingerprint was built. For the published FPs/AD_CFP.mat and FPs/AD_NFP.mat:
%   - A424 parcellation (Glasser cortex + subcortex + cerebellum), HCP-style preprocessing incl. ICA-FIX and
%     global (mean grayordinate) signal regression;
%   - the same feature definition and AA level as the fingerprint (e.g. AA-24 NRS from funaNRS for AA24_NRSv);
%   - CHANGE scores: week-1 minus baseline (pre-treatment), not single-session connectivity;
%   - the outcome as in the original analysis (percent improvement; treatment contrast when comparing arms).
% Applied to differently defined features, the score is still computed but is not meaningful.

if nargin < 5 || isempty(k), k = 10; end
dim = size(features);
n = dim(end);
X = reshape(features, [], n);
P = posW(:); N = negW(:);
P(isnan(P)) = 0; N(isnan(N)) = 0;
X(isnan(X)) = 0;
score = (P' * X - N' * X)';

r = NaN; pval = NaN; pred = nan(n, 1);
if nargin < 4 || isempty(behav)
    return
end
behav = behav(:);
ok = ~isnan(behav);
[r, pval] = corr(score(ok), behav(ok));

if nargout > 3
    idx = find(ok);
    k = min(k, numel(idx));
    folds = mod(randperm(numel(idx)), k) + 1;
    for f = 1:k
        te = idx(folds == f); tr = idx(folds ~= f);
        b = polyfit(score(tr), behav(tr), 1);
        pred(te) = b(1) * score(te) + b(2);
    end
end
end
