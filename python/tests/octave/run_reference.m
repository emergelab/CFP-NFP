% Runs the repository's own MATLAB functions (getPMs/) under Octave on Example/TS.mat with recorded folds and
% permutations (shims/ replace the Statistics/Bioinformatics toolbox calls), writing reference.mat for
% test_parity.py. Called as: octave --eval "repo='...'; work='...'; here='<this folder>'; run(...)"
pkg load statistics
addpath(fullfile(here, 'shims'));
addpath(fullfile(repo, 'getPMs'));
global CFP_FOLDS CFP_FI CFP_PERMS CFP_PI
load(fullfile(repo, 'getPMs', 'atlasindex.mat')); C = atlasindex.A424.AAC;
load(fullfile(work, 'inputs.mat'));
nc = numel(ts); levels = [5 20 43 137];
ref = struct();
for z = 0:1
  ni = zeros(424, numel(levels), nc); ne = ni; tn = cell(1, nc);
  for c = 1:nc
    [ni(:,:,c), ne(:,:,c), tn{c}] = funaNRS(ts{c}, C, levels, 0, z);
  end
  tt = zeros(size(tn{1},1), numel(levels), nc); for c = 1:nc, tt(:,:,c) = tn{c}; end
  ref.(sprintf('ni_z%d', z)) = ni; ref.(sprintf('ne_z%d', z)) = ne;
  ref.(sprintf('tnrs_z%d', z)) = tt(1:1275, :, :);
end
call = 0;
for li = 1:3
  m = max(C(:, levels(li))); midx = triu(ones(m)); midx = logical(midx(:)); F = zeros(m, m, nc);
  for c = 1:nc, A = zeros(m); A(midx) = ref.tnrs_z0(1:sum(midx), li, c); F(:,:,c) = A; end
  for model = {'wcpm', 'cpm'}
    call = call + 1;
    CFP_FOLDS = {folds_single(call,:)}; CFP_FI = 0;
    [r, p, pos, neg] = funPM(F, behav(:), model{1}, 10, 0.05, 0, 1, 1e-6, 1);
    ref.(sprintf('pm_%s_L%d', model{1}, levels(li))) = struct('r', r, 'pos', pos, 'neg', neg);
    CFP_FOLDS = {folds_single(call,:)}; CFP_FI = 0;
    [r, p, pos, neg] = funPM2(F, behav(:), model{1}, 10, 0.05, 0, 1, 1e-6, 1);
    ref.(sprintf('pm2_%s_L%d', model{1}, levels(li))) = struct('r', r, 'pos', pos, 'neg', neg);
  end
end
m = max(C(:, 20)); midx = triu(ones(m)); midx = logical(midx(:)); F = zeros(m, m, nc);
for c = 1:nc, A = zeros(m); A(midx) = ref.tnrs_z0(1:sum(midx), 2, c); F(:,:,c) = A; end
CFP_FOLDS = num2cell(folds_iter, 2); CFP_FI = 0; CFP_PERMS = num2cell(perms, 2); CFP_PI = 0;
[r, p, pos, neg] = funpPM(F, behav(:), 5, 'wcpm', 10, 0.05, 0, 1, 0, 1e-6, 1);
ref.ppm_wcpm_L20 = struct('r', r, 'p', p, 'pos', pos, 'neg', neg);
save('-v7', fullfile(work, 'reference.mat'), 'ref');
