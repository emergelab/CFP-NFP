function [R_both, pval_both, pos_features, neg_features, coef_features] =...
    funpPM(features, behav, iterations, model, k, thresh,...
    spearman, s, d, v_alpha, lambda)

% features = input in 2D, 3D, 4D with subjects in last dimension 
% set defaults
    if ~exist('iterations', 'var')
        iterations = 200;
    end
        
    if ~exist('model', 'var')
        model = 'cpm'; % wcpm (weighted) or rcpm (ridge)
    end
        
    if ~exist('k', 'var')
        k = 10; % number of cross validations (CV)
    end
    
    if ~exist('thresh', 'var')
        thresh = 0.05; % p value to select features.
    end
    
    if ~exist('spearman', 'var')
        spearman = 0; % 0 (pearson) or 1 (spearman)
    end

    if ~exist('s', 'var')
        s = 1; % sparsity 0 to 1; 1 -> features retained in all CV
    end
        
    if ~exist('d', 'var')
        d = 0; % density 0 to 1; 0 -> features retained in any iteration
    end
    
    if ~exist('v_alpha', 'var') || length(v_alpha) ~= 1 
        v_alpha = 1e-6;
    end
    
    if ~exist('lambda', 'var') || length(lambda) ~= 1 
        lambda = '';
    end
    
 % run PMs
 if (iterations == 0)
    [R_both, pval_both, pos_features, neg_features, coef_features] =...
        funPM(features, behav, model, k, thresh, spearman, s,...
        v_alpha, lambda);
 else
     
% convert to 2D matrix
dim = size(features);
if length(dim) > 3 && dim(4) > 1
    features = reshape(features,[],dim(4));
elseif length(dim) > 2 && dim(3) > 1
    features = reshape(features,[],dim(3));
end

% get variables
    num_sub = size(features,2);
    num_features = size(features,1);
% empty matrices
    R_both_true = nan(1, iterations);
    pos_features_true = zeros(num_features, iterations);
    neg_features_true = zeros(num_features, iterations);
    coef_total_true = zeros(num_features, iterations);

    R_both_rand = nan(1, iterations);
% loop through iterations
    for rr = 1:iterations
    fprintf('\n Iteration # %i out of %i\n',rr,iterations);

    [R_both_true(rr), ~, pos_features_true(:,rr), neg_features_true(:,rr),...
        coef_total_true(:,rr)] = funPM(features, behav, model, k,...
        thresh, spearman, s, v_alpha, lambda);

    rand_behav = behav(randperm(num_sub));
    [R_both_rand(rr),~,~,~,~] = funPM(features, rand_behav, model, k,...
        thresh, spearman, s, v_alpha, lambda);
    end
        
    R_both = nanmean(R_both_true);
    R_both_rand(isnan(R_both_rand)) = 0;
    try
        pval_both = (sum(R_both_rand > R_both)+1)/(iterations+1);
    catch
        pval_both = 1;
    end
    pos_features = nansum(pos_features_true > 0,2) >= d*iterations;
    pos_features = pos_features .* nansum(pos_features_true,2);
    if max(pos_features(:)) ~= 0
        pos_features = pos_features./max(pos_features(:));
    end
    neg_features = nansum(neg_features_true > 0,2) >= d*iterations;
    neg_features = neg_features .* nansum(neg_features_true,2);
    if max(neg_features(:)) ~= 0    
        neg_features = neg_features./max(neg_features(:));
    end
    coef_features = nansum(coef_total_true,2);
    if max(coef_features(:)) ~= 0    
        coef_features = coef_features./max(coef_features(:));
    end
    
% reshape features as needed
if length(dim) > 3 && dim(4) > 1
    pos_features = reshape(pos_features,dim(1),dim(2),dim(3));
    neg_features = reshape(neg_features,dim(1),dim(2),dim(3));   
    coef_features = reshape(coef_features,dim(1),dim(2),dim(3));   
elseif length(dim) > 2 && dim(3) > 1
    pos_features = reshape(pos_features,dim(1),dim(2));
    neg_features = reshape(neg_features,dim(1),dim(2));   
    coef_features = reshape(coef_features,dim(1),dim(2));   
end

 end
 
end



