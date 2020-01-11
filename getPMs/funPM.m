function [R_both, pval_both, pos_features, neg_features, coef_features] =...
    funPM(features, behav, model, k, thresh, spearman, s,...
    v_alpha, lambda)

% features = input in 2D, 3D, 4D with subjects in last dimension 
% set defaults
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
        
    if ~exist('v_alpha', 'var') || length(v_alpha) ~= 1 
        v_alpha = 1e-6;
    end
    
    if ~exist('lambda', 'var') || length(lambda) ~= 1 
        lambda = '';
    end
    
% convert to 2D matrix
dim = size(features);
if length(dim) > 3 && dim(4) > 1
    features = reshape(features,[],dim(4));
elseif length(dim) > 2 && dim(3) > 1
    features = reshape(features,[],dim(3));
end

% remove NaN
nanidx = isnan(behav);
behav(nanidx) = [];
features(:,nanidx) = [];

% get variables
num_subj = size(features,2);
num_features = size(features,1);
if k > num_subj
    k = num_subj;
end
% get cross validation indices
cvidx = crossvalind('Kfold',num_subj,k);

% ---------------------------------------
pos_features = nan(num_features,k);
neg_features = nan(num_features,k);
r2w = zeros(num_features,k);
coef_features = zeros(num_features,1);
behav_pred_both = nan(num_subj,1);
        
for leftout = 1:k
                        
        % leave out subjects from features and behavior
        testidx = (cvidx == leftout);
        test_features = features(:,testidx);
            
        train_features = features(:,~testidx);
        train_behav = behav(~testidx);
            
        % correlate all features with behavior         
        [r,p] = corr(train_features',train_behav,'rows','complete');
                        
        % set threshold and define masks                    
        pos_mask = (r > 0 & p < thresh);
        neg_mask = (r < 0 & p < thresh);
        mask = (p < thresh);
            
        % update features
        pos_features(:, leftout) = pos_mask;
        neg_features(:, leftout) = neg_mask;
            
        % get features' weight
        r(isnan(r)) = 0;
        r2w(:, leftout) = r.^2 .* (pos_mask + neg_mask);
            
     % build and test
     if strcmp(model,'cpm') || strcmp(model,'wcpm')
        if strcmp(model,'cpm')
        % get sum of all features
        train_sumpos = sum(train_features(pos_mask,:),1)';
        train_sumneg = sum(train_features(neg_mask,:),1)';
        elseif strcmp(model,'wcpm')                
        % get sum of all features
        train_sumpos = sum(train_features(pos_mask,:).*r2w(pos_mask,...
            leftout),1)';
        train_sumneg = sum(train_features(neg_mask,:).*r2w(neg_mask,...
            leftout),1)';          
        end
            
        train_sumboth = train_sumpos-train_sumneg;
           
        % build model on TRAIN subjs
        fit_both = polyfit(train_sumboth, train_behav,1);
            
        % get sum of all features in TEST subjs
        test_sumpos = sum(test_features(pos_mask,:),1)';
        test_sumneg = sum(test_features(neg_mask,:),1)';
        test_sumboth = test_sumpos-test_sumneg;
            
        % run model on TEST subjs  
        behav_pred_both(testidx) = fit_both(1)*test_sumboth + fit_both(2);
            
     elseif strcmp(model,'rcpm')
        % build model on TRAIN subs
        if isempty(lambda) 
            [fit_coef, fit_info] = lasso(train_features(mask, :)',...
                train_behav, 'Alpha',v_alpha, 'CV', 10);
            idxLambda1SE = fit_info.Index1SE;
            coef = fit_coef(:,idxLambda1SE);
            coef0 = fit_info.Intercept(idxLambda1SE);
%            lambda_total(leftout) = fit_info.Lambda(idxLambda1SE);
        else
            [coef, fit_info] = lasso(train_features(mask, :)',...
                train_behav, 'Alpha',v_alpha, 'Lambda', lambda);
            coef0 = fit_info.Intercept;
        end

        % run model on TEST sub with the best lambda parameter
        behav_pred_both(testidx) = test_features( mask, :)'*coef+coef0;
        
        coef_features( mask, leftout) = coef;
%        coef0_total(:, leftout) = coef0;
     end
     
end
                
% weight features
pos_features = mean(r2w,2) .* (nansum(pos_features,2)>=s*k); 
neg_features = mean(r2w,2) .* (nansum(neg_features,2)>=s*k);
coef_features = (pos_features~=0 | neg_features~=0) .*mean(coef_features,2);
if max(pos_features(:)) > 0
    pos_features = pos_features./max(pos_features(:));
end
if max(neg_features(:)) > 0
    neg_features = neg_features./max(neg_features(:));
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

%compare predicted and observed scores
if spearman == 0
   [R_both, pval_both] = corr(behav_pred_both,behav,'rows','complete');
else
   [R_both, pval_both] = corr(behav_pred_both,behav,'rows','complete',...
                'Type','Spearman');
end
        
end

