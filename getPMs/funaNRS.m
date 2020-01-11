function [niNRS,neNRS,tNRS] = funaNRS(TS,C,levels,dCor,zscore)

% TS = time series or correlations matrix
% C = nodal affiliations; % whole-brain (AAC) or cortical (AAcC)
% levels = [20 43]; % [] to run all levels
% dCor = 1; % if correlation matrix
% zscore = 1; % to normalize the connectome

[n,l] = size(C);

if isempty(levels)
    levels = 1:l;
end

l = length(levels);
m = max(C(:,levels(end)));

%calculate Pearson networks
if dCor == 1
    M = TS(1:n,1:n);
else
    TS = TS(:,1:n);
    M = corr(TS);
end
M(1:size(M,1)+1:end) = 0; %zero diag
M(M > 0.99) = 0.99; M(M < -0.99) = -0.99; %clip extreme corr
% remove Infs and NaNs
ni = isnan(M) | isinf(M);
if any(ni(:)); M(ni)=0; end
M = atanh(M);

% for normalizing matrices
if zscore == 1
 mu = mean(M(:)); sigma = std(M(:));
 M = (M - mu) / sigma; 
end

niNRS = nan(n,l); neNRS = nan(n,l); tNRS = nan(sum(sum(triu(ones(m)))),l);
i = 0;
for ii = levels
    Ci = C(:,ii);
    i = i+1;
    m = max(Ci);
 if m == n % FC; faster than the loop
     niNRS(:,i) = mean(M,2); %i.e., nS
     neNRS(:,i) = 0;
     idx = triu(ones(m)); idx = logical(idx(:));
     tNRS(1:sum(idx),i) = M(idx);
 else
    NRS = nan(m);    
    for s = 1:m
        sindx = (Ci==s);
        % nodal internal strength
        niNRS(sindx,i) = mean(M(sindx,sindx),2);
        % nodal external strength
        neNRS(sindx,i) = mean(M(sindx,~sindx),2);
        % internetworks
        for t = 1:m
            tindx = (Ci==t);
            NRS(t,s) = mean(mean(M(tindx,sindx),2));
        end
    end

    idx = triu(ones(m)); idx = logical(idx(:));
    tNRS(1:sum(idx),i) = NRS(idx);
end
end

end

