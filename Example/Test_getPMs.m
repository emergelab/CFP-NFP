%% FILE SETUP -------------------------------------------------------------
addpath('../getPMs','../getPMs/CircularGraph')
load('atlasindex.mat')
workingatlas = 'A424'; 
aac = 'AA'; modularity = [aac 'C'];
C = atlasindex.(workingatlas).(modularity); 
[nn,~] = size(C);

% Load Data
load('TS.mat')

fprintf('\n Data Loaded\n');    

%% Variables
dCor = 0;
zscore = 0;
levels = [5 20 43 137]; %[5 20 43 137] for AA or [5 19 116] for AAc

cases = TS.cases; ncases = length(cases);
behav = TS.behav;

fprintf('\n Variables Loaded\n');

%% GetNRS LOOP

% loop through all cases
for c = 1:ncases
  
  cname = cases{c};
  fprintf('Running Subject: \t%i\t\t%s\n',c,cname); 
  
  % get NRS
  [niNRS.(cname), neNRS.(cname), tNRS.(cname)] = funaNRS(TS.(cname), C,...
      levels, dCor, zscore);
  
end % end cases loop 

    fprintf('\n GetNRS is Done\n');   
    
    
%% Compile NRS
i = 0;
for ii = levels
    i = i+1;
    Ci = C(:,ii);
    m = max(Ci);
    aa = strcat(aac, num2str(m), '_NRS');
    midx = triu(ones(m)); midx = logical(midx(:));
    b_all_tNRS.(aa) = zeros(m,m,ncases);
    fprintf('Running Level %i %s\n',ii,aa); 
    for c = 1:ncases
    
        bcname = cases{c};
        NRS = zeros(m);
        NRS(midx) = tNRS.(bcname)(1:sum(midx),i);
        b_all_tNRS.(aa)(:,:,c) = NRS;
        
    end
    
    % remove effect of covariates
  if exist('covars','var')
    NRS = reshape(b_all_tNRS.(aa),[],ncases)';
    NRS = NRS(:,midx); res = nan(size(NRS));
    for n = 1:sum(midx)
        mdl = fitlm(covars,NRS(:,n));
        res(:,n) = mdl.Residuals.Raw;
    end
    br_all_tNRS.(aa) = zeros(m*m,ncases);
    br_all_tNRS.(aa)(midx,:) = res';
    br_all_tNRS.(aa) = reshape(br_all_tNRS.(aa),m,m,ncases);
  end
end

    fprintf('\n Compile NRS is Done\n'); 

    
    

%% Run NRS-PM
k = 10;
spearman = 0;
thresh = 0.05;
iterations = 0; % change to 200
model = 'wcpm';
s = 1; d = 0;
v_alpha = 1e-6;
lambda = 1;

nrsR = table; i = 0;
for ii = levels
    i = i+1;
    Ci = C(:,ii);
    m = max(Ci);
    aa = strcat(aac, num2str(m), '_NRS');
    nrsR.Modules(i,1) = m;
    nrsR.Properties.RowNames{i} = aa;
    fprintf('Running Level %i %s\n',ii,aa); 
    
    [nrsR.r(i,1), nrsR.p(i,1), nrsPosEdges.(aa), nrsNegEdges.(aa),...
        nrsCoef.(aa)] = funpPM(b_all_tNRS.(aa), behav, iterations,...
        model, k, thresh, spearman, s, d, v_alpha, lambda);

end

    fprintf('\n NRS-PM is Done\n');
nrsR    
    
%% plot edges
% Set parameters
nm = 'AA150'; % AA150 for AA & AAc126
lv = 43; % level to plot
laa = 'AA50'; % AA7, AA24 or AA50 for AA & AAc6 or AAc22 for AAc
caa = 'col7'; % col7 or col24 for AA & col6 or col22 for AAc
llim = 0; % 0 if NRS, 1 if FC
pct = 100; % e.g. 2.5 will plot top 2.5% of postive and negative edges
bin = 0; % 0 (weighted edges), 1 (remove edges' weight)
sp = -1.49; % -1.45 to -1.55; change if any edge is outside the circle
lw = .5; % lowest line weight


idx = atlasindex.(workingatlas).(['sorted_' nm]);

% get sorting index
sV = atlasindex.(workingatlas).(modularity)(idx,lv);
[sidx,ia,~] = unique(sV,'stable');
lg = length(sidx);
% get sorted labels
sL = atlasindex.(workingatlas).(['labels_' laa])(idx);
lb = sL(ia);
if llim == 1
    lb = atlasindex.(workingatlas).(['lb_' laa '_' nm])(sidx,:);
end

% get colors
sC = atlasindex.(workingatlas).([caa '_' nm])(idx,:);
col = sC(ia,:);

mname = [aac num2str(lg) '_NRS'];

% get matrices
pM = triu(nrsPosEdges.(mname));
pM = pM + transpose(triu(pM,1));
nM = triu(nrsNegEdges.(mname)); 
nM = nM + transpose(triu(nM,1)); 

M = sum(pM,2) - sum(nM,2); % treshold based on nodal strength
tidx = (M >= prctile(M,(100-pct))); pM(~tidx,~tidx) = 0;
pM = triu(pM(sidx,sidx));
%    pM(pM < prctile(pM(pM > 0),pct)) = 0; % threshold based on edges
if (bin == 1)
    pM = double(pM > 0);
end

tidx = (M <= prctile(M,pct)); nM(~tidx,~tidx) = 0;
nM = triu(nM(sidx,sidx));
%    nM(nM < prctile(nM(nM > 0),pct)) = 0;
if (bin == 1)
    nM = double(nM > 0);
end

% generate figures
figure('Position',[350 100 800 400]);
subplot(1,2,1);
circularGraph(nM,'Colormap',col,'Label',lb,'LineWidth',lw,...
    'StartPoint',sp);

subplot(1,2,2);
circularGraph(pM,'Colormap',col,'Label',lb,'LineWidth',lw,...
    'StartPoint',sp);
    
%% Compile niNRS
i = 0;
for ii = levels
    i = i+1;
    Ci = C(:,ii);
    m = max(Ci);
    aa = strcat(aac, num2str(m), '_niNRS');
    b_all_niNRS.(aa) = zeros(nn,ncases);
    fprintf('Running Level %i %s\n',ii,aa); 

    for c = 1:ncases
    
        bcname = cases{c};
        b_all_niNRS.(aa)(:,c) = niNRS.(bcname)(:,i);
        
    end
    
end

    fprintf('\n Compile niNRS is Done\n');
      
%% Run NPM on niNRS
k = 10;
spearman = 0;
thresh = 0.05;
iterations = 0;
model = 'wcpm';
v_alpha = 1e-6;
lambda = 1;
s = 1; d = 1;

ninrsR = table; ninrsPosNodes = table; ninrsNegNodes = table; ninrsCoef = table;
i = 0;
for ii = [5 137]
    i = i+1;
    Ci = C(:,ii);
    m = max(Ci);
    aa = strcat(aac, num2str(m), '_niNRS');
    ninrsR.Modules(i,1) = m;
    ninrsR.Properties.RowNames{i} = aa;
    fprintf('Running Level %i %s\n',i,aa); 
    
    [ninrsR.r(i,1), ninrsR.p(i,1), ninrsPosNodes.(aa),...
        ninrsNegNodes.(aa), ninrsCoef.(aa)] = funpPM(b_all_niNRS.(aa),...
        behav, iterations, model, k, thresh, spearman, s, d, v_alpha, lambda);

end

    fprintf('\n NPM on niNRS is Done\n');
 ninrsR   

%% Compile neNRS
i = 0;
for ii = levels(1:end-1)
    i = i+1;
    Ci = C(:,ii);
    m = max(Ci);
    aa = strcat(aac, num2str(m), '_neNRS');
    b_all_neNRS.(aa) = zeros(nn,ncases);
    fprintf('Running Level %i %s\n',ii,aa); 

    for c = 1:ncases
    
        bcname = cases{c};
        b_all_neNRS.(aa)(:,c) = neNRS.(bcname)(:,i);

    end
    
end

    fprintf('\n Compile neNRS is Done\n');
    
%% Run NPM on neNRS
k = 10;
spearman = 0;
thresh = 0.05;
iterations = 0;
model = 'wcpm';
v_alpha = 1e-6;
lambda = 1;
s = 1; d = 0;

nenrsR = table; nenrsPosNodes = table; nenrsNegNodes = table; nenrsCoef = table;
i = 0;
for ii = 5
    i = i+1;
    Ci = C(:,ii);
    m = max(Ci);
    aa = strcat(aac, num2str(m), '_neNRS');
    nenrsR.Modules(i,1) = m;
    nenrsR.Properties.RowNames{i} = aa;
    fprintf('Running Level %i %s\n',i,aa); 
    
    [nenrsR.r(i,1), nenrsR.p(i,1), nenrsPosNodes.(aa),...
        nenrsNegNodes.(aa), nenrsCoef.(aa)] = funpPM(b_all_neNRS.(aa),...
        behav, iterations, model, k, thresh, spearman, s, d, v_alpha, lambda);

end

    fprintf('\n NPM on neNRS is Done\n');
nenrsR
   