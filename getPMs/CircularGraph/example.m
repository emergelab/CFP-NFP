%% Circular Graph Examples
% Copyright 2016 The MathWorks, Inc.

%% 1. Adjacency matrix of 1s and 0s
% Create an example adjacency matrix made up of ones and zeros.
rng(0);
x = rand(50);
thresh = 0.93;
x(x >  thresh) = 1;
x(x <= thresh) = 0;

%%
% Call CIRCULARGRAPH with only the adjacency matrix as an argument.
circularGraph(x);

%%
% Click on a node to make the connections that emanate from it more visible
% or less visible. Click on the 'Show All' button to make all nodes and
% their connections visible. Click on the 'Hide All' button to make all
% nodes and their connections less visible.

%% 2. Supply custom properties
% Create an example adjacency matrix made up of various values and supply
% custom properties.
rng(0);
x = rand(20);
thresh = 0.93;
x(x >  thresh) = 1;
x(x <= thresh) = 0;
for i = 1:numel(x)
  if x(i) > 0
    x(i) = rand(1,1);
  end
end

%%
% Create custom node labels
myLabel = cell(length(x));
for i = 1:length(x)
  myLabel{i} = num2str(round(1000000*rand(1,1)));
end

%%
% Create custom colormap
figure;
myColorMap = lines(length(x));

circularGraph(x,'Colormap',myColorMap,'Label',myLabel);


%% create sorted labels
sidx = atlasindex.MMP_BNAC.sortAA6BNAC;

gsidx = atlasindex.MMP_BNAC.sorted471AA6;
lg = length(gsidx);

lbAA6 = atlasindex.MMP_BNAC.s424AA6;
glbAA6 = atlasindex.MMP_BNAC.s471AA6;

% create colormap
AA6col = atlasindex.MMP_BNAC.AA6RGB(sidx,:);

gAA6col = atlasindex.MMP_BNAC.gAA6RGB(gsidx,:);

% sort matrix
%sm = M(sidx,sidx);

gsm = M; gsm(lg,lg) = 0; gsm = gsm(gsidx,gsidx);

M = squeeze(corrPosEdges.p0_05);

tidx = (sum(M,2) > prctile(sum(M,2),99));
M(~tidx,~tidx) = 0;
gsm = M; gsm(lg,lg) = 0; gsm = gsm(gsidx,gsidx);

% generate figure
%circularGraph(sm,'Colormap',AA6col,'Label',lbAA6);

circularGraph(gsm,'Colormap',gAA6col,'Label',glbAA6);



%% create sorted labels
sidx = atlasindex.MMP_BNAC.sortedAA22;

%gsidx = atlasindex.MMP_BNAC.sorted471AA6;
%lg = length(gsidx);

lbAA22 = atlasindex.MMP_BNAC.s22AA22;
%glbAA6 = atlasindex.MMP_BNAC.s471AA6;

% create colormap
AA22col = atlasindex.MMP_BNAC.AA22RGB(sidx,:);

%gAA6col = atlasindex.MMP_BNAC.gAA6RGB(gsidx,:);

% sort matrix
sm = M(sidx,sidx);

%gsm = M; gsm(lg,lg) = 0; gsm = gsm(gsidx,gsidx);

%M = squeeze(corrPosEdges.p0_001);

%tidx = (sum(M,2) > prctile(sum(M,2),99));
%M(~tidx,~tidx) = 0;
%gsm = M; gsm(lg,lg) = 0; gsm = gsm(gsidx,gsidx);

% generate figure
circularGraph(sm,'Colormap',AA22col,'Label',lbAA22,'LineWidth',1.5);

%circularGraph(gsm,'Colormap',gAA6col,'Label',glbAA6);
