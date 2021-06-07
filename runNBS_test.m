%%	RUNNBS
% RUN NBS analyzes the output of hopf_particleswarm script.  Specifically,
% runNBS applies the directed network-based statistic (NBS) to the
% effective connectivity matrices generated by hopf_particleswarm.
% The NBS utilizes an version of cluster-based statistical methods adapted
% to network topology.  This allows the NBS to drastically reduce the
% multiple comparison problem when comparing links, thus dramatically
% improving its sensitivity to distributed effects.  However, it is no more
% sensitive to focal effects than any other method.


%%	SETUP

%% Set paths & directories

clear; clc; close all;

% Shuffle random seed.  Necessary in array parallelization to avoid
% repeating same random seed across arrays.
rng('shuffle');

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{4,1} = strjoin(path{1}, '/');
path{1,1} = strjoin(path{1}(1:end-2),'/');
path{2,1} = fullfile(path{1},'MATLAB');

% Set required subdirectories
path{5,1} = fullfile(path{3},'OCD','Results','EC');
path{6,1} = fullfile(path{3},'Atlases','AAL');

% Add relevant paths
addpath(genpath(fullfile(path{2}, 'BCT', 'NBS')));
addpath(fullfile(path{2},'spm12'));
addpath(genpath(fullfile(path{2}, 'mArrow3')));
addpath(genpath(fullfile(path{2}, 'permutationTest')));
addpath(genpath(fullfile(path{3}, 'Functions')));
addpath(genpath(fullfile(path{3}, 'LEICA', 'Functions')));


%% Settings

% Define files to load
fileName = 'LEICA90*.mat';
dirName = 'SubjectwithGroupPrior';
fileName = dir(fullfile(path{5}, dirName, fileName));
fileName = fileName.name;

% Analysis Settings
intercept = false;	% determines whether to use intercept term in NBS

% Load network labels
label_ROI = load(fullfile(path{6},'AAL_labels.mat'));
label_ROI = string(label_ROI.label90);			% Convert labels to strings
label_ROI = strip(LR_version_symm(label_ROI));	% Convert labels to symmetric format

% Generate MNI coordinates
coords_ROI = load(fullfile(path{6}, 'aal_cog.txt'), 'aal_cog');
coords_ROI = LR_version_symm(coords_ROI);
origin = [65 45.5 35];							% center origin
MNIscale = 5.5/10;								% scale MNI coordinates
sphereScale = 2.5;								% scale sphere size
coords_ROI = MNIscale*coords_ROI;			% scale MNI coordinates
clear label90 MNIscale

% Set brain rendering parameters
cortex.file = fullfile(path{3},'OCD','Data','MNI152_T1_2mm_brain_mask.nii');	% file containing cortical atlas
cortex.color = [0.9 0.9 0.9];					% color for cortical rendering
cortex.transparency = 0.1;						% set to 1 for opaque cortex
cortex.val = 0.3;								% set isonormal line spacing
cortex.view = [-90 90];							% set camera angle
rdux = 0.7;						% proportion of surface faces to keep (<1)

% set color index
cind.node = [1 0 0; 0 0 1];
cind.conn = cind.node;

% Set number of brain images to generate
render = 'layered';         % Number of images to generate per threshold.
							% 'compare' to compare thresholds in single figure
							% 'single' for single figure per threshold
							% 'multiple' for single figure per network
                            % 'layered' for two-level figure: single figure
                            %           per threshold + individual networks

% Set filename for saved figures
fN = strsplit(fileName, '_');


%% Load and format EC data

% Load EC data
load(fullfile(path{5}, dirName, fileName), 'EC','I','N','condName','combs','mecDist','memberships','h','labels');

% Convert index format
ind = zeros(max(N.subjects), N.conditions);
n = 0;
for c = 1:N.conditions
	ind(1:N.subjects(c),c) = I(n+(1:N.subjects(c)),c);
    n = n + N.subjects(c);
end

% Set formatting indices
ind = logical(ind);
i = logical(I);

% Check if need intercept term
if intercept == true
	I = horzcat(ones(size(I,1),1), I);
end 


%% Run node strength analyses

% Declare strength array
stronk = nan(N.ROI, N.conditions, max(N.subjects), 2);

% Compute in-strength and out-strength
stronk(:,:,:,1) = squeeze(sum(EC, 2));	% In-strength: sum over rows
stronk(:,:,:,2) = squeeze(sum(EC, 1));	% Out-strength: sum over columns

% Run comparisions
strength.in = robustTests(stronk(:,1,1:N.subjects(1),1), stronk(:,2,1:N.subjects(2),1), N.ROI, 'p',0.05, 'testtype','permutation');
strength.out = robustTests(stronk(:,1,:,2), stronk(:,2,:,2), N.ROI, 'p',0.05, 'testtype','permutation');

% Tabulate FDR results
strength.summary = table(strength.in.FDR, strength.out.FDR, 'RowNames',label_ROI, 'VariableNames',{'In','Out'});

% Format labels for strength results
instr = horzcat(squeeze(stronk(:,1,:,1))', squeeze(stronk(:,2,:,1))');
outstr = horzcat(squeeze(stronk(:,1,:,2))', squeeze(stronk(:,2,:,2))');
lbl = {repmat(label_ROI, [2 1]), cat(1, repmat(string(labels.Properties.VariableNames{1}),[N.ROI 1]), repmat(string(labels.Properties.VariableNames{2}),[N.ROI 1]))};
cg = cat(1, repmat(["b";"r"],[N.ROI 1])); %, repmat("b",[N.ROI 1]));	% cg = ["r", "b"];
[r, c] = find(strength.summary{:,:});

% Visualize strength results
S = figure('Position', [0 0 1280 1024]);
ax = subplot(3, nnz(strength.summary{:,:}), [1 nnz(strength.summary{:,:})]);
boxplot(ax, instr, lbl, 'PlotStyle','compact', 'Notch','on', 'ColorGroup',cg); hold on;
title('In-Strength');
scatter(r(c==1), (max(instr,[],'all','omitnan')+0.1).*ones(1,nnz(c==1)), 36, 'r', '*');
ax = subplot(3, nnz(strength.summary{:,:}), [1+nnz(strength.summary{:,:}) 2*nnz(strength.summary{:,:})]);
boxplot(ax, outstr, lbl, 'PlotStyle','compact', 'Notch','on', 'ColorGroup',cg); hold on;
title('Out-Strength');
scatter(r(c==2), (max(outstr,[],'all','omitnan')+0.1).*ones(1,nnz(c==2)), 36, 'r', '*');
for n = 1:nnz(strength.summary{:,:})
	f = figure; hold on;
	hg{1} = histogram(squeeze(stronk(r(n),1,:,c(n))), 'Normalization','probability');
	hg{2} = histogram(squeeze(stronk(r(n),2,:,c(n))), 'Normalization','probability');
	sz = min(hg{1}.BinWidth, hg{2}.BinWidth);
	close(f);
	
	figure(S);
	ax = subplot(3, nnz(strength.summary{:,:}), 2*nnz(strength.summary{:,:})+n);
	histogram(squeeze(stronk(r(n),1,:,c(n))), 'Normalization','Probability', 'BinWidth',sz, 'FaceAlpha',0.5); hold on;
	histogram(squeeze(stronk(r(n),2,:,c(n))), 'Normalization','Probability', 'BinWidth',sz, 'FaceAlpha',0.5);
	legend(labels.Properties.VariableNames);
	title(['Modeled ', strength.summary.Properties.VariableNames{c(n)}, '-Strength of ', label_ROI{r(n)}]);
end
clear n ax r c lbl f cg % stronk instr outstr

% Save strength figure
% savefig(S, fullfile(path{5}, dirName, strcat(fN{1}, '_strength')), 'compact');
% saveas(S, fullfile(path{5}, dirName, strcat(fN{1}, '_strength')), 'png');


%% Set parameter sweep for NBS

% Format EC data to NxNxM
C = EC; clear EC;
EC = zeros(N.ROI, N.ROI, sum(N.subjects));
for c = 1:N.conditions
	EC(:,:,i(:,c)) = squeeze(C(:,:,c,ind(:,c)));
end
clear C c ind n i

% Set arrays for parameter sweeps
tstat = 3.5:0.5:6;
contrast = [-1,1; 1,-1];
strcont = {strjoin([labels.Properties.VariableNames(1), '<', labels.Properties.VariableNames(2)]), strjoin([labels.Properties.VariableNames(1), '>', labels.Properties.VariableNames(2)])};

% Set storage arrays
storarray = nan(length(tstat), size(contrast,1));
nbs = cell(length(tstat), size(contrast,1));


%% Analyze network using directed network-based statistic (NBS)

% Set directed NBS parameters
%	C: connectivity matrices, arranged in a 3D (NxNxM) tensor
%	GLM: structure defining the GLM for analyzing connectivity tensor
GLM.perms = 10000;		% number of permutations (scalar)
GLM.X = I;				% design matrix (observations x conditions)
GLM.test = 'ttest';		% type of statistical test to run ('ttest' or 'ftest')
STATS.size = 'Extent';	% method to measure size of cluster significance ('Intensity' or 'Extent')
STATS.alpha = 0.05;		% significance level to test

% Sweep parameters
for c = 1:size(contrast,1)
	GLM.contrast = contrast(c,:);	% contrast vector (1 x conditions)
	for t = 1:length(tstat)
		STATS.thresh = tstat(t);		% threshold test statistic
		nbs{t, c} = NBSdirected(EC, GLM, STATS);	% Run NBS
		
		storarray(t, c) = ~isempty(nbs{t,c});
	end
end
clear t c


%% Visualize NBS results

% locate significant components
j = h{:,'IC'};
if numel(j) > 1
	i = j{1};
	for k = 2:numel(j)
		i = union(i, j{k});
	end
	i = unique(i);
else
	i = j{:};
end
clear j

% Find average distance between groups
k(1) = combs(1,1); k(2) = combs(1,2);
scomb = nchoosek(1:sum(N.subjects(k)), 2);
scomb(scomb(:,2)<=N.subjects(k(1)),:) = [];
scomb(scomb(:,1)>N.subjects(k(1)),:) = [];
scomb(:,2) = scomb(:,2) - N.subjects(k(1));
d = nan(N.ROI, N.ROI, size(scomb,1));
for s = 1:size(scomb, 1)
	d(:,:,s) = pdist2(squeeze(EC(:,:, scomb(s,1))), squeeze(EC(:,:, scomb(s,2))), mecDist);
end
clear s

[thresh, cont] = find(storarray);
if ~isempty(thresh) && strcmpi(render, 'compare')
	[thresh,~] = unique(thresh);
	F = figure('Position', [0 0 1280 1024]);
	for t = 1:numel(thresh)
		
		% Compute mean distance map
		map = nbs(thresh(t),:);
		for m = 1:numel(map)
			if iscell(map{m})
				for n = 1:numel(map{m})
					map{m}{n} = map{m}{n}.*mean(d, 3, 'omitnan')./10;	% this would be an excellent place to apply recursive programming
				end
			else
				map{m} = map{m}.*mean(d, 3, 'omitnan')./10;
			end
		end
		
		% Plot glass brain
		subplot(3, numel(thresh), t); hold on;
		plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,i),2)), coords_ROI, origin, sphereScale, [], map, cind, strcont, strength.summary, rdux);
		title(['Threshold ', num2str(tstat(t))]);
		
		% Set up subplots
		ax(1) = subplot(3, numel(thresh), t+(numel(thresh)));	% Plot significant connections as binarized connectivity map
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title('NBS Networks');
			yticks([]); xticks([]);
			pbaspect([1 1 1]);
			% Calculate scatter Marker width in points
			currentunits = get(ax(1),'Units');
			set(ax(1), {'Color', 'Units'}, {'k', 'Points'});
			axpos = get(ax(1),'Position');
			set(ax(1), 'Units', currentunits); hold on;
			markerWidth(1) = 1/diff(xlim(ax(1)))*axpos(3);
		ax(2) = subplot(3, numel(thresh), t+2*(numel(thresh)));	% distance map
			colormap(ax(2),cool); hold on
			imagesc(ax(2), mean(d, 3, 'omitnan')); colorbar; hold on
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title('Mean EC Distance');
			yticks([]); xticks([]);
			pbaspect([1 1 1]);
			% Calculate scatter Marker width in points
			currentunits = get(ax(2),'Units');
			set(ax(2), {'Units'}, {'Points'});
			axpos = get(ax(2),'Position');
			set(ax(2), 'Units', currentunits); hold on;
			markerWidth(2) = 1/diff(xlim(ax(2)))*axpos(3);
		
		for m = 1:numel(map)
			
			if iscell(map{m})
				% scale transparency by strength
				for n = 1:numel(map{m})
					a(n) = sum(map{m}{n},'all');
				end
				a = a./max(a,[],'all','omitnan');
				[~,ai(m)] = max(a);
				
				% Find & plot significant connections
				for n = 1:numel(map{m})
					[y, x] = find(map{m}{n});
					scatter(ax(1), x, y, markerWidth(1)^2, [1 1 1].*a(n), 'filled', 's');
				end

				% Highlight mean EC matrix
				for n = 1:numel(map{m})
					[sconns(:,1), sconns(:,2)] = find(nbs{thresh(t),m}{n});	% extract significant connections
					s(n,m) = scatter(ax(2), sconns(:,2), sconns(:,1), markerWidth(2)^2, cind.conn(m,:).*a(n), 's');
					clear sconns
				end
			else
				% Plot significant connections
				[y, x] = find(map{m});
				scatter(ax(1), x, y, markerWidth(1)^2, [1 1 1], 'filled', 's');

				% Highlight mean EC matrix
				[sconns(:,1), sconns(:,2)] = find(nbs{thresh(t),m}{n});	% extract significant connections
				s(1,m) = scatter(ax(2), sconns(:,2), sconns(:,1), markerWidth(2)^2, cind.conn(m,:).*a(n), 's');
				clear sconns
				ai = 1;
			end
		end
		% Plot legend in mean EC matrix
		legend(ax(2), [s(ai(1),1),s(ai(2),2)], strcont, 'Location','bestoutside');
		clear ai
	end
	
	% Save as PNG file, MATLAB figure
	% saveas(F, fullfile(path{5}, dirName, strjoin({fN{1},'NBS'},'_')), 'png');
	% savefig(F, fullfile(path{5}, dirName, strjoin({fN{1},'NBS'},'_')), 'compact');
	
elseif ~isempty(thresh) && strcmpi(render, 'single')
	[thresh,~] = unique(thresh);
	for t = 1:length(thresh)
		
		% Compute mean distance map
		map = nbs(thresh(t),:);
		for m = 1:numel(map)
			if iscell(map{m})
				for n = 1:numel(map{m})
					map{m}{n} = map{m}{n}.*mean(d, 3, 'omitnan')./10;	% this would be an excellent place to apply recursive programming
				end
			else
				map{m} = map{m}.*mean(d, 3, 'omitnan')./10;
			end
		end
		
		% Plot glass brain
		F(t) = figure('Position', [0 0 1280 1024]);
		subplot(2, 4, [3 4 7 8]); hold on;
		plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,i),2)), coords_ROI, origin, sphereScale, [], map, cind, strcont, strength.summary, rdux);
		title(['NBS Networks, Threshold ', num2str(tstat(t))]);
		
		% Set up subplots
		ax(1) = subplot(2, 4, [5 6]); set(ax(1),'Color','k'); hold on;	% Plot significant connections as binarized connectivity map				
			title(['NBS Networks, Threshold ', num2str(tstat(t))]);
			yticks(1:N.ROI); yticklabels(label_ROI); xticks([]);
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			pbaspect([1 1 1]);
			% Calculate scatter Marker width in points
			currentunits = get(ax(1),'Units');
			set(ax(1), {'Color', 'Units'}, {'k', 'Points'});
			axpos = get(ax(1),'Position');
			set(ax(1), 'Units', currentunits); hold on;
			markerWidth(1) = 1/diff(xlim(ax(1)))*axpos(3);
		ax(2) = subplot(2, 4, [1 2]); colormap(ax(2),cool); hold on;	% distance map
			imagesc(ax(2), mean(d, 3, 'omitnan')); colorbar; hold on
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Mean Distance, ', condName{combs(1,1)}, ' to ', condName{combs(1,2)}, ' EC']);
			yticks(1:N.ROI); yticklabels(label_ROI); xticks([]);
			pbaspect([1 1 1]);
			% Calculate scatter Marker width in points
			currentunits = get(ax(2),'Units');
			set(ax(2), {'Units'}, {'Points'});
			axpos = get(ax(2),'Position');
			set(ax(2), 'Units', currentunits); hold on;
			markerWidth(2) = 1/diff(xlim(ax(2)))*axpos(3);
		
		for m = 1:numel(map)
			
			if iscell(map{m})
				% scale transparency by strength
				for n = 1:numel(map{m})
					a(n) = sum(map{m}{n},'all');
				end
				a = a./max(a,[],'all','omitnan');
				
				% Compute significant connections
				for n = 1:numel(map{m})
					[y, x] = find(map{m}{n});
					scatter(ax(1), x, y, markerWidth(1)^2, [1 1 1].*a(n), 'filled', 's');
				end

				% Highlight mean EC matrix
				for n = 1:numel(map{m})
					[sconns(:,1), sconns(:,2)] = find(nbs{thresh(t),m}{n});	% extract significant connections
					s(n,m) = scatter(ax(2), sconns(:,2), sconns(:,1), markerWidth(2)^2, cind.conn(m,:).*a(n), 's');
					clear sconns
				end
			else
				% Plot significant connections
				[y, x] = find(map{m});
				scatter(ax(1), x, y, markerWidth(1)^2, [1 1 1], 'filled', 's');

				% Highlight mean EC matrix
				[sconns(:,1), sconns(:,2)] = find(nbs{thresh(t),m}{n});	% extract significant connections
				s(1,m) = scatter(ax(2), sconns(:,2), sconns(:,1), markerWidth(2)^2, cind.conn(m,:).*a(n), 's');
				clear sconns
			end
		end
		% Plot legend in mean EC matrix
		legend(ax(2), s(1,:), strcont, 'Location','bestoutside', 'Orientation','horizontal');	
		
		% Save as PNG file, MATLAB figure
		% saveas(F(t), fullfile(path{5}, dirName, strcat(fN{1},"_NBS_Threshold", string(join(strsplit(num2str(tstat(t)),'.'),'')))), 'png');
	end
	
	% Save as MATLAB figure
	% savefig(F, fullfile(path{5}, dirName, strcat(fN{1},"_NBS")), 'compact');
	
elseif ~isempty(thresh) && strcmpi(render, 'multiple')
	imgs = cell(numel(unique(thresh)), numel(unique(cont)));
	
	for t = 1:length(thresh)
		ind = [thresh(t), cont(t)];
		for f = 1:numel(nbs{ind(1), ind(2)})
			F(f) = figure('Position', [0 0 1280 1024]);

			% Plot significant connections
			ax(1) = subplot(2, 4, [1 2]); colormap(ax(1),bone); hold on;
			imagesc(full(nbs{ind(1), ind(2)}{f})); colorbar;
			title(['NBS Network ', num2str(f)]);
			yticks(1:N.ROI); yticklabels(label_ROI);
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			pbaspect([1 1 1]);

			% Extract significant connections
			[sconns(:,1), sconns(:,2)] = find(nbs{ind(1), ind(2)}{f});
			
			% Calculate scatter marker width in points
			ax(2) = subplot(2, 4, [5 6]); colormap(ax(2),cool); hold on;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			currentunits = get(ax(2),'Units');
			set(ax(2), {'Units'}, {'Points'});
			axpos = get(ax(2),'Position');
			set(ax(2), 'Units', currentunits); hold on;
			markerWidth(2) = 1/diff(xlim(ax(2)))*axpos(3);

			% Highlight mean EC matrix
			imagesc(mean(d, 3, 'omitnan')); colorbar; hold on
			scatter(sconns(:,2), sconns(:,1), 10, 'g', 's');
			title(['Mean Distance, ', condName{combs(1,1)}, ' to ', condName{combs(1,2)}, ' EC']);
			yticks(1:N.ROI); yticklabels(label_ROI); xticks([]);
			pbaspect([1 1 1]);

			% Render in SPM
			map = nbs{ind(1),ind(2)}{f}.*mean(d, 3, 'omitnan')./10;
			ax = subplot(2, 4, [3 4 7 8]); hold on
			plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,i),2)), coords_ROI, origin, sphereScale, [], map, cind, strcont, strength.summary, rdux);

			% Title figure
			sgtitle(['Threshold ', num2str(tstat(ind(1))), ', ', strcont{ind(2)}]);

			% Save as PNG file
			% saveas(F(f), fullfile(path{5}, dirName, strcat(fileName(1:8), 'T', char(join(strsplit(num2str(tstat(ind(1))),'.'), '')), '_C', char(join(strsplit(num2str(contrast(ind(2),:))), '')), '_N', num2str(f))), 'png');
			clear sconns
		end

		% Place figure(s) in cell array
		imgs{ind(1), ind(2)} = F;

		% Save F
		% savefig(F, fullfile(path{5}, dirName, strcat(fileName(1:8), 'T', char(join(strsplit(num2str(tstat(ind(1))),'.'),'')), '_C', char(join(strsplit(num2str(contrast(ind(2),:))), '')))));
		clear F
	end
	
	% Save as MATLAB figure
	% savefig(F, fullfile(path{5}, dirName, strcat(fN{1},"_NBS")), 'compact');
end
clear K k n f t scomb ax ind nsig s thresh cont m n bin x y F fN


%% Save results

% Save figure
if exist('F', 'var')
	fN = strsplit(fileName, '_');
	% saveas(F, fullfile(path{5}, dirName, strjoin({fN{1},'NBS'},'_')), 'png');
	% savefig(F, fullfile(path{5}, dirName, strjoin({fN{1:end-1},'NBS'},'_')), 'compact');
end
clear c C F S

% Save results
% save(fullfile(path{5}, dirName, fileName(1:end-4)), 'strength','nbs','STATS','GLM','tstat','contrast','storarray', '-append');

