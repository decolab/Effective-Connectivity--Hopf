%%	SETUP

% Clear
clear; close all; clc;


%% Set paths & directories

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
path{5,1} = fullfile(path{3},'Data');
path{6,1} = fullfile(path{3},'Functions');
path{7,1} = fullfile(path{3},'Results','EC');

% Add relevant paths
addpath(genpath(path{6}));


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_CIC_EC';

% Load data
load(fullfile(path{7}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_', S{2}, '_EC');
clear loadFile S


%% Reset paths & directories

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{4,1} = strjoin(path{1}, '/');
path{1,1} = strjoin(path{1}(1:end-2),'/');
path{2,1} = fullfile(path{1},'MATLAB');

% Set required subdirectories
path{5,1} = fullfile(path{3},'Data');
path{6,1} = fullfile(path{4},'Functions');
path{7,1} = fullfile(path{3},'Results','EC');



%% Analyze network using network-based statistic (NBS)

% Add path to NBS
addpath(fullfile(path{2}, 'BCT', 'NBS', 'directed')),

% Load network labels
load(fullfile(path{5}, 'Atlases', 'AAL', 'AAL_labels.mat'));

% Convert labels to strings
labelROI = strings(size(label90,1), 1);
for k = 1:size(label90,1)
	labelROI(k,1) = string(label90(k,:));
end
clear k label90

% Convert labels to mirrored version
labelROI = LR_version_symm(labelROI);

% Convert EC to 3D format
EC3 = nan(N.ROI, N.ROI, sum(N.subjects));
k = 0;
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		k = k+1;
		EC3(:,:,k) = EC(:,:,c,s);
	end
end


%% Analyze network using network-based statistic (NBS)

% Add path to NBS
addpath(fullfile(path{2}, 'BCT', 'NBS', 'directed')),

% Load network labels
load(fullfile(path{5}, 'Atlases', 'AAL', 'AAL_labels.mat'));

% Convert labels to strings
labelROI = strings(size(label90,1), 1);
for k = 1:size(label90,1)
	labelROI(k,1) = string(label90(k,:));
end
clear k label90

% Convert labels to mirrored version
labelROI = LR_version_symm(labelROI);

% Convert EC to 3D format
EC3 = nan(N.ROI, N.ROI, sum(N.subjects));
k = 0;
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		k = k+1;
		EC3(:,:,k) = EC(:,:,c,s);
	end
end
clear c s k

% Set GLM parameters
GLM.perms = 5000;					% number of permutations (scalar)
GLM.X = horzcat(ones(size(I,1),1), I);	% design matrix (observations x conditions)
GLM.contrast = [0, 1, -1];			% set contrast element of intercept term to zero

% Set statistical test to run
conttype = {[0, 1, -1], [1, -1]};		% set contrast element of intercept term to zero
testtype = {'ttest', 'ftest'};			% GLM.test accepts 'ttest' or 'ftest'
stattype = {'Extent', 'Intensity'};		% STATS.size accepts 'Extent' or 'Intensity'
STATS.alpha = 0.05;						% significance threshold

% Preallocate storage array
nbs = cell(numel(1:2*4), numel(conttype), numel(testtype), numel(stattype));

% Run NBS over test types, statistical sizes, and multiple threshold
for c = 1:numel(conttype)
	GLM.contrast = conttype{c};
	GLM.X = GLM.X(:,c:end);
	
	for test = 1:numel(testtype)
		GLM.test = testtype{test};

		for stat = 1:numel(stattype)
			STATS.size = stattype{stat};

			for thresh = 1:2*4
				STATS.thresh = 0.5*thresh;			% test statistic threshold
				nbs{thresh, c, test, stat} = NBSdirected(EC3, GLM, STATS);
			end
		end
	end
end
clear stat test thresh EC3

% Save results
save(fullfile(path{7}, fileName), 'nbs', '-append');


%% Compare group-level degree distributions

% Add paths
addpath(fullfile(path{2}, 'BCT'));

% Preallocate storage arrays
densities = nan(N.conditions, max(N.subjects));
degree.in = nan(N.ROI, N.conditions, max(N.subjects), 2);
degree.out = nan(N.ROI, N.conditions, max(N.subjects), 2);
degree.global.h = nan(3,2);
degree.global.p = nan(3,2);

% Convert h, p to tables
degree.global.h = array2table(degree.global.h, 'VariableNames',{'In','Out'}, 'RowNames',{'Unnormalized', 'Normalized', 'Density'});
degree.global.p = array2table(degree.global.p, 'VariableNames',{'In','Out'}, 'RowNames',{'Unnormalized', 'Normalized', 'Density'});

% Compute densities and strengths
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		% Compute densities
		densities(c,s) = density_dir(EC(:,:,c,s));
		
		% Unnormalized strengths
		degree.in(:,c,s,1) = sum(EC(:,:,c,s), 2);
		degree.out(:,c,s,1) = sum(EC(:,:,c,s), 1);
		
		% Density-normalized degrees
		% degree.in(:,c,s,2) = degree.in(:,c,s,1) ./ densities(c,s);
		% degree.out(:,c,s,2) = degree.out(:,c,s,1) ./ densities(c,s);
	end
end
clear s c

% Compare unnormalized strength distributions
s11 = degree.in(:,1,:,1); s12 = degree.in(:,2,:,1);	% in-strengths
s21 = degree.out(:,1,:,1); s22 = degree.out(:,2,:,1);	% out-strengths
[degree.global.h{'Unnormalized','In'}, degree.global.p{'Unnormalized','In'}] = kstest2(s11(:), s12(:));		% unnormalized in-strengths
[degree.global.h{'Unnormalized','Out'}, degree.global.p{'Unnormalized','Out'}] = kstest2(s21(:), s22(:));	% unnormalized out-strengths

% Compare normalized strength distributions
s11 = degree.in(:,1,:,2); s12 = degree.in(:,2,:,2);	% in-strengths
s21 = degree.out(:,1,:,2); s22 = degree.out(:,2,:,2);	% out-strengths
[degree.global.h{'Normalized','In'}, degree.global.p{'Normalized','In'}] = kstest2(s11(:), s12(:));		% normalized in-strengths
[degree.global.h{'Normalized','Out'}, degree.global.p{'Normalized','Out'}] = kstest2(s21(:), s22(:));	% normalized out-strengths
clear s11 s12 s21 s22

% Compare densities
[degree.global.h{'Density','In'}, degree.global.p{'Density','Out'}] = kstest2(densities(1,:), densities(2,:));

% Visualize densities
D(1) = figure;
histogram(densities(1,:)); hold on;
histogram(densities(2,:));
title('Connection Densities of Hopf-Based Effective Connectivity');
xlabel('Density');
ylabel('Counts');
legend(labels.Properties.VariableNames);
	
% Display strength distributions
D(2) = figure;
for w = 1:2
	subplot(2,2, 2*w-1);
	for c = 1:N.conditions
		histogram(degree.in(:,c,:,w)); hold on;
	end
	title('In-Strength Distribution of Hopf-Based Effective Connectivity');
	xlabel('In-Strength');
	ylabel('Counts');
	legend(labels.Properties.VariableNames);

	subplot(2,2, 2*w);
	for c = 1:N.conditions
		histogram(degree.out(:,c,:,w)); hold on;
	end
	title('Out-Strength Distribution of Hopf-Based Effective Connectivity');
	xlabel('Out-Strength');
	ylabel('Counts');
	legend(labels.Properties.VariableNames);
end
clear s c

% Save figure
savefig(D, fullfile(path{7}, fileName), 'compact');

% Save results
save(fullfile(path{7}, fileName), 'degree', '-append');


%% Compare node-level degree distributions

% Preallocate storage arrays
degree.ROI.p = nan(N.ROI, 2, 2);
degree.ROI.h = nan(N.ROI, 2, 2);

% Test for significant differences between patients and controls
for roi = 1:N.ROI
	% for w = 1:2
		cin = squeeze(degree.in(roi,1,:,w)); cin = cin(~isnan(cin));
		cout = squeeze(degree.out(roi,1,:,w)); cout = cout(~isnan(cout));
		pin = squeeze(degree.in(roi,2,:,w)); pin = pin(~isnan(pin));
		pout = squeeze(degree.out(roi,2,:,w)); pout = pout(~isnan(pout));
		
		[degree.ROI.h(roi, w, 1), degree.ROI.p(roi, w, 1)] = kstest2(cin, pin);
		[degree.ROI.h(roi, w, 2), degree.ROI.p(roi, w, 2)] = kstest2(cout, pout);
	% end
end
clear roi w point cin cout pin pout

% Run FDR multiple comparison correction
degree.ROI.FDR = zeros(N.ROI, 2, 2);
% for w = 1:2
	for point = 1:2
		[ind] = FDR_benjHoch(degree.ROI.p(:,w,point), pval.target);
		degree.ROI.FDR(ind, w, point) = 1;
	end
% end
degree.ROI.FDR = logical(degree.ROI.FDR);
clear ind w point

% Run Bonferroni multiple comparison correction
degree.ROI.Bonferroni = (degree.ROI.p < (pval.target/N.ROI));
degree.ROI.Bonferroni = logical(degree.ROI.Bonferroni);

% Run Dunn-Sidak multiple comparison correction
alpha = 1-(1-pval.target)^(1/N.ROI);
degree.ROI.Sidak = (degree.ROI.p < alpha);
degree.ROI.Sidak = logical(degree.ROI.Sidak);
clear alpha

% Locate significantly different distributions
sig{1,1} = find(degree.ROI.FDR(:,1,1));	% unnormalized, in-strength
sig{1,2} = find(degree.ROI.FDR(:,1,2));	% unnormalized, out-strength
sig{2,1} = find(degree.ROI.FDR(:,2,1));	% normalized, out-strength
sig{2,2} = find(degree.ROI.FDR(:,2,2));	% normalized, out-strength

% Visualize entropies
D(3) = figure;
% Mean Un-Normalized In-Strength
subplot(2,2,1);
p = bar(mean(degree.in(:,:,:,1),3,'omitnan')); hold on;
errorbar((1:N.ROI)-0.15, mean(degree.in(:,1,:,1),3,'omitnan'), std(degree.in(:,1,:,1),0,3,'omitnan'), '.b');
errorbar((1:N.ROI)+0.15, mean(degree.in(:,2,:,1),3,'omitnan'), std(degree.in(:,2,:,1),0,3,'omitnan'), '.r');
scatter(sig{1,1}, 6.2*ones(numel(sig{1,1}),1), '*k');
xticklabels(labelROI);
ylabel('Mean In-Strength');
title('Mean In-Strength Per ROI');
legend(p, labels.Properties.VariableNames);
clear p

% Mean Normalized In-Strength
subplot(2,2,3);
p = bar(mean(degree.in(:,:,:,2),3,'omitnan')); hold on;
errorbar((1:N.ROI)-0.15, mean(degree.in(:,1,:,2),3,'omitnan'), std(degree.in(:,1,:,2),0,3,'omitnan'), '.b');
errorbar((1:N.ROI)+0.15, mean(degree.in(:,2,:,2),3,'omitnan'), std(degree.in(:,2,:,2),0,3,'omitnan'), '.r');
scatter(sig{2,1}, 18*ones(numel(sig{2,1}),1), '*k');
xticklabels(labelROI);
ylabel('Mean Normalized In-Strength');
title('Mean Normalized In-Strength Per ROI');
legend(p, labels.Properties.VariableNames)
clear p

% Mean Out-Strength
subplot(2,2,2);
p = bar(mean(degree.out(:,:,:,1),3,'omitnan')); hold on;
errorbar((1:N.ROI)-0.15, mean(degree.out(:,1,:,1),3,'omitnan'), std(degree.out(:,1,:,1),0,3,'omitnan'), '.b');
errorbar((1:N.ROI)+0.15, mean(degree.out(:,2,:,1),3,'omitnan'), std(degree.out(:,2,:,1),0,3,'omitnan'), '.r');
legend(labels.Properties.VariableNames);
scatter(sig{1,2}, 6.2*ones(numel(sig{1,2}),1), '*k');
xticklabels(labelROI);
ylabel('Mean Out-Strength');
title('Mean Out-Strength Per ROI');
legend(p, labels.Properties.VariableNames);
clear p

% Mean Normalized Out-Strength
subplot(2,2,4);
p = bar(mean(degree.out(:,:,:,2),3,'omitnan')); hold on;
errorbar((1:N.ROI)-0.15, mean(degree.out(:,1,:,2),3,'omitnan'), std(degree.out(:,1,:,2),0,3,'omitnan'), '.b');
errorbar((1:N.ROI)+0.15, mean(degree.out(:,2,:,2),3,'omitnan'), std(degree.out(:,2,:,2),0,3,'omitnan'), '.r');
legend(labels.Properties.VariableNames);
scatter(sig{2,2}, 18*ones(numel(sig{2,2}),1), '*k');
xticklabels(labelROI);
ylabel('Mean Normalized Out-Strength');
title('Mean Normalized Out-Strength Per ROI');
legend(p, labels.Properties.VariableNames);
clear p

% Save figure
savefig(D, fullfile(path{7}, fileName), 'compact');
clear D

% Save results
save(fullfile(path{7}, fileName), 'degree', '-append');


