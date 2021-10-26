%%	HOPF PARTICLESWARM
% HOPF PARTICLESWARM uses the particle swarm optimization algorithm to fit
% an effective connectivity model to empirical fMRI data.  The effective
% connectivity model is based on the Hopf oscillator.  The particle
% swarm method seeks the connectivity strengths which produce the best
% match between simulated data based on effective connectivity and the
% empirical data.  This optimization method is run on each subject
% sequentially, resulting in subject-specific connectivity maps.
% 
% This script may be broken into two distinct sections:
%	1)	Vectorizing the nonzero structural connections
%	2)	Optimizing the effective connectivity vector
%	3)	Repopulating the connectivity matrix


%%	SETUP

% Clear
clear; close all; clc;

% Set global variables
global activation ent dsig bfilt afilt a W N T omega G co aType cfType C;


%% Set paths & directories

% Shuffle random seed.  Necessary in array parallelization to avoid
% repeating same random seed across arrays.
rng('shuffle');

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, "/");
path{3,1} = strjoin(path{1}(1:end-1),"/");
path{4,1} = strjoin(path{1},"/");
path{1,1} = strjoin(path{1}(1:end-2),'/');
path{2,1} = fullfile(path{1},'MATLAB');

% Add relevant paths
fpath{1,1} = fullfile(path{3},'Functions');
fpath{2,1} = fullfile(path{3},'LEICA','Functions');
fpath{3,1} = fullfile(path{3},'EC','Functions');
fpath{4,1} = fullfile(path{4},'Functions');
fpath{5,1} = fullfile(path{2},'permutationTest');
fpath{6,1} = fullfile(path{2},'BCT');
fpath{7,1} = fullfile(path{2},'mArrow3');
fpath{8,1} = fullfile(path{2},'spm12');
for k = 1:numel(fpath)-1
	addpath(genpath(fpath{k}));
end
addpath(fpath{numel(fpath)});
clear fpath k

% Set required subdirectories
path{5,1} = fullfile(path{3},'UCLA','Data');
path{6,1} = fullfile(path{3},'UCLA','Results','LEICA');
path{7,1} = fullfile(path{3},'UCLA','Results','EC');


%% Set file names & load data

% Define file to analyze
fileName = 'GMR/Control_ICs/LE_COS_ICA_All_wideband_k1_Iteration1';

% Set comparisons for visualizations
ttypes = ["kstest2", "permutation"];
spaces = ["dFC", "IC"];     
ttype = "kstest2";          % test on which to base visualization
space = "IC";               % space in which to base visualization: 'ROI' to fit by region, 'ICA' to fit by assembly
typeoffit = 'Subject with Group Prior';
mecDist = 'seuclidean';
cfType = 'eucEntro';

% Load data
load(fullfile(path{6}, fileName), 'cortex','labels_ROI','coords_ROI','origin','sphereScale','rdux', 'activities', 'entro','memberships','I','h', 'co','W','T','N','comps','labels','dFC', 'aType');
N.comp = N.IC; N = rmfield(N, 'IC');
cortex.file = fullfile(path{3},'Atlases','AAL',cortex.file);

% Reset N.fig
N.fig = 1;

% File to save
fList = dir(fullfile(path{7}, strcat(fileName, '_*')));		% Get file list
nIter = numel(fList);										% Find number of previous iterations
if nIter == 0
	nIter = 1;
end
fileName = fullfile(erase(typeoffit,' '), strcat(fileName, '_EC_', cfType, '_Iteration', num2str(nIter)));
clear nIter fList


%% Compute structural connectivity

% Set connectivity parameters
G = 0.2;				% global coupling weight

% Set structural connectivity
d = load(fullfile(path{5}, 'formattedUCLA.mat'));
C = d.ADJ_average; clear d      % Define file to analyze
Cnorm = C/max(C,[],'all','omitnan')*G;


%% Set parameters

% Set condition names
if ~exist('condName', 'var')
	condName = I.Properties.VariableNames;
end

% Temporal parameters
T.dt = 0.1*T.TR/2;

% Noise parameters
dsig = sqrt(T.dt)*0.02;		% set amplitude for Gaussian random noise

% Spatial parameters
a = -0.0*ones(N.ROI, 2);

% Set narrowband filter
[bfilt, afilt] = narrowband(T.TR);
params.filt.bfilt = bfilt;
params.filt.afilt = afilt;

% Set design matrix
I = I{:,:};


%% Compute effective connectivity

% Find number of nonzero connections, approximate structural connectivity
% per condition
nvars = 0;		% number of nonzero connections
for i = 1:N.ROI
	for j=1:N.ROI
		if (C(i,j)>0 || j == N.ROI-i+1)
			nvars = nvars+1;
			xinit(1,nvars) = Cnorm(i, j);
		end
	end
end
clear i j

% Set optimization parameters
lb = zeros(1, nvars);		% sets lower connectivity bound to zero (cannot have negative connectivity)
ub = G*ones(1, nvars);		% sets upper connectivity bound to G
initpop = (G/10)*randn(20,nvars)+repmat(xinit,20,1);	% sets initial value(s) for connectivity matrix
options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop, 'MaxTime',90000, 'Display','iter');


% Switch between fitting methods
switch typeoffit
	case 'Subject'
		% Preallocate EC matrices
		EC = nan(N.ROI, N.ROI, N.conditions, max(N.subjects));	% EC matrices
		alpha = nan(N.ROI, 2, N.conditions, max(N.subjects));	% bifurcation parameters
		fval = nan(max(N.subjects), N.conditions);              % cost function

		% Compute EC per subject
		for c = 1:N.conditions
			for s = 1:N.subjects(c)

				% Display current status
				disp(['Computing EC of condition ', num2str(c),  ', subject ', num2str(s), '.']);

				% Compute omega (intrinsic frequency per node) per condition
				omega = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);

				% Set activation
				ent = squeeze(entro.IC(:,s,c));
				activation = squeeze(activities.subj{s,c});

				% Optimize connectivity for each condition
				[x, fval(s,c)] = particleswarm(@NLDhopf, nvars, lb, ub, options);
				display(['Optimized Particle Swarm: ', num2str(fval(s,c))]);

				% Repopulate connectivity matrix with optimized values
				nn=0;
				for i = 1:N.ROI
					for j = 1:N.ROI
						if (C(i,j)>0 || j == N.ROI-i+1)	% enforces finite values on anti-diagonal
							nn = nn+1;
							EC(i,j,c,s) = x(nn);
						end
					end
				end
			end
		end
		
	case 'Group Mean'
		% Preallocate EC matrices
		EC = zeros(N.ROI, N.ROI, N.conditions);		% EC matrices
		alpha = nan(N.ROI, 2, N.conditions);		% bifurcation parameters
		fval = nan(N.conditions, 1);				% cost function
		
		% Compute mean entropy & activation
		for c = 1:N.conditions
			
			% Extract relevant entropy, activations
			ent = entro.mIC{:,condName{c}};
			activation = squeeze(activities.cond{c});
			
			% Display current status
			disp(['Computing EC of condition ', num2str(c),  '.']);
			
			% Compute omega (intrinsic frequency per node) per condition
			omega = nan(N.ROI, 2, N.subjects(c));
			for s = 1:N.subjects(c)
				omega(:,:,s) = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);
			end
			omega = squeeze(mean(omega, 3, 'omitnan'));

			% Optimize connectivity for each condition
			[x, fval(c)] = particleswarm(@NLDhopf, nvars, lb, ub, options);
			display(['Optimized Particle Swarm: ', num2str(fval(c))]);

			% Repopulate connectivity matrix with optimized values
			nn=0;
			for i = 1:N.ROI
				for j = 1:N.ROI
					if (C(i,j)>0 || j == N.ROI-i+1)	% enforces finite values on anti-diagonal
						nn = nn+1;
						EC(i,j,c) = x(nn);
					end
				end
			end
		end
		clear nn i j s c
		
	case 'Subject with Group Prior'
		% Preallocate EC matrices
		EC = zeros(N.ROI, N.ROI, N.conditions, max(N.subjects));	% EC matrices
		alpha = nan(N.ROI, 2, N.conditions, max(N.subjects));		% bifurcation parameters
		fval = nan(max(N.subjects), N.conditions);		% cost function

		% Compute priors
		for c = 1:N.conditions
			
			% Extract relevant entropy, activations
			ent = entro.mIC{:,condName{c}};
			
			% Display current status
			disp(['Computing prior for condition ', num2str(c),  '.']);
			
			% Compute omega (intrinsic frequency per node) per condition
			omega = nan(N.ROI, 2, N.subjects(c));
			for s = 1:N.subjects(c)
				omega(:,:,s) = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);
			end
			omega = squeeze(mean(omega, 3, 'omitnan'));
			
			% Set activation, entropy
			activation = squeeze(activities.cond{c});
			ent = squeeze(entro.mIC{:,condName{c}});

			% Optimize connectivity for each condition
			[x, fval(c)] = particleswarm(@NLDhopf, nvars, lb, ub, options);
			display(['Optimized Particle Swarm: ', num2str(fval(c))]);
			
			% Find number of nonzero connections, approximate structural connectivity
			% per condition
			nvars = 0;		% number of nonzero connections
			for i = 1:N.ROI
				for j=1:N.ROI
					if (C(i,j)>0 || j == N.ROI-i+1)
						nvars = nvars+1;
						xinit(1,nvars) = x(nvars);
					end
				end
			end
			clear i j

			% Set optimization parameters
			lb = zeros(1, nvars);		% sets lower connectivity bound to zero (cannot have negative connectivity)
			ub = G*ones(1, nvars);		% sets upper connectivity bound to G
			initpop = (G/10)*randn(20,nvars)+repmat(xinit,20,1);	% sets initial value(s) for connectivity matrix
			options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop, 'MaxTime',90000, 'Display','iter');
			
			% Compute EC for each subject
			for s = 1:N.subjects(c)

				% Display current status
				disp(['Computing EC of condition ', num2str(c),  ', subject ', num2str(s), '.']);

				% Compute omega (intrinsic frequency per node) per condition
				omega = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);

				% Set activation
				activation = squeeze(activities.subj{s,c});
				ent = squeeze(entro.IC(:,s,c));

				% Optimize connectivity for each condition
				[x, fval(s,c)] = particleswarm(@NLDhopf, nvars, lb, ub, options);
				display(['Optimized Particle Swarm: ', num2str(fval(s,c))]);

				% Repopulate connectivity matrix with optimized values
				nn=0;
				for i = 1:N.ROI
					for j = 1:N.ROI
						if (C(i,j)>0 || j == N.ROI-i+1)	% enforces finite values on anti-diagonal
							nn = nn+1;
							EC(i,j,c,s) = x(nn);
						end
					end
				end
			end
            EC(:,:,c,N.subjects(c)+1:max(N.subjects)) = nan(N.ROI, N.ROI, max(N.subjects)-N.subjects(c));
		end
end
clear x i j s c a ng nn omega entro lb ub initpop options nvars xinit afilt bfilt Isubdiag sig dsig conn T activation Cnorm timeseries sc90

% Save results
save(fullfile(path{7}, fileName));


%% Display EC of controls, OCD

switch typeoffit
	case {'Subject', 'Subject with Group Prior'}
		
		% Preallocate arrays
		mEC = nan(N.ROI, N.ROI, N.conditions);
		sEC = nan(N.ROI, N.ROI, N.conditions);
		
		% Display mean EC for each condition
		F(1) = figure;
		for c = 1:N.conditions
			
			% Compute mean, standard deviation of EC per group
			mEC(:,:,c) = mean(squeeze(EC(:,:,c, 1:N.subjects(c))), 3, 'omitnan');
			sEC(:,:,c) = std(squeeze(EC(:,:,c, 1:N.subjects(c))), [], 3, 'omitnan');
			
			% Plot mean EC of each condition
			subplot(2, N.conditions, 2*c-1); colormap jet
			imagesc(mEC(:,:,c)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Mean of ', condName{c}]);
			yticks(1:N.ROI); yticklabels(labels_ROI); xticks([]);
			pbaspect([1 1 1]);

			% Plot standard deviation EC of each condition
			subplot(2, N.conditions, 2*c); colormap jet
			imagesc(sEC(:,:,c)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Standard Deviation of ', condName{c}]);
			yticks([]); xticks([]); pbaspect([1 1 1]);
		end
		sgtitle('Effective Connectivity');
		clear c

		% Display distance between mean ECs
		F(2) = figure;
		for c = 1:size(comps, 1)
			k(1) = comps(c,1); k(2) = comps(c,2);
			
			% Plot mean
			subplot(size(comps,1), 2, c); colormap jet
			imagesc(pdist2(squeeze(mEC(:,:,k(1))), squeeze(mEC(:,:,k(2))), mecDist)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Mean: ', condName{k(1)}, ' vs. ' condName{k(2)}]);
			yticks(1:N.ROI); yticklabels(labels_ROI); xticks([]);
			pbaspect([1 1 1]);
			
			% Plot standard deviation
			subplot(size(comps,1), 2, 2*c); colormap jet
			imagesc(pdist2(squeeze(sEC(:,:,k(1))), squeeze(sEC(:,:,k(2))), mecDist)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Standard Deviation: ', condName{k(1)}, ' vs. ' condName{k(2)}]);
			yticks([]); xticks([]);
			pbaspect([1 1 1]);
		end
		sgtitle('Distances between Summary Statistics');

		% Display mean distance between ECs
		F(3) = figure;
		for c = 1:size(comps, 1)
			k(1) = comps(c,1); k(2) = comps(c,2);
			scomb = nchoosek(1:sum(N.subjects(k)), 2);
			scomb(scomb(:,2)<=N.subjects(k(1)),:) = [];
			scomb(scomb(:,1)>N.subjects(k(1)),:) = [];
			scomb(:,2) = scomb(:,2) - N.subjects(k(1));
			d = nan(N.ROI, N.ROI, size(scomb,1));
			
			for s = 1:size(scomb, 1)
				d(:,:,s) = pdist2(squeeze(EC(:,:, scomb(s,1))), squeeze(EC(:,:, scomb(s,2))), mecDist);
			end
			
			% Plot mean of distance
			subplot(size(comps,1), 2, c); colormap jet
			imagesc(mean(d, 3, 'omitnan')); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Mean: ' condName{k(1)} ' vs. ' condName{k(2)}]);
			yticks(1:N.ROI); yticklabels(labels_ROI); xticks([]);
			pbaspect([1 1 1]);
			
			% Plot standard deviation of distance
			subplot(size(comps,1), 2, c+1); colormap jet
			imagesc(std(d, 0, 3, 'omitnan')); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Standard Deviation: ' condName{k(1)} ' vs. ' condName{k(2)}]);
			yticks([]); xticks([]); pbaspect([1 1 1]);
		end
		sgtitle('Summary Statistics of Distances');
		
	case 'Group Mean'
		% Display mean EC for each condition
		F(1) = figure;
		for c = 1:N.conditions
			subplot(1, N.conditions, c);
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			colormap jet
			imagesc(squeeze(EC(:,:,c))); colorbar;
			title(['EC of ', condName{c}, ' Mean Entropy']);
			yticks(1:N.ROI); yticklabels(labels_ROI);
			pbaspect([1 1 1]);
		end
		clear c

		% Display distance mean differences between ECs in each condition
		F(2) = figure;
		for c = 1:size(comps, 1)
			subplot(1, size(comps,1), c);
			k(1) = comps(c,1); k(2) = comps(c,2);
			mdist = abs(squeeze(EC(:,:,k(1))) - squeeze(EC(:,:,k(2))));	% pdist2(squeeze(EC(:,:,k(1))), squeeze(EC(:,:,k(2))), mecDist);

			xlim([1 N.ROI]); ylim([1 N.ROI]);
			colormap jet
			imagesc(mdist); colorbar;
			title(['Distance between Mean EC of Conditions ', num2str(k(1)), ' and ' num2str(k(2))]);
			yticks(1:N.ROI); yticklabels(labels_ROI);
			pbaspect([1 1 1]);
		end
end
clear k s d c m


%% Strength Analysis

% set color index
cind.node = [1 0 0; 0 0 1];
cind.conn = [1 0 1; 0 1 1];

% Set strength index
index = ["In", "Out"];

% Run strength analysis
[strength, vn] = netStrength(EC, I, comps, labels, labels_ROI, index, N);

% Visualize strength results
S = netStrengthVis(strength, comps, labels, I, labels_ROI, index, N, cind, vn);
clear vn

% Save strength results
save(fullfile(path{7}, fileName), 'strength', '-append');


%% NBS analysis

% Analysis Settings
intercept = false;	% determines whether to use intercept term in NBS

% Check if need intercept term
if intercept == true
	I = horzcat(ones(size(I,1),1), I);
end

% Set arrays for parameter sweeps
tstat = 4:0.5:6;

% Allocate contrast matrices
cont = zeros((size(comps,1)+size(I,2))*2, size(I,2));
strcont = cell((size(comps,1)+size(I,2))*2, 1);

% Set all-way contrasts
c = 2*(size(comps,1)+1:size(comps,1)+N.conditions);
cont(c-1, :) = -ones(size(I,2), size(I,2)) + diag(2*ones(N.conditions,1));
cont(c, :) = ones(size(I,2), size(I,2)) - diag(2*ones(N.conditions,1));
for c = 2*(size(comps,1)+1:size(comps,1)+N.conditions)
    strcont{c-1} = strjoin([labels((c-2*size(comps,1))/2), '> ALL']);
    strcont{c} = strjoin([labels((c-2*size(comps,1))/2), '< ALL']);
end

% Set pairwise contrasts
for c = 2*(1:size(comps,1))
    cont(c-1, comps(c/2,:)) = [1 -1];
    cont(c, comps(c/2,:)) = [-1 1];
    strcont{c-1} = strjoin([labels(comps(c/2,1)), '>', labels(comps(c/2,2))]);
    strcont{c} = strjoin([labels(comps(c/2,1)), '<', labels(comps(c/2,2))]);
end

% Run NBS analysis
[nbs, STATS, GLM, storarray] = runNBS(EC, cont, I, N, tstat);
NBS = layered(cortex, nbs, memberships, ttype, h{strcmpi(spaces, space)}(ttype,:), storarray, N, labels, cont, cind, [8 8], coords_ROI, labels_ROI, origin, sphereScale, strcont, rdux);


%% Save results

% Add strength figures to F
for s = 1:numel(S)
	F(numel(F)+1) = S(s);
end
clear s S

% Add NBS figures to F
for n = 1:numel(NBS)
	F(numel(F)+1) = NBS(n);
end
clear n NBS

% Save figure
savefig(F, fullfile(path{7}, typeoffit(~isspace(typeoffit)), fileName), 'compact');
clear F

% Save results
save(fullfile(path{7}, fileName), 'h','mEC','sEC','nbs','strength', '-append');