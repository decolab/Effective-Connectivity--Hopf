%%	HOPF PARTICLESWARM
% HOPF PARTICLESWARM uses the particle swarm optimization algorithm to fit
% an effective connectivity model to empirical fMRI data.  The effective
% connectivity model is based on the Hopf oscillator.  The particle
% swarm method seeks the connectivity strengths which produce the best
% match between simulated data based on effective connectivity and the
% empirical data.  This opt+imization method is run on each subject
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
fpath{3,1} = fullfile(path{4});
fpath{4,1} = fullfile(path{2},'permutationTest');
fpath{5,1} = fullfile(path{2},'BCT');
fpath{6,1} = fullfile(path{2},'mArrow3');
fpath{7,1} = fullfile(path{2},'spm12');
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
fileName = 'GMR/Control_ICs/LE_COS_ICA_All_wideband_k1_Iteration2';
% fileName = 'Common_ICs/LE_COS_ICA_importIC_ControlvOCD_wideband_k1_Iteration1';

% Set comparisons for visualizations
spaces = "IC";     
ttype = "permutation";      % test on which to base visualization
space = "IC";               % space in which to base visualization: 'ROI' to fit by region, 'ICA' to fit by assembly
typeoffit = 'Subject with Group Prior';
mecDist = 'seuclidean';
cfType = 'eucEntro';	% 'maxEntro';

% Load data
load(fullfile(path{6}, fileName), 'ROI', 'BOLD', 'cortex','origin', 'activities', 'entro','memberships','I','h', 'co','W','T','N','comps','groups','dFC', 'aType', 'SC');

% Set cortex structure
cortex.file = strsplit(cortex.file, filesep);
cortex.file = char(fullfile(path{3},'Atlases','AAL',cortex.file(end)));
if ~isfield(cortex, 'scale') && exist('sphereScale', 'var')
    cortex.scale = sphereScale; clear sphereScale
end
if ~isfield(cortex, 'redux') && exist('rdux', 'var')
    cortex.redux = rdux; clear rdux
end

% Reset N.fig
N.fig = 1;

% File to save
fList = dir(fullfile(path{7}, erase(typeoffit,' '), strcat(fileName, '_EC_', cfType, '_*.mat')));		% Get file list
nIter = numel(fList);										% Find number of previous iterations
fileName = fullfile(erase(typeoffit,' '), strcat(fileName, '_EC_', cfType, '_Iteration', num2str(nIter+1)));
clear nIter fList


%% Compute structural connectivity

% Set connectivity parameters
G = 0.2;				% global coupling weight

% Set structural connectivity
C = SC;
C = C/max(C,[],'all','omitnan')*G;


%% Set parameters

% set color index
cind.node = [1 0 0; 0 0 1];
cind.conn = [1 0 1; 0 1 0];

% Set condition names
if ~exist('groups', 'var')
	groups = I.Properties.VariableNames;
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
			xinit(1,nvars) = C(i, j);
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
			ent = entro.mIC{:,groups{c}};
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
			ent = squeeze(entro.mIC{:,groups{c}});

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
clear x i j s c a ng nn omega ent lb ub initpop options nvars xinit afilt bfilt Isubdiag sig dsig conn activation Cnorm timeseries sc90

% Save results
save(fullfile(path{7}, fileName));

% Display goodness-of-fit metrics
fDims = [0 0 32.46 22.7]; fUnits = 'centimeters';
[F] = Fig5(BOLD, entro, C, EC, N, T, dFC, groups, aType, W, co, fDims, fUnits);
savefig(F, fullfile(path{7}, fileName), 'compact');


%% Display EC of controls, OCD

switch typeoffit
	case {'Subject', 'Subject with Group Prior'}
		
		% Preallocate arrays
		mEC = nan(N.ROI, N.ROI, N.conditions);
		sEC = nan(N.ROI, N.ROI, N.conditions);
		
		% Display mean EC for each condition
		F(numel(F)+1) = figure;
		for c = 1:N.conditions
			
			% Compute mean, standard deviation of EC per group
			mEC(:,:,c) = mean(squeeze(EC(:,:,c, 1:N.subjects(c))), 3, 'omitnan');
			sEC(:,:,c) = std(squeeze(EC(:,:,c, 1:N.subjects(c))), [], 3, 'omitnan');
			
			% Plot mean EC of each condition
			subplot(2, N.conditions, 2*c-1); colormap jet
			imagesc(mEC(:,:,c)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Mean of ', groups{c}]);
			yticks(1:N.ROI); yticklabels(ROI{:,'Label'}); xticks([]);
			pbaspect([1 1 1]);

			% Plot standard deviation EC of each condition
			subplot(2, N.conditions, 2*c); colormap jet
			imagesc(sEC(:,:,c)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Standard Deviation of ', groups{c}]);
			yticks([]); xticks([]); pbaspect([1 1 1]);
		end
		sgtitle('Effective Connectivity');
		clear c

		% Display distance between mean ECs
		F(numel(F)+1) = figure;
		for c = 1:size(comps, 1)
			k(1) = comps(c,1); k(2) = comps(c,2);
			
			% Plot mean
			subplot(size(comps,1), 2, c); colormap jet
			imagesc(pdist2(squeeze(mEC(:,:,k(1))), squeeze(mEC(:,:,k(2))), mecDist)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Mean: ', groups{k(1)}, ' vs. ' groups{k(2)}]);
			yticks(1:N.ROI); yticklabels(ROI{:,'Label'}); xticks([]);
			pbaspect([1 1 1]);
			
			% Plot standard deviation
			subplot(size(comps,1), 2, 2*c); colormap jet
			imagesc(pdist2(squeeze(sEC(:,:,k(1))), squeeze(sEC(:,:,k(2))), mecDist)); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Standard Deviation: ', groups{k(1)}, ' vs. ' groups{k(2)}]);
			yticks([]); xticks([]);
			pbaspect([1 1 1]);
		end
		sgtitle('Distances between Summary Statistics');

		% Display mean distance between ECs
		F(numel(F)+1) = figure;
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
			title(['Mean: ' groups{k(1)} ' vs. ' groups{k(2)}]);
			yticks(1:N.ROI); yticklabels(ROI{:,'Label'}); xticks([]);
			pbaspect([1 1 1]);
			
			% Plot standard deviation of distance
			subplot(size(comps,1), 2, c+1); colormap jet
			imagesc(std(d, 0, 3, 'omitnan')); colorbar;
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			title(['Standard Deviation: ' groups{k(1)} ' vs. ' groups{k(2)}]);
			yticks([]); xticks([]); pbaspect([1 1 1]);
		end
		sgtitle('Summary Statistics of Distances');
		
	case 'Group Mean'
		% Display mean EC for each condition
		F(numel(F)+1) = figure;
		for c = 1:N.conditions
			subplot(1, N.conditions, c);
			xlim([1 N.ROI]); ylim([1 N.ROI]);
			colormap jet
			imagesc(squeeze(EC(:,:,c))); colorbar;
			title(['EC of ', groups{c}, ' Mean Entropy']);
			yticks(1:N.ROI); yticklabels(ROI{:,'Label'});
			pbaspect([1 1 1]);
		end
		clear c

		% Display distance mean differences between ECs in each condition
		F(numel(F)+1) = figure;
		for c = 1:size(comps, 1)
			subplot(1, size(comps,1), c);
			k(1) = comps(c,1); k(2) = comps(c,2);
			mdist = abs(squeeze(EC(:,:,k(1))) - squeeze(EC(:,:,k(2))));	% pdist2(squeeze(EC(:,:,k(1))), squeeze(EC(:,:,k(2))), mecDist);

			xlim([1 N.ROI]); ylim([1 N.ROI]);
			colormap jet
			imagesc(mdist); colorbar;
			title(['Distance between Mean EC of Conditions ', num2str(k(1)), ' and ' num2str(k(2))]);
			yticks(1:N.ROI); yticklabels(ROI{:,'Label'});
			pbaspect([1 1 1]);
		end
end
clear k s d c m
savefig(F, fullfile(path{7}, fileName), 'compact');


%% Strength Analysis

% Set strength index
index = ["In", "Out"];

% Run strength analysis
[strength, vn] = netStrength(EC, I, comps, groups, ROI{:,'Label'}, index, N);

% Visualize strength results
S = netStrengthVis(strength, comps, groups, I, ROI{:,'Label'}, index, N, cind, vn);
clear vn

% Add strength figures to F
for s = 1:numel(S)
	F(numel(F)+1) = S(s);
end
clear s S

% Save strength results
save(fullfile(path{7}, fileName), 'strength', '-append');
savefig(F, fullfile(path{7}, fileName), 'compact');


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
if size(comps,1) > 1
    cont = zeros((size(comps,1)+size(I,2))*2, size(I,2));
    strcont = cell((size(comps,1)+size(I,2))*2, 1);

    % Set all-way contrasts
    c = 2*(size(comps,1)+1:size(comps,1)+N.conditions);
    cont(c-1, :) = -ones(size(I,2), size(I,2)) + diag(2*ones(N.conditions,1));
    cont(c, :) = ones(size(I,2), size(I,2)) - diag(2*ones(N.conditions,1));
    for c = 2*(size(comps,1)+1:size(comps,1)+N.conditions)
        strcont{c-1} = strjoin([groups((c-2*size(comps,1))/2), '> ALL']);
        strcont{c} = strjoin([groups((c-2*size(comps,1))/2), '< ALL']);
    end
else
    cont = zeros(size(comps,1)*2, size(I,2));
    strcont = cell(size(comps,1)*2, 1);
end

% Set pairwise contrasts
for c = 2*(1:size(comps,1))
    cont(c-1, comps(c/2,:)) = [1 -1];
    cont(c, comps(c/2,:)) = [-1 1];
    strcont{c-1} = strjoin([groups(comps(c/2,1)), '>', groups(comps(c/2,2))]);
    strcont{c} = strjoin([groups(comps(c/2,1)), '<', groups(comps(c/2,2))]);
end

% Run NBS analysis
[nbs, STATS, GLM, storarray] = runNBS(EC, cont, I, N, tstat);
nbsNames = nbsNames(nbs, ROI{:,'Full'});      % get names of nodes in significant components

% Save results
save(fullfile(path{7}, fileName), 'h','mEC','sEC','nbs','storarray','nbsNames','strength', '-append');


%% Visualize NBS networks

% Set subplot dimensions
fInds{1} = 1;
fInds{2} = 2;

% Figure dimensions
fDim = [0 256 1240 768; 0 0 1240 512];
fPan = [1 2];

% Find results with paired contrasts
pair = findseq(storarray, 2);

% Find results with single contrasts
sing = findseq(storarray, 1);
sing(ismember(sing, pair, 'rows'),:) = [];

% Visualize paired NBS results
NBS{1} = layered(cortex, nbs(:,unique(pair(:,2))), memberships, h, pair(mod(pair(:,2),2)==1,:), N, cont(unique(pair(:,2)),:), cind, fPan, fDim, fInds, ROI, origin, strcont);
NBS{1} = NBS{1}(isgraphics(NBS{1}));

% Visualize unpaired NBS results
NBS{2} = singlefig(cortex, nbs, memberships, h, sing, N, cont, cind, fPan, fDim, fInds, ROI, origin, strcont);
NBS{2} = NBS{2}(isgraphics(NBS{2}));

% Add NBS figures to F
for s = 1:numel(NBS)
    for n = 1:numel(NBS{s})
        F(numel(F)+1) = NBS{s}(n);
    end
end
clear n s NBS c

% Save figure
savefig(F, fullfile(path{7}, fileName), 'compact');
clear F