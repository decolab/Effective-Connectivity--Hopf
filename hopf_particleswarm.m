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
global activation dsig bfilt afilt C a W N T omega G co;


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
path{6,1} = fullfile(path{3},'Results','LEICA');
path{7,1} = fullfile(path{3},'Results','EC');
path{8,1} = fullfile(path{3},'Functions');
path{9,1} = fullfile(path{3},'LEICA','Functions');
path{10,1} = fullfile(path{4},'Functions');

% Add relevant paths
addpath(genpath(path{8}));
addpath(genpath(path{9}));
addpath(genpath(path{10}));


%% Set file names & load data

% Define files to load
loadFile = 'MICA90_CIC_EXP_Iteration1';

% Load data
load(fullfile(path{6}, loadFile));

% Reset N.fig
N.fig = 1;


%% Reset paths & directories

clear path;

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{4,1} = strjoin(path{1}, '/');
path{1,1} = strjoin(path{1}(1:end-2),'/');
path{2,1} = fullfile(path{1},'MATLAB');
path{5,1} = fullfile(path{3},'Data');
path{6,1} = fullfile(path{3},'Results','LEICA');
path{7,1} = fullfile(path{3},'Results','EC');

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_', S{2}, '_EC');
fList = dir(fullfile(path{7}, strcat(fileName, '_*')));		% Get file list
nIter = numel(fList);										% Find number of previous iterations
fileName = strcat(fileName, '_Iteration', num2str(nIter));	% Edit fileName
clear loadFile S nIter fList


%% Compute structural connectivity

% Set connectivity parameters
G = 0.2;				% global coupling weight

% Set structural connectivity (why?)
C = load(fullfile(path{5}, 'sc90.mat'));
C = C.sc90;
Cnorm = C/max(C,[],'all','omitnan')*G;


%% Set parameters

% Set condition names
if ~exist('condName', 'var')
	condName = labels.Properties.VariableNames;
end

% Determine whether to fit ROIs or ICA components
tofit = 'ICA';  % 'ROI' to fit by region, 'ICA' to fit by assembly

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


%% Compute effective connectivity per subject

% Preallocate EC matrices
EC = zeros(N.ROI, N.ROI, N.conditions, max(N.subjects));	% EC matrices
alpha = nan(N.ROI, 2, N.conditions, max(N.subjects));		% bifurcation parameters
fval = nan(max(N.subjects), N.conditions);		% cost function

% Set design matrix
I = zeros(sum(N.subjects), N.conditions);
I(1:N.subjects(1), 1) = 1;
I(1+N.subjects(1):sum(N.subjects), 2) = 1;

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

% Compute EC per condition
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		
		% Display current status
		disp(['Computing EC of condition ', num2str(c),  ', subject ', num2str(s), '.']);
		
		% Compute omega (intrinsic frequency per node) per condition
		omega = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);

		% Set activation
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
clear x i j s c a ng nn omega lb ub initpop options nvars xinit afilt bfilt Isubdiag sig dsig conn T activation Cnorm timeseries sc90

% Save results
save(fullfile(path{7}, fileName));


%% Display EC of controls, OCD

% Open figure
F = figure;

% Display EC for each condition
for c = 1:N.conditions
	subplot(1, N.conditions, c);
	xlim([1 N.ROI]);
	ylim([1 N.ROI]);
	imagesc(mean(squeeze(EC(:,:,c,:)), 3)); colorbar;
	title(['Mean EC of ', condName{c}]);
end
clear C

% Save figure
savefig(F, fullfile(path{7}, fileName));
clear F

