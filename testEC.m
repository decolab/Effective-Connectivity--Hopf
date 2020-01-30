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
global activation dsig bfilt afilt W N T omega;


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
path{7,1} = fullfile(path{3},'LEICA','Functions');
path{8,1} = fullfile(path{3},'Results','LEICA');

% Add relevant paths
addpath(genpath(pwd));
addpath(genpath(path{6,1}));
addpath(genpath(path{7,1}));


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_CIC_Assemblies';

% Load data
load(fullfile(path{8}, loadFile));

% File to save
S = strsplit(loadFile, '_');
fileName = strcat(S{1}, '_', S{2}, '_EC');
clear loadFile S

% Reset N.fig
N.fig = 1;


%% Reset paths & directories

% Find general path (enclosing folder of current directory)
path{1} = strsplit(pwd, '/');
path{3,1} = strjoin(path{1}(1:end-1),'/');
path{4,1} = strjoin(path{1}, '/');
path{1,1} = strjoin(path{1}(1:end-2),'/');
path{2,1} = fullfile(path{1},'MATLAB');

% Set required subdirectories
path{5,1} = fullfile(path{3},'Data');
path{6,1} = fullfile(path{3},'Results','LEICA');
path{7,1} = fullfile(path{4},'Functions');
path{8,1} = fullfile(path{3},'Results','EC');


%% Compute structural connectivity

% Set connectivity parameters
G = 0.2;				% global coupling weight

% Set structural connectivity (why?)
C = load(fullfile(path{5}, 'sc90.mat'));
C = C.sc90;
C = C/max(max(C))*G;


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

% Set narrowband filter
[bfilt, afilt] = narrowband(T.TR);
params.filt.bfilt = bfilt;
params.filt.afilt = afilt;


%% Compute EC, a per subject

% Preallocate EC matrices
EC = nan(N.ROI, N.ROI, max(N.subjects), N.conditions);	% EC matrices
alpha = nan(N.ROI, max(N.subjects), N.conditions);		% bifurcation parameters
fval = nan(max(N.subjects), N.conditions);				% KS distances

% Generate initial values
a = -0.05*ones(1,N.ROI);
vC = reshape(C, [1, N.ROI^2]);
xinit = horzcat(a, vC);
clear a vC

% Optimize alpha and connectivity in each subject
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		disp(['Computing EC for condition ',num2str(c), ', subject ', num2str(s)])
		
		% Compute omega (intrinsic frequency per node) per condition
		omega = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);
		
		% Set activation
		activation = squeeze(activities.subj.TS{s,c});
		
		% Set optimization parameters
		nvars = N.ROI*(1+N.ROI);
		lb = -G*ones(1, N.ROI^2);		% sets lower connectivity bound to -G
		ub = G*ones(1, N.ROI^2);		% sets upper connectivity bound to G
		initpop = (G/10)*randn(40, nvars)+repmat(xinit,40,1);	% sets initial value(s) for connectivity matrix
		options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop, 'MaxTime',2700000, 'Display','iter');
		
		% Optimize connectivity for each condition
		[x, fval(s,c)] = particleswarm(@NLDhopf, nvars, lb, ub, options);
		display(['Optimized EC fval: ', num2str(fval)]);
		
		% Repopulate connectivity matrix with optimized values
		alpha(:,s,c) = x(1:N.ROI)';
		eC = reshape(x(N.ROI+1:end), [N.ROI, N.ROI]);
		EC(:,:,s,c) = eC;
	end
end
clear x i j s c a ng nn omega lb ub initpop options nvars xinit afilt bfilt Isubdiag sig dsig conn T activation Cnorm timeseries sc90 eC



