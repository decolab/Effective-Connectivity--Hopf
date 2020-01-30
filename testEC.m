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
global activation dsig bfilt afilt EC C a W N T omega G c;


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
path{7,1} = fullfile(path{3},'Results','LEICA');
path{8,1} = fullfile(path{3},'Results','EC');

% Add relevant paths
addpath(genpath(path{6}));


%% Set file names & load data

% Define files to load
loadFile = 'LEICA90_CIC_Assemblies';

% Load data
load(fullfile(path{7}, loadFile));

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
Cnorm = C/max(max(C))*G;


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
a = -0.05*ones(N.ROI, 2);

% Set narrowband filter
[bfilt, afilt] = narrowband(T.TR);
params.filt.bfilt = bfilt;
params.filt.afilt = afilt;


%% Compute EC, a per subject

% Preallocate EC matrices
EffC = nan(N.ROI, N.ROI, max(N.subjects), N.conditions);	% EC matrices
alpha = nan(N.ROI, max(N.subjects), N.conditions);			% bifurcation parameters

% Vectorize structural connectivity
xinit = reshape(Cnorm, [1, N.ROI^2]);

% Optimize each subject
for c = 1:N.conditions
	for s = 1:N.subjects(c)
		%% Compute effective connectivity per subject
		disp(['Computing EC for condition ',num2str(c), ', subject ', num2str(s)])
        
		% Compute omega (intrinsic frequency per node) per condition
		omega = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);
		
		% Set activation
		activation = squeeze(activities.subj.TS{s,c});
		
		% Set optimization parameters
		% lb = zeros(1, nvars);		% sets lower connectivity bound to zero (cannot have negative connectivity)
		lb = -G*ones(1, N.ROI^2);		% sets lower connectivity bound to -G
		ub = G*ones(1, N.ROI^2);		% sets upper connectivity bound to G
		initpop = (G/10)*randn(20,N.ROI^2)+repmat(xinit,20,1);	% sets initial value(s) for connectivity matrix
		options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop, 'MaxTime',2700000, 'Display','iter');
		
		% Optimize connectivity for each condition
		[x, fval] = particleswarm(@NLDhopf, N.ROI^2, lb, ub, options);
		display(['Optimized EC fval: ', num2str(fval)]);
		
		% Repopulate connectivity matrix with optimized values
		EC = reshape(x, [N.ROI, N.ROI]);
		EffC(:,:,s,c) = EC;
% 		nn=0;
% 		for i = 1:N.ROI
% 			for j = 1:N.ROI
% 				nn = nn+1;
% 				EC(i,j) = x(nn);
% 			end
% 		end


        %% Compute a per subject

        % Set optimization parameters
        lb = -ones(1, N.ROI);		% sets lower connectivity bound to -1
        ub = ones(1, N.ROI);		% sets upper connectivity bound to 1
        initpop = (G/10)*randn(2*20, N.ROI)+repmat(a',20,1);	% sets initial value(s) for connectivity matrix
        options = optimoptions('particleswarm', 'InitialSwarmMatrix',initpop, 'MaxTime',2700000, 'Display','iter');

        % Optimize bifurcation parameter for each condition
        [alpha(:,s,c), fval] = particleswarm(@NLDhopf_a, N.ROI, lb, ub, options);
        display(['Optimized a fval: ', num2str(fval)]);
	end
end
clear x i j s c a ng nn omega lb ub initpop options nvars xinit afilt bfilt Isubdiag sig dsig conn T activation Cnorm timeseries sc90 EC



