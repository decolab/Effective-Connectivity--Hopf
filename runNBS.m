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
path{5,1} = fullfile(path{3},'Results','EC');
path{6,1} = fullfile(path{3},'Data','Atlases','AAL');

% Add relevant paths
addpath(genpath(fullfile(path{2}, 'BCT', 'NBS'))),


%% Load and format EC data

% Define files to load
fileName = 'LEICA90_CIC_EC_Iteration2';

% Load EC data
load(fullfile(path{5}, fileName), 'EC','I','N');
ind = logical(I);
clear loadFile

% Format EC data to NxNxM
C = EC; clear EC;
EC = zeros(N.ROI, N.ROI, sum(N.subjects));
for c = 1:N.conditions
	EC(:,:,ind(:,c)) = squeeze(C(:,:,c,ind(:,c)));
	ind(1:N.subjects(c),:) = [];
end
clear C c ind


%% Load & format network labels

% Load network labels
load(fullfile(path{6},'AAL_labels.mat'));

% Convert labels to strings
labels = strings(size(label90,1), 1);
for k = 1:size(label90,1)
	labels(k,1) = string(label90(k,:));
end
clear k label90

% Convert labels to mirrored version
labels = LR_version_symm(labels);


%% Analyze network using directed network-based statistic (NBS)

% Set directed NBS parameters
UI.method.ui = 'Run NBS';	% Options: 'Run NBS' or 'Run FDR'.
UI.test.ui = 't-test';		% Options: 'One Sample', 't-test', 'F-test'.
UI.size.ui = 'Extent';		% Options: 'Extent' or 'Intensity'.  'Intensity' more sensitive to focal effects, 'Extent' to distributed effects.
UI.thresh.ui = '3.1';		% Test statistic threshold.
UI.perms.ui = '5000';		% Number of permutations for computing null distribution.
UI.alpha.ui = '0.05';		% Significance level for p-test.
UI.contrast.ui = '[-1,1]';	% Options: '[-1,1]', '[1,1]', '[1,-1]'
UI.design.ui = I;			% Design matrix (defines which connectome in which group)
UI.matrices.ui = EC;		% Connectome matrices
UI.node_coor.ui = fullfile(path{6},'aal_cog.txt');	% Node coordinates in space
UI.node_label.ui = labels;	% Node labels

% Run NBS
NBSrun(UI);
global nbs;

% Save results
% save(fullfile(path{5}, fileName), 'nbs', '-append');

