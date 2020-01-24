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

%% Set global variables
global activation dsig bfilt afilt C a W N T omega G;


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
loadFile = 'LEICA90_PhaseAssemblies_Subj_CIC';

% Load data
load(fullfile(path{7}, loadFile));
clear loadFile

% File to save
fileName = 'EC90_CIC';


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
path{8,1} = fullfile(path{3},'Results','EC','Subjects');


%% Compute structural connectivity

% Set connectivity parameters
G = 0.2;				% global coupling weight

% Set structural connectivity (why?)
D = load(fullfile(path{5}, 'sc90.mat'));
C = D.sc90; clear D;
Cnorm = C/max(max(C))*G;


%% Reformat activation, mixing matrices

% % Concatenate ROI activations
% TS = cat(3, TS{1}, TS{2});
% clear TS_AAL
% 
% % Set activation index
% I{:,'Controls'} = zeros(length(I{:,'Controls'}),1);
% I{:,'OCD'} = I{:,'Controls'};
% I{1:size(activations.subj{1}, 3), 'Controls'} = 1;
% I{size(activations.subj{1},3)+1 : size(activations.subj{1},3)+size(activations.subj{2}, 3),'OCD'} = 1;

% Set mixing matrix
W = W.concat;


%% Set parameters

% Determine whether to fit ROIs or ICA components
tofit = 'ICA';  % 'ROI' to fit by region, 'ICA' to fit by assembly

% Temporal parameters
T.TR = 2;				% Time to Repetition (seconds)
T.dt = 0.1*T.TR/2;

% Noise parameters
sig = 0.02;
dsig = sqrt(T.dt)*sig;	% set amplitude for Gaussian random noise

% Spatial parameters
a = -0.0*ones(N.ROI, 2);
Isubdiag = find(tril(ones(N.ROI),-1));	% 

% Set narrowband filter
[bfilt, afilt] = narrowband(T.TR);
params.filt.bfilt = bfilt;
params.filt.afilt = afilt;


%% Compute effective connectivity per subject

% Preallocate EC matrices
EC = zeros(N.ROI, N.ROI, N.condition);
fval = nan(sum(N.subjects),1);

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

% Set design matrix
I = zeros(sum(N.subjects), N.condition);

% Compute EC per condition
for c = 1:N.condition
	
	% Compute EC per subject
	for nsub = 1:N.subjects(c)
		
		% Fill index
		I(sum(N.subjects(1:c-1))+nsub, c) = 1;
		
		% Display current status
		disp(['Computing EC of condition ', num2str(c), ', subject ', num2str(nsub)]);

		% Compute omega (intrinsic frequency per node) per condition
		omega = findomega(squeeze(TS{c}(:,:,nsub)), N.ROI, 1, T, afilt, bfilt);
		
		% Set activation
		activation = squeeze(activations.subj{c}(:,:,nsub));
		
		% Optimize connectivity for each condition
		[x, fval(sum(N.subjects(1:c-1))+nsub)] = particleswarm(@NLDhopf_subj, nvars, lb, ub, options);
		display(['Optimized Particle Swarm: ', num2str(fval(sum(N.subjects(1:c-1))+nsub))]);
		
		% Repopulate connectivity matrix with optimized values
		nn=0;
		for i=1:N.ROI
			for j=1:N.ROI
				if (C(i,j)>0 || j == N.ROI-i+1)
					nn=nn+1;
					EC(i,j, sum(N.subjects(1:c-1))+nsub) = x(nn);
				end
			end
		end
	end
end
I = logical(I);
clear x i j ng nn omega lb ub initpop options nvars xinit a afilt bfilt Isubdiag sig dsig conn T c

% Save results
save(fullfile(path{8}, fileName));


%% Analyze network using network-based statistic (NBS)

% Add NBS folder to path
addpath(genpath(fullfile(path{2}, 'BCT')));

% % Save EC matrices as files
% for nsub = 1:sum(N.subjects)
% 	E = squeeze(EC(:,:,nsub));
% 	save(fullfile(path{8}, [fileName, num2str(nsub)]), 'E');
% end
% clear E c

% Add path to NBS
addpath(fullfile(path{2}, 'BCT', 'NBS')),

% Load network labels
load(fullfile(path{5},'AAL_labels.mat'));

% Convert labels to strings
labels = strings(size(label90,1), 1);
for k = 1:size(label90,1)
	labels(k,1) = string(label90(k,:));
end
clear k

% Convert labels to mirrored version
labels = LR_version_symm(labels);

% Set NBS parameters
UI.method.ui = 'Run NBS';	% Options: 'Run NBS' or 'Run FDR'.
UI.test.ui = 't-test';		% Options: 'One Sample', 't-test', 'F-test'.
UI.size.ui = 'Extent';		% Options: 'Extent' or 'Intensity'.  'Intensity' more sensitive to focal effects, 'Extent' to distributed effects.
UI.thresh.ui = '3.1';		% Test statistic threshold.
UI.perms.ui = '5000';		% Number of permutations for computing null distribution.
UI.alpha.ui = '0.05';		% Significance level for p-test.
UI.contrast.ui = '[1,1]';	% Options: '[-1,1]', '[1,1]', '[1,-1]'
UI.design.ui = I;			% Design matrix (defines which connectome in which group)
UI.matrices.ui= fullfile(path{8}, 'Subjects');	% Directory containing connectome matrices
UI.node_coor.ui = '/Volumes/GoogleDrive/My Drive/OCD/Data/aal_cog.txt';	% Node coordinates in space
UI.node_label.ui = labels;		% Node labels

% Run NBS
NBSrun(UI);
global nbs;

% Save results
save(fullfile(path{8}, fileName), 'nbs', '-append');



