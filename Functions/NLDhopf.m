function fval = NLDhopf(x)
%%	


%% Load parameters
global activation G dsig bfilt afilt C a W N T omega;


%% Setup

% Repopulate connectivity matrix
wC = reshape(x, [N.ROI, N.ROI]);


%% Run nonlinear Hopf model

% Initialize z (x = z(:,1), y = z(:,2))
z = 0.1*ones(N.ROI, 2);

% initialize z (move past initial transients)
BOLD = hopf(z, wC, a, omega, dsig, T, N.ROI);

% Compute connectivity
[~, timeserietotal] = phasesync(BOLD, N.ROI, T.scan, bfilt, afilt);

% compute simulated assembly activations
projection =  W*timeserietotal;


%% Compute KS distance between simulated & empirical assembly activations

% compute KS distance
ksdist = nan(N.assemblies, 1);
for ass = 1:N.assemblies
	[~, ~, ksdist(ass)] = kstest2(activation(ass,:), projection(ass,:));
end

% Compute mean of KS distance
fval = mean(ksdist);
% display(['Current KS Distance: ', num2str(fval)]);
