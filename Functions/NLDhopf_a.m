function fval = NLDhopf_a(A)
%%	


%% Load parameters
global EC activation dsig bfilt afilt W N T omega;


%% Run nonlinear Hopf model

% Initialize z (x = z(:,1), y = z(:,2))
z = 0.1*ones(N.ROI, 2);

% initialize z (move past initial transients)
BOLD = hopf(z, EC, A', omega, dsig, T, N.ROI);

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
