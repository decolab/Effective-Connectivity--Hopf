function fval = NLDhopf(x)
%%	


%% Load parameters
global activation dsig bfilt afilt W N T omega co;


%% Setup

% Repopulate connectivity matrix
a = repmat(x(1:N.ROI)', [1,2]);
wC = reshape(x(N.ROI+1:end), [N.ROI, N.ROI]);


%% Run nonlinear Hopf model

% Initialize z (x = z(:,1), y = z(:,2))
z = 0.1*ones(N.ROI, 2);

% initialize z (move past initial transients)
BOLD = hopf(z, wC, a, omega, dsig, T, N.ROI);

% Compute connectivity
[~, timeserietotal, ~, ~] = phasesync(BOLD, N.ROI, T.scan, bfilt, afilt);

% compute simulated assembly activations
projection =  W*timeserietotal;


%% Compute KS distance between simulated & empirical assembly activations

% compute entropy, mean KS distance between assembly distributions
entro = nan(N.assemblies, 2);
% ksdist = nan(N.assemblies, 1);
for ass = 1:N.assemblies
	% Compute entropy of each assembly
	entro(ass, 1) = HShannon_kNN_k_estimation(activation(ass,:), co);
	entro(ass, 2) = HShannon_kNN_k_estimation(projection(ass,:), co);
	
	% Compute KS distance between assembly activation distributions
	%[~, ~, ksdist(ass)] = kstest2(activation(ass,:), projection(ass,:));
end
% fval = mean(ksdist);	% mean KS distance between activtion distributions

% Use KS distance between entropy distributiosn as cost function
[~, ~, fval] = kstest2(entro(:,1), entro(:,2));

% Display cost function
% display(['Current KS Distance: ', num2str(fval)]);