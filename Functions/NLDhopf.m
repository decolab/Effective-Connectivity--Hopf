function fval = NLDhopf(x)
%%	


%% Load parameters
global activation dsig bfilt afilt C W N T omega co a;


%% Setup

% Repopulate nonzero elements of connectivity matrix
nvars = 0;	% number of nonzero connections
wC = C;		% convert 
for i = 1:N.ROI
	for j=1:N.ROI
		if (C(i,j)>0 || j == N.ROI-i+1)
			nvars = nvars+1;
			wC(i,j) = x(nvars);
		end
	end
end
clear i j nvars

% If generating random initial connectivity matrix
% wC = reshape(x, [N.ROI, N.ROI]);	% if not fitting bifurcation parameters
% a = repmat(x(1:N.ROI)', [1,2]);	% if fitting bifurcation parameters
% wC = reshape(x(N.ROI+1:end), [N.ROI, N.ROI]);	% if fitting bifurcation parameters


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
entro = nan(N.IC, 2);
% ksdist = nan(N.IC, 1);
for ass = 1:N.IC
	% Compute entropy of each assembly
	entro(ass, 1) = HShannon_kNN_k_estimation(activation(ass,:), co);
	entro(ass, 2) = HShannon_kNN_k_estimation(projection(ass,:), co);
	
	% Compute KS distance between assembly activation distributions
	%[~, ~, ksdist(ass)] = kstest2(activation(ass,:), projection(ass,:));
end
% fval = mean(ksdist);	% mean KS distance between activtion distributions

% Use difference between sum of entropies as cost function
fval = sum(entro, 1, 'omitnan');
fval = abs(fval(1) - fval(2));

% Use Euclidean distance between total entropies as cost function
% fval = pdist(entro');

% Use KS distance between entropy distributions as cost function
% [~, ~, fval] = kstest2(entro(:,1), entro(:,2));

% Display cost function
% display(['Current KS Distance: ', num2str(fval)]);