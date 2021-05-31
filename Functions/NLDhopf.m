function fval = NLDhopf(x)
%%	


%% Load parameters
global activation entro dsig bfilt afilt C W N T omega co a aType cfType;


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

% Compute dFC
[~, dFC] = phasesync(BOLD, N.ROI, T.scan, bfilt, afilt, aType);

% compute simulated assembly activations
projection =  W*dFC;


%% Compute cost function

% Compute cost function
switch cfType
	case 'ksActive'			% mean KS distance between activation distributions
		ksdist = nan(N.comp, 1);
		for ass = 1:N.comp
			[~, ~, ksdist(ass)] = kstest2(activation(ass,:), projection(ass,:));
		end
		fval = mean(ksdist);
	case 'entroSum'			% Use difference between sum of entropies as cost function
		entroSim = nan(N.comp, 1);
		for ass = 1:N.comp
			entroSim(ass) = HShannon_kNN_k_estimation(projection(ass,:), co);
		end
		fval = sum([entro entroSim], 1, 'omitnan');
		fval = abs(fval(1) - fval(2));
	case 'eucEntro'			% Use Euclidean distance between total entropies as cost function
		entroSim = nan(N.comp, 1);
		for ass = 1:N.comp
			entroSim(ass) = HShannon_kNN_k_estimation(projection(ass,:), co);
		end
		entroSim = horzcat(entro, entroSim)';
		fval = pdist(entroSim);
	case 'ksEntro'			% Use KS distance between entropy distributions as cost function
		entroSim = nan(N.comp, 1);
		for ass = 1:N.comp
			entroSim(ass) = HShannon_kNN_k_estimation(projection(ass,:), co);
		end
		[~, ~, fval] = kstest2(entro, entroSim);
end

% Display cost function
% display(['Current KS Distance: ', num2str(fval)]);