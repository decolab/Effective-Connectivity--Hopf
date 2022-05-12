function fval = NLDhopf(x)
%%	


%% Load parameters
global activation ent dsig bfilt afilt C W N T omega co a aType cfType;


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
		ksdist = nan(N.IC, 1);
		for ass = 1:N.IC
			[~, ~, ksdist(ass)] = kstest2(activation(ass,:), projection(ass,:));
		end
		fval = mean(ksdist);
	case 'entroSum'			% Use difference between sum of entropies as cost function
		entroSim = nan(N.IC, 1);
		for ass = 1:N.IC
			entroSim(ass) = HShannon_kNN_k_estimation(projection(ass,:), co);
		end
		fval = sum([ent entroSim], 1, 'omitnan');
		fval = abs(fval(1) - fval(2));
	case 'eucEntro'			% Use Euclidean distance between total entropies as cost function
		entroSim = nan(N.IC, 1);
		for ass = 1:N.IC
			entroSim(ass) = HShannon_kNN_k_estimation(projection(ass,:), co);
		end
		entroSim = horzcat(ent, entroSim)';
		fval = pdist(entroSim);
	case 'ksEntro'			% Use KS distance between entropy distributions as cost function
		entroSim = nan(N.IC, 1);
		for ass = 1:N.IC
			entroSim(ass) = HShannon_kNN_k_estimation(projection(ass,:), co);
		end
		[~, ~, fval] = kstest2(ent, entroSim);
	case 'maxEntro'			% Use maximum difference between individual elements as cost function
		entroSim = nan(N.IC, 1);
		for ass = 1:N.IC
			entroSim(ass) = HShannon_kNN_k_estimation(projection(ass,:), co);
		end
		fval = max(abs(ent-entroSim));
end