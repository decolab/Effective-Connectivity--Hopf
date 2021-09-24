%%	EC_hopf_ICprob
% EC_hopf_ICprob optimizes an effective connectivity model with a gradient
% descent algorithm.  After simulating the expected brain activity, this
% algorithm adjusts each link in the connectivity matrix by a scalar
% multiple of the error between simulation and empirical activation time
% traces.  This continues for a set number of interations or until
% convergence, whichever comes first.  The KS distance between the
% simulated and empirical data are then computed for each potential
% connectivity matrix in order to find the optimal fit.

% This script may be broken into three distinct sections:
%	1)	Extracting the empirical assembly activation time series
%	2)	Computing the effective connectivity using gradient descent
%	3)	Detecting which connectivity matrix and global coupling value
%		provides the optimal fit per condition.


%% Load data

% Load connectivity
load GC_GenCog.mat;
C=GrC;
C=C/max(max(C))*0.2;

% Load metadata (data labels)
load metadata.mat;
Group1=ageMatchedCtlForSZ;

% Load empirical data
load  empiricalICprob_Schizos.mat;	% IC probabilities
load ucla_time_series.mat;			% BOLD time series


%% Set up parameters

% Set data parameters
N.ROI = 82;			% Number of ROIs
T.scan = 152;		% Time of scan
T.TR = 2;			% Time of repetition (seconds)


%% Compute parameters

% Set narrowband filter
[bfilt, afilt] = narrowband(T.TR);
params.filt.bfilt = bfilt;
params.filt.afilt = afilt;

% Compute omega
omega = nan(N.ROI, N.ROI, N.condition);
for ng = 1:N.condition
	N.subject = N.subjects(ng);
	omega(:,:,ng) = findomega(TS_AAL(:,ng), N, T, params);
end
N = rmfield(N, 'subject');

% Compress inputs to NLDhopf
conn.C = C;
params.sim.a = a;
params.sim.omega = omega;
params.sim.sig = sig;


%% Set parameters

% Set model parameters
T.dt = 0.1*T.TR/2;			% 
sig = 0.02;				% 
dsig = sqrt(T.dt)*sig;	% set amplitude of white noise

% set working point of Hopf model (bifurcation parameter)
a = -0.0*ones(N,2);

% set number of iterations
ITER = 1:200;


%% Compute empirical data

% Preallocate storage arrays
FCphases.data = nan(N.ROI, N.ROI, N.condition);
projection.data = cell(N.condition,1);

% % Compute time series phase data & leading eigenvector time series
for ng = 1:N.condition
	
	% preallocate arrays
	timeserietotaldata = cell(1, N.subjects(ng));		% concatenated time series vector
	iPHmean = nan(N.ROI, N.ROI, N.subjects(ng));	% time average of phase connectivity
	
	for nsub = 1:N.subjects(ng)
		signaldata = squeeze(TS_AAL{nsub, ng});
		[~, timeserietotaldata{nsub}, ~, iPHmean(:,:,nsub)] = phasesync(signaldata, N, T.scan, bfilt, afilt);
	end
	
	% Compute FC phases, assembly activation
	timeserietotaldata = cell2mat(timeserietotaldata);	% concatenate subject time series
	FCphases.data(:,:,ng) = squeeze(mean(iPHmean),3);	% find mean instantaneous phase connectivity over subjects
	projection.data{ng} = W.PH * timeserietotaldata;	% recover assembly activations
end

% Remove unnecesary data
clear signaldata ngroup nsub timeserietotaldata iPHmean


%% Find effective connectivity

% set model parameters
Cnew = C;					% connectivity
G = 0:0.02:2;				% coupling parameters to test

% set coupling count variable
iwe = 1;

% Set Hopf model parameters
wC = we*Cnew;	% weighted connectivity
z = 0.1*ones(N,2);		% Hopf variable z; x = z(:,1), y = z(:,2)

% Preallocate data arrays
Phaseserror = nan(numel(G), N.condition);	% mean squared error
Ceff = nan(N, N, numel(G), N.condition);	% effective connectivity per coupling strength

% Obtain effective connectivity, MSE for each coupling strength and each
% condition
for ng = 1:N.condition
	for we = G
		disp(['Testing global coupling ', num2str(we)]);

		% Test each coupling strength many times
		for iter = ITER

			% Simulate Hopf model
			[BOLD.hopf.data] = hopf(T, N, wC, a, z);
			
			% Compute simulated BOLD signal
			[~, ~, iPH.sim, ~] = phasesync(BOLD.hopf.data, N, T, bfilt, afilt);
			
			% Compute mean instantaneous phase connectivity matrix of simulated data
			FCphases.sim = squeeze(mean(iPH.sim),3);
			
			% Adjust connectivity to minimize gap between simulated, empirical
			% data (gradient descent?)
			for i=1:N
				for j=i+1:N
					if (C(i,j)>0 || j==N/2+i)
						
						% Minimze gap between empirical, simulated data
						Cnew(i,j) = Cnew(i,j) + 0.01*(FCphases.data(i,j,ng)-FCphases.sim(i,j));
						
						% Impose non-negativity
						if Cnew(i,j)<0
							Cnew(i,j)=0;
						end
						
						% Impose symmetry (why?)
						Cnew(j,i) = Cnew(i,j);
					end
				end
			end
			
			% Compute normalized connectivity, error scores
			Cnew = Cnew/max(max(Cnew))*0.2;					% normalize connectivity
			D = abs(FCphases.data(:,:,ng)-FCphases.sim).^2;	% compute squared error
			MSE = sum(D(:))/numel(FCphases.sim);			% mean squared error
			if MSE<0.01
				break;
			end
		end
		
		% Record result of each coupling strength
		Phaseserror(iwe, ng) = MSE;	% mean squared error
		Ceff(:,:,iwe,ng) = Cnew;	% effective connectivity
	end
end
clear FCphases D MSE Cnew


%% Find KS distance between simulated and empirical connectivity

% Preallocation storage arrays for KS distance
KSfit = nan(numel(G), N.condition);					% KS distance
ksdist = nan(numel(N.assemblies.ph), 1);	% KS distance per assembly

% set coupling count variable
iwe = 1;

% Set Hopf simulation parameters
z = 0.1*ones(N.ROI, 2); % --> x = z(:,1), y = z(:,2)

% Detect coupling with minimal KS distance per condition
for ng = 1:N.condition
	for we = G
		
		% Set Hopf simulation parameters
		wC = we*Ceff(:,:,iwe,ng);
		
		% Preallocate timeseries cell array
		timeserietotal = cell(1,N.subjects(ng));
		
		% Compute subjectwise 
		for nsub = 1:N.subjects(ng)
			% 
			[BOLD.hopf.sim] = hopf(z, a, wC, omega, dsig, T, N);
			T.sim = size(BOLD.phase.sim,2);
			% Compute total phase synchrony time series
			[~, timeserietotal{nsub}, ~, ~] = phasesync(BOLD.hopf.sim, N, T.sim, bfilt, afilt);
		end
		
		% Convert time series to array form
		timeserietotal = cell2mat(timeserietotal);
		
		% Compute assembly activations from mixing matrix
		projection.sim =  W.PH * timeserietotal;
		
		% Compute KS distance between 
		for ass = 1:N.assemblies.ph
			[~, ~, ksdist(ass)] = kstest2(projection.data{ng}(ass,:), projection.sim(ass,:));
		end
		KSfit(iwe, ng) = mean(ksdist);

		% update count
		iwe = iwe+1;
	end
end
clear iwe ng ksdist timeserietotal wC z t_all

save generative_phases.mat Ceff G omega sig;



