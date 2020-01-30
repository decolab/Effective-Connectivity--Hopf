function [BOLD] = hopf(z, wC, a, omega, dsig, T, nROI)
%% HOPF simulates the Hopf oscillator on a network given by 
%	Simulate Hopf model

% Preallocate storage array for BOLD time series
BOLD = zeros(T.scan, nROI);

% Compute total weight incoming to each node
sumC = repmat(sum(wC,2),1,2);	% for sum Cij*xj

% Set counter
nn = 0;

% Compute first 3000 time steps in order to eliminate transients
for t = 0:T.dt:1000
	suma = wC*z - sumC.*z;	% compute total external influence on node
	zz = z(:,end:-1:1);		% flip z in order to simulate complex conjugate
	dz = a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma;	% Hopf model
	z = z + T.dt*dz + dsig*randn(nROI, 2);		% Euler's method plus random noise
end

% actual modeling (x = BOLD signal (Interpretation), y = some other oscillation)
for t = 0:T.dt:((T.scan-1)*T.TR)

	% Hopf model
	suma = wC*z - sumC.*z;	% remove total influence from node
	zz = z(:,end:-1:1);		% flip z in order to simulate complex conjugate
	dz = a.*z + zz.*omega - z.*(z.*z+zz.*zz) + suma;	% Hopf model
	z = z + T.dt*dz + dsig*randn(nROI,2);	% Euler's method plus random noise

	% Store real part as BOLD signal
	if abs(mod(t,T.TR))<0.01
		nn = nn+1;
		BOLD(nn,:) = z(:,1)';
	end
end

% Convert BOLD signal to standard format
BOLD = BOLD';


end

