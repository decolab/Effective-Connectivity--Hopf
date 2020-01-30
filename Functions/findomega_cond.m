function [omega] = findomega(TC_AAL, N, T, params)
%%	FINDOMEGA computes the average natural frequency of the offending 
% 

%% Set parameters

% Set filter coefficients
afilt = params.filt.afilt;
bfilt = params.filt.bfilt;


%% Extract power spectra

% Preallocate subject power spectra
PowSpect = nan(floor(T.scan/2), N.ROI, N.subject);

% Extract subject power spectra
for nsub = 1:N.subject

	% Extract frequency, signal
	signaldata = squeeze(TC_AAL{nsub});
	freq = (0:T.scan/2-1)/T.scan*T.TR;

	% Compute power spectra per subject, per condition, per ROI
	for seed=1:N.ROI
		x = detrend(demean(signaldata(seed,:)));
		ts = zscore(filtfilt(bfilt,afilt,x));
		pw = abs(fft(ts));
		PowSpect(:, seed, nsub) = pw(1:floor(T.scan/2)).^2/(T.scan/T.TR);
	end
end

% Extract mean power per ROI per condition
Power_Areas = squeeze(mean(PowSpect,3,'omitnan'));
for seed=1:N.ROI
	Power_Areas(:,seed) = gaussfilt(freq, Power_Areas(:,seed)', 0.01);
end

% Find frequencies of maximum power for each ROI
[~, index] = max(Power_Areas,[],1);
freq = freq(index);

% Compute natural frequencies of each area in each condition
omega = repmat(2*pi*freq',1,2);
omega(:,1) = -omega(:,1);


