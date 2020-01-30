function [omega] = findomega(signaldata, nROI, nSubj, T, afilt, bfilt)
%%	FINDOMEGA computes the average natural frequency of each ROI over a condition
% 


%% Extract power spectra

% Preallocate subject power spectra
PowSpect = nan(floor(T.scan/2), nROI, nSubj);

% Extract frequency
freq = (0:T.scan/2-1)/T.scan*T.TR;

% Compute power spectra per subject, per condition, per ROI
for nsub = 1:nSubj
	for seed = 1:nROI
		x = detrend(demean(signaldata(seed,:,nsub)));
		ts = zscore(filtfilt(bfilt,afilt,x));
		pw = abs(fft(ts));
		PowSpect(:, seed, nsub) = pw(1:floor(T.scan/2)).^2/(T.scan/T.TR);
	end
end

% Extract power per ROI
Power_Areas = squeeze(mean(PowSpect,3));
for seed=1:nROI
	Power_Areas(:,seed) = gaussfilt(freq, Power_Areas(:,seed)', 0.01);
end

% Find frequencies of maximum power for each ROI
[~, index] = max(Power_Areas,[],1);
freq = freq(index);

% Compute natural frequencies of each area in each condition
omega = repmat(2*pi*freq',1,2);
omega(:,1) = -omega(:,1);


