function [strength, vn] = netStrength(EC, I, comps, labels, labels_ROI, index, N)
%UNTITLED2 Summary of this function goes here
%   Detailed explanation goes here


% Compute in-strength and out-strength arrays
for i = 1:numel(index)
    strength.(index(i)).array = squeeze(sum(EC, i, 'omitnan'));
end

% Run comparisions
vn = cell(size(comps,1), 1);
for c = 1:size(comps,1)
    for i = 1:numel(index)
        [strength.(index(i)).h(:,c), strength.(index(i)).p(:,c), strength.(index(i)).tstat(:,c), strength.(index(i)).FDR(:,c), strength.(index(i)).Bonferroni(:,c), strength.(index(i)).Sidak(:,c)] = ...
        robustTests(strength.(index(i)).array(:,I(:,comps(c,1))), strength.(index(i)).array(:,I(:,comps(c,2))), N.ROI, 'p',0.05, 'testtype','permutation');
    end
    vn{c} = [labels{comps(c,1)}, ' v. ', labels{comps(c,2)}];
end
vn = string(vn);

% Run multiple-comparison corrections for strength
if size(comps,1) > 1
    for i = 1:numel(index)
        p_dum = reshape(strength.(index(i)).p, [size(comps,1)*N.ROI, 1]);		% reshape p-value array
        [f, B, S] = mCompCorr([], p_dum, 0.05);									% run mutliple comparison tests
        strength.(index(i)).FDR = reshape(f, [N.ROI, size(comps,1)]);			% False-Discovery Rate
        strength.(index(i)).Bonferroni = reshape(B, [N.ROI, size(comps,1)]);	% Bonferroni threshold
        strength.(index(i)).Sidak = reshape(S, [N.ROI, size(comps,1)]);			% Sidak threshold
    end
end
clear f B S p_dum

% Convert comparisons to table
if size(comps,1) > 1
    for i = 1:numel(index)
        strength.(index(i)).p = array2table(strength.(index(i)).p, 'RowNames',labels_ROI, 'VariableNames',vn);						% p-value
        strength.(index(i)).h = array2table(strength.(index(i)).h, 'RowNames',labels_ROI, 'VariableNames',vn);						% uncorrected hypothesis test
        strength.(index(i)).tstat = array2table(strength.(index(i)).tstat, 'RowNames',labels_ROI, 'VariableNames',vn);				% test statistic (t-statistic or Hodges' G)
        strength.(index(i)).FDR = array2table(strength.(index(i)).FDR, 'RowNames',labels_ROI, 'VariableNames',vn);					% False-Discovery Rate
        strength.(index(i)).Bonferroni = array2table(strength.(index(i)).Bonferroni, 'RowNames',labels_ROI, 'VariableNames',vn);	% Bonferroni threshold
        strength.(index(i)).Sidak = array2table(strength.(index(i)).Sidak, 'RowNames',labels_ROI, 'VariableNames',vn);				% Sidak threshold
    end
end

end

