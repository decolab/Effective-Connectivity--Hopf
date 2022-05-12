function [F] = Fig5(BOLD, entro, SC, EC, N, T, dFC, groups, aType, W, co, fDims, fUnits)

%% Figure V: Hopf Model Fit
%	This function illustrates the goodness of fit for the Hopf model.
% These measures are at the group level, not individual.  Specific
% visualized measures include:
%	1)	Empirical group-level metrics
%	2)	Group-level metric comparisions per component
%		2.a)  Pre-fit
%		2.b)  Post-fit
%	3)	Joint entropy
%		3.a)  Pre-fit
%		3.b)  Post-fit
%	4)	Distance from simulated to empirical metrics
%	5)	Connectivity matrics
%		5.a)  Strructural atlas
%		5.b)  Estimated connectivity


%% SETUP

if ~isfield(N, 'IC')
    N.IC = N.comp;
end


%% Simulate entropy

% Reshape EC
EC = permute(EC, [4 1 2 3]);

% Define filter for Hopf simulation
z = 0.1*ones(N.ROI, 2);     % x = z(:,1), y = z(:,2)
T.dt = 0.1*T.TR/2;
[bfilt, afilt] = narrowband(T.TR);

% Convert BOLD to array
b = BOLD;
BOLD = nan(N.ROI, T.scan, 3);

% Simulate BOLD, entropy for each subject
entro.ini = nan(N.IC, max(N.subjects), N.conditions);
entro.fit = nan(N.IC, max(N.subjects), N.conditions);
for c = 1:N.conditions
    for s = 1:N.subjects(c)
        
        % empirical BOLD signal
        BOLD(:,:,1) = b{s,c};
        
        % Compute subject-level intrisic frequencies
        omega = findomega_subj(dFC.subj{s,c}, N.ROI, T, afilt, bfilt);
        
        % Simulate oscillator time series & entropy (include role of connectivity)
        BOLD(:,:,2) = hopf(z, SC, -0.0*ones(N.ROI, 2), omega, sqrt(T.dt)*0.02, T, N.ROI);
        BOLD(:,:,3) = hopf(z, squeeze(EC(s,:,:,c)), -0.0*ones(N.ROI, 2), omega, sqrt(T.dt)*0.02, T, N.ROI);
        
        % Compute leading eigenvector dFC
        [~, LEdFC.ini] = phasesync(BOLD(:,:,2), N.ROI, T.scan, bfilt, afilt, aType);
        [~, LEdFC.fit] = phasesync(BOLD(:,:,3), N.ROI, T.scan, bfilt, afilt, aType);
        
        % Compute IC activities
        projection.ini =  W*LEdFC.ini;
        projection.fit =  W*LEdFC.fit;
        
        % Compute simulated entropy pre- and post-fit
        for ic = 1:N.IC
            entro.ini(ic, s, c) = HShannon_kNN_k_estimation(projection.ini(ic,:), co);
            entro.fit(ic, s, c) = HShannon_kNN_k_estimation(projection.fit(ic,:), co);
        end
    end
end
clear s c b ic projection LEdFC BOLD


%% Plot component-level match between empirical and simulated entropies

% Generate labels, color coding
for c = 1:N.conditions
    l(1+(c-1)*N.IC:c*N.IC) = repmat(groups(c),[N.IC 1]);
end
l = string(repmat(l', [2 1]));																					% labels empirical vs. simulated data	%	cat(1, repmat('Control',[N.IC 1]), repmat('Patient',[N.IC 1]))
l(:,2) = strcat(repmat("C ", [2*N.conditions*N.IC 1]), num2str(repmat([1:N.IC]', [2*N.conditions 1])));			% labels component
l(:,3) = string(cat(1, repmat('Empirical',[N.conditions*N.IC 1]), repmat('Simulated',[N.conditions*N.IC 1])));	% labels patient vs. control
lbl = {l(:,2), l(:,1), l(:,3)}; clear l
cg = repmat(1:2*N.conditions, [1 N.IC])';	% just label each box in order it appears on the plot.  Unbelievably stupid: colorGroup should be ordered the same way as the groups

% Open figure
F(1) = figure('Units', fUnits, 'Position', fDims); hold on;

% Plot pre-fit match
e = nan(max(N.subjects), N.IC*N.conditions*2);
for c = 1:N.conditions
    e(:,1+(c-1)*N.IC : c*N.IC) = entro.IC(:,:,c)';
    e(:,1+(c-1+N.conditions)*N.IC : (c+N.conditions)*N.IC) = entro.ini(:,:,c)';
end
% e = horzcat(entro.IC(:,:,1)',entro.IC(:,:,2)', entro.ini(:,:,1)',entro.ini(:,:,2)');
ax(1,1) = subplot(2, 1, 1);
boxplot(ax(1), e, lbl, 'BoxStyle','filled', 'Notch','on', 'ColorGroup',cg, 'LabelVerbosity','minor', 'FactorSeparator',1); hold on;
a = findobj(ax(1), 'Type','text');
b = findall(ax(1), 'Tag','Box');
xlabel('Substate'); ylabel('Entropy');
title('Pre-Fit Entropy', 'FontSize',12);

% Plot post-fit match
e = nan(max(N.subjects), N.IC*N.conditions*2);
for c = 1:N.conditions
    e(:,1+(c-1)*N.IC : c*N.IC) = entro.IC(:,:,c)';
    e(:,1+(c-1+N.conditions)*N.IC : (c+N.conditions)*N.IC) = entro.fit(:,:,c)';
end
% e = horzcat(entro.IC(:,:,1)',entro.IC(:,:,2)', entro.fit(:,:,1)',entro.fit(:,:,2)');
ax(2,1) = subplot(2, 1, 2);
boxplot(ax(2), e, lbl, 'BoxStyle','filled', 'Notch','on', 'ColorGroup',cg, 'LabelVerbosity','minor', 'FactorSeparator',1); hold on;
a(:,2) = findobj(ax(2), 'Type','text');
b(:,2) = findall(ax(2), 'Tag','Box');
xlabel('Substate'); ylabel('Entropy');
title('Post-Fit Entropy', 'FontSize',12);
clear e

% Reshape text and box arrays for ease of manipulation
a = flip(a); b = flip(b);
a = reshape(a, [2*N.conditions*N.IC,3,2]);

% Add legends to boxplots
lgnd = cell(4,1);
for i = 1:2*N.conditions
    lgnd{i} = strcat(string(a(i,2).String), ", ", string(a(i,3).String));
end
legend(ax(1), b(1:2*N.conditions, 1), lgnd, 'Location','eastoutside'); %, 'Location','northoutside', 'Orientation','horizontal');
legend(ax(2), b(1:2*N.conditions, 2), lgnd, 'Location','eastoutside'); %, 'Location','northoutside', 'Orientation','horizontal');

% Edit boxplot labels
set(a, 'HorizontalAlignment', 'left', 'Fontsize',10);
delete(a(:, 2:end, 1:end));
clear i j e ttl lgnd ylbl c s a b


%% Display pre-, post-fit joint entropies for groups

% Generate labels, color coding
% cg = zeros(N.conditions*max(N.subjects), 1);
for c = 1:N.conditions
    l(1+(c-1)*max(N.subjects) : c*max(N.subjects)) = repmat(groups(c),[max(N.subjects) 1]);
%     cg(1+(c-1)*max(N.subjects) : c*max(N.subjects)) = (c-1)*ones(max(N.subjects),1);
end
l = string(repmat(l', [2 1]));																					% labels empirical vs. simulated data	%	cat(1, repmat('Control',[N.IC 1]), repmat('Patient',[N.IC 1]))
% l(:,2) = strcat(repmat("C ", [2*N.conditions*max(N.subjects) 1]), num2str(repmat([1:max(N.subjects)]', [2*N.conditions 1])));			% labels component
l(:,3) = string(cat(1, repmat('Empirical',[N.conditions*max(N.subjects) 1]), repmat('Simulated',[N.conditions*max(N.subjects) 1])));	% labels patient vs. control
lbl = {l(:,1), l(:,3)}; clear l
% cg = repmat(cg, [2 1]);
% cg = repmat(1:2*N.conditions, [1 N.IC])';	% just label each box in order it appears on the plot.  Unbelievably stupid: colorGroup should be ordered the same way as the groups
% lbl = {string(repmat(cat(1, repmat('Control',[max(N.subjects) 1]), repmat('Patient',[max(N.subjects) 1])), [2 1]))	% labels empirical vs. simulated data
% 		string(cat(1, repmat('Empirical',[2*max(N.subjects) 1]), repmat('Simulated',[2*max(N.subjects) 1])))};		% labels patient vs. control		
cg  = vertcat(zeros(N.conditions*max(N.subjects),1), ones(N.conditions*max(N.subjects),1));

% Compare joint entropies for pre-fit simulations
F(2) = figure('Units', fUnits, 'Position', fDims); hold on;
ax(1,2) = subplot(2, 2, 1);
e = nan(N.IC, max(N.subjects)*N.conditions*2);
for c = 1:N.conditions
    e(:, 1+(c-1)*max(N.subjects) : c*max(N.subjects)) = entro.IC(:,:,c);
    e(:, 1+(c-1+N.conditions)*max(N.subjects) : (c+N.conditions)*max(N.subjects)) = entro.ini(:,:,c);
end
e = sum(e);
boxplot(ax(1,2), e, lbl, 'Notch','on', 'ColorGroup',cg, 'Colors','br', 'LabelVerbosity','minor', 'FactorSeparator',1); hold on;
a = findobj(ax(1,2), 'Type','text');
b = findall(ax(1,2), 'Tag','Box');
ylabel('Entropy');
title('Pre-Fit Joint Entropy', 'FontSize',12);

% Compare joint entropies for post-fit simulations
ax(2,2) = subplot(2, 2, 3);
e = nan(N.IC, max(N.subjects)*N.conditions*2);
for c = 1:N.conditions
    e(:, 1+(c-1)*max(N.subjects) : c*max(N.subjects)) = entro.IC(:,:,c);
    e(:, 1+(c-1+N.conditions)*max(N.subjects) : (c+N.conditions)*max(N.subjects)) = entro.fit(:,:,c);
end
e = sum(e);
% e = sum(horzcat(entro.IC(:,:,1),entro.IC(:,:,2), entro.fit(:,:,1),entro.fit(:,:,2)));
boxplot(ax(2,2), e, lbl, 'Notch','on', 'ColorGroup',cg, 'Colors','br', 'LabelVerbosity','minor', 'FactorSeparator',1); hold on;
a(:,2) = findobj(ax(2,2), 'Type','text');
b(:,2) = findall(ax(2,2), 'Tag','Box');
ylabel('Entropy');
title('Post-Fit Joint Entropy', 'FontSize',12);

% Reshape text and box arrays for ease of manipulation
a = flip(a); b = flip(b);
a = reshape(a, [2*N.conditions,2,2]);

% Add legends to boxplots
lgnd = cell(2,1);
for i = 1:2
    lgnd{i} = strcat(string(a(i,2).String));
end
legend(ax(1,2), b(1:2, 1), lgnd, 'Location','southeast'); %, 'Location','northoutside', 'Orientation','horizontal');
legend(ax(2,2), b(1:2, 2), lgnd, 'Location','southeast'); %, 'Location','northoutside', 'Orientation','horizontal');

% Edit boxplot labels
set(a, 'HorizontalAlignment', 'left', 'Fontsize',10);
delete(a(:, 2:end, 1:end));
clear i j e ttl lgnd ylbl c s a b


%% Display pre-, post-fit Euclidean distances for groups

% Calculate Euclidean distances per subject
e = nan(max(N.subjects), 2*N.conditions);
for c = 1:N.conditions
    for s = 1:N.subjects(c)
        e(s, 2*c-1) = pdist(horzcat(entro.IC(:,s,c), entro.ini(:,s,c))');  % pre-fit distances between entropies
        e(s, 2*c) = pdist(horzcat(entro.IC(:,s,c), entro.fit(:,s,c))');  % post-fit distances between entropies
    end
    l(2*c-1:2*c) = repmat(groups(c), [1 2]);
end
lbl = {l repmat(["Pre-Fit" "Post-Fit"], [1 N.conditions])}; clear l
cg = repmat(1:2, [1 N.conditions])';
clear c s

% Plot Euclidean distances per subject for pre-, post-fit simulations
ax(3,2) = subplot(2, 2, [2 4]);
boxplot(ax(3,2), e, lbl, 'Notch','on', 'ColorGroup',cg, 'Colors','br', 'LabelVerbosity','minor', 'FactorSeparator',1); hold on;
a = findobj(ax(3,2), 'Type','text');
b = findall(ax(3,2), 'Tag','Box');
ylabel("Euclidean Distance");
title("Distance: Empirical to Simulated Entropies");

% Reshape text and box arrays for ease of manipulation
a = flip(a); b = flip(b);
a = reshape(a, [2*N.conditions,2]);

% Add legends to boxplots
lgnd = cell(2,1);
for i = 1:2
    lgnd{i} = strcat(string(a(i,2).String));
end
legend(ax(3,2), b(1:2, 1), lgnd, 'Location','northwest'); %, 'Location','northoutside', 'Orientation','horizontal');

% Edit boxplot labels
set(a, 'HorizontalAlignment', 'left', 'Fontsize',10);
delete(a(:, 2:end, 1:end));
clear i j e ttl lgnd ylbl c s a b


%% Plot connectivity matrices

% % Structural atlas
% F(7) = figure;
% imagesc(SC); axis square; colorbar; hold on
% title('Structural Atlas', 'FontSize',12);
% 
% % Effective connectivity
% F(8) = figure;
% imagesc(squeeze(EC(1,:,:,1))); axis square; colorbar; hold on
% title('Estimated EC', 'FontSize',12);


%%

% % Plot pre-fit, post-fit Euclidean distances for components
% e = nan(N.IC, N.conditions, 2);
% ttl = ["Pre-Fit Entropy Match", "Post-Fit Entropy Match"];
% % Calculate Kolmogorov-Smirnov distances per component
% for c = 1:N.conditions
%     for ic = 1:N.IC
%         [~,~,e(ic,c,1)] = kstest2(entro.IC(ic,:,c), entro.ini(ic,:,c));	% pre-fit distances between entropies
%         [~,~,e(ic,c,2)] = kstest2(entro.IC(ic,:,c), entro.fit(ic,:,c));	% post-fit distances between entropies
%     end
% end
% % Plot Kolmogorov-Smirnov distances per subject for pre-, post-fit simulations
% for i = 1:numel(ttl)
%     m = squeeze(median(e,1,'omitnan'));           % mean empirical-simulation distance
%     ax(2+i) = subplot(2,9, (9*i)-[3 2]);
%     bar(e(:,:,i)); title(ttl(i), 'FontSize',12); hold on
%     for c = 1:N.conditions
%         plot(0:max(N.IC), m(c,i)*ones(1, max(N.IC)+1), ['-', col(c)]);
%     end
%     xlabel("Components"); ylabel("Kolmogorov-Smirnov Distance");
%     legend("Control", "Patient", "median");
% end
% clear m
% 
% % Plot pre-fit, post-fit inter-median distances
% % e(:,[1 3], 2) = abs(squeeze(median(entro.IC,2,'omitnan')) - squeeze(median(entro.ini,2,'omitnan')));
% % e(:,[2 4], 2) = abs(squeeze(median(entro.IC,2,'omitnan')) - squeeze(median(entro.fit,2,'omitnan')));
% % e = e(:,[1 3 2 4],:);
% % lgnd = ["Control, Pre-Fit", "Control, Post-Fit", "Patient, Pre-Fit", "Patient, Post-Fit"];
% % ylbl = ["Kolmogorov-Smirnov Distance", "Inter-Median Distance"];

clear e ax ind lgnd i j p