function [S] = netStrengthVis(strength, comps, labels, I, labels_ROI, index, N, cind, vn)
%UNTITLED3 Summary of this function goes here
%   Detailed explanation goes here


% Visualize strength results for each comparison (if any)
t = 0;
for c = 1:size(comps,1)
    % Tabulate FDR results
    FDR = table(strength.In.FDR{:,vn(c)}, strength.Out.FDR{:,vn(c)}, 'RowNames',labels_ROI, 'VariableNames',{'In','Out'});
    
    % Check if any significant results: if not, do not plot
    if nnz(FDR{:,:}) > 0
        t = t+1;
        
        % Find significant differences
        [r, col] = find(FDR{:,:});
        
        % Preallocate instrength, outstrength indices
        instr = nan(N.ROI, size(I,1));
        outstr = instr;
        cg = [];
        
        % Isolate strength results for boxplot
        for c2 = 1:size(comps,2)
            instr(:, I(:,comps(c,c2))) = squeeze(strength.In.array(:, I(:,comps(c,c2))));
            outstr(:, I(:,comps(c,c2))) = squeeze(strength.Out.array(:, I(:,comps(c,c2))));
            cg = vertcat(cg, repmat(string(labels{comps(c,c2)}), [nnz(I(:,comps(c,c2))), 1]));
        end
        instr(:, ~sum(I(:,comps(c,:)),2)) = [];
        outstr(:, ~sum(I(:,comps(c,:)),2)) = [];
        
        % Format data, group labels for boxplot
        lbl = repmat(labels_ROI, [size(instr,2), 1]);
        cg = repmat(cg, [N.ROI, 1]);
        lbl = {lbl, cg};
        instr = reshape(instr, [numel(instr), 1]);
        outstr = reshape(outstr, [numel(outstr), 1]);

        % Visualize in-strength results
        S(t) = figure('Position', [0 0 1280 1024]);
        ax = subplot(3, nnz(FDR{:,:}), [1 nnz(FDR{:,:})]);
        boxplot(ax, instr, lbl, 'PlotStyle','compact', 'Notch','on', 'ColorGroup',cg, 'Colors',cind.node); hold on;
        title(['In-Strength: ', vn{c}]);
        scatter(r(col==1), (max(instr,[],'all','omitnan')+0.1).*ones(1,nnz(col==1)), 36, 'r', '*');

        % Visualize out-strength results
        ax = subplot(3, nnz(FDR{:,:}), [1+nnz(FDR{:,:}) 2*nnz(FDR{:,:})]);
        boxplot(ax, outstr, lbl, 'PlotStyle','compact', 'Notch','on', 'ColorGroup',cg, 'Colors',cind.node); hold on;
        title(['Out-Strength: ', vn{c}]);
        scatter(r(col==2), (max(outstr,[],'all','omitnan')+0.1).*ones(1,nnz(col==2)), 36, 'r', '*');

        % Visualize significant nodes
%         for n = 1:nnz(FDR{:,:})
%             f = figure; hold on;
%             hg{1} = histogram(strength.(index(col(n))).array(r(n),1:N.subjects(comps(c,2))), 'Normalization','probability');
%             hg{2} = histogram(strength.(index(col(n))).array(r(n),1:N.subjects(comps(c,1))), 'Normalization','probability');
%             sz = min(hg{1}.BinWidth, hg{2}.BinWidth);
%             close(f);
% 
%             figure(S(t));
%             ax = subplot(3, nnz(FDR{:,:}), 2*nnz(FDR{:,:})+n);
%             for c2 = 1:size(comps,2)
%                 histogram(strength.(index(col(n))).array(r(n),comps(c,c2),1:N.subjects(comps(c,c2))), 'Normalization','probability', 'BinWidth',sz, 'FaceAlpha',0.5, 'FaceColor',cind.node(c2,:)); hold on;
%             end
%             legend(labels(comps(c,:)));
%             title(['Modeled ', FDR.Properties.VariableNames{col(n)}, '-Strength of ', labels_ROI{r(n)}]);
%         end
    else
        S = [];
    end
end
clear n ax r c lbl f cg instr outstr i c2 col t index vn


end

