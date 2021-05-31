

%% Extract individual networks

[ind(:,1), ind(:,2)] = find(storarray);
fN = strsplit(fileName, '_');

for t = 1:max(unique((ind)))
    
    % Open figure
    F(t,1) = figure('Position', [0 0 1240 1024]);
    
    % Plot significant connections as binarized connectivity map
    ax(2,1) = subplot(1, 5, 3:5); hold on;
    yticks(1:N.ROI); xticks([]);
    set(ax(2,1),{'Color', 'FontSize', 'YTickLabel'},{'k', 6, label_AAL90});
    xlim([0.5 N.ROI+0.5]); ylim([0.6 N.ROI+0.5]);
    title("Significant Connections", 'FontSize',12);
    pbaspect([1 1 1]);
        
%     % Plot mean distance between patient & control
%     ax(3,1) = subplot(2, 2, 4); colormap(ax(3,1),cool);
%     imagesc(ax(3,1), mean(d, 3, 'omitnan')); colorbar; hold on
%     yticks(1:N.ROI); xticks([]);
%     set(ax(3,1),{'Color', 'FontSize', 'YTickLabel'},{'k', 4, label_AAL90});
%     xlim([1 N.ROI+0.5]); ylim([0.5 N.ROI]);
%     title(['Mean Distance, ', condName{combs(1,1)}, ' to ', condName{combs(1,2)}, ' EC'], 'FontSize',12);
%     pbaspect([1 1 1]);
    
    % Calculate scatter Marker width in points
    currentunits = get(ax(2,1),'Units');
    set(ax(2,1), {'Units'}, {'Points'});
	axpos = get(ax(2,1),'Position');
	set(ax(2,1), 'Units', currentunits); hold on;
	markerWidth = 1/diff(xlim(ax(2,1)))*axpos(3)-1;
    
    % Compute mean distance map & get number of components per threshold
    map = nbs(ind(t,1),:); s = scatter([],[]); lgnd  = 0;
    for c = 1:numel(map)
        if iscell(map{c})
            for n = 1:numel(map{c})
                map{c}{n} = map{c}{n}.*mean(d, 3, 'omitnan')./10;	% this would be an excellent place to apply recursive programming
                a(n) = sum(map{c}{n},'all');                        % scale transparency by strength
            end
            [~,~,a] = find(a);
            a = a./max(a,[],'all','omitnan');
            
            for n = 1:numel(map{c})
                % Plot significant connections
                [y, x] = find(map{c}{n});
                col = cind.conn + (n-1)/(max(numel(map{c})-1,1)).*[0 0 1; 0 1 0];
                % col = [a(n) 1-a(n) 0; 0 1-a(n), a(n)];
                s(end+1,1) = scatter(ax(2,1), x, y, markerWidth(1)^2, col(c,:), 'filled', 's');
                lgnd(end+1,1) = c;
                
%                 % Highlight mean EC matrix
%                 [sconns(:,1), sconns(:,2)] = find(nbs{thresh(t),c}{n});	% extract significant connections
%                 m(c) = scatter(ax(3,1), sconns(:,2), sconns(:,1), markerWidth(2)^2, cind.node(c,:), 's');
%                 clear sconns
            end
        else
            map{c} = map{c}.*mean(d, 3, 'omitnan')./10;
            
            % Plot significant connections
            [y, x] = find(map{c});
            s(end+1,1) = scatter(ax(2,1), x, y, markerWidth(1)^2, cind.node(c,:), 'filled', 's');
            lgnd(end+1,1) = c;

%             % Highlight mean EC matrix
%             [sconns(:,1), sconns(:,2)] = find(nbs{thresh(t),c}{n});	% extract significant connections
%             m(c) = scatter(ax(3,1), sconns(:,2), sconns(:,1), markerWidth(2)^2, cind.node(c,:), 's');
%             clear sconns
        end
    end
    delete(s(1,:)); s(1,:) = []; lgnd(1,:) = [];
    
    % Render overall network in SPM
    ax(1,1) = subplot(1, 5, [1 2]); hold on;
    plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,i),2)), coords_AAL90, origin, sphereScale, [], map, cind, strcont, strength.summary, rdux);
    sgtitle(F(t,1), ['Threshold: t-statistic = ', num2str(tstat(ind(t,1)))]);
    title(ax(1,1), "All Components", 'FontSize',12);
    
    % Plot legend in mean EC matrix
    legend(ax(2,1), s(:,1), strcont(lgnd), 'Location','bestoutside', 'Color','w', 'FontSize',10);   %  'Orientation','horizontal',
    % legend(ax(3,1), m, strcont, 'Location','bestoutside', 'Color','w', 'FontSize',10);
    clear s lgnd


    %% Set up figure for individual components
    F(t,2) = figure('Position', [0 0 1280 1024]);
    T = tiledlayout(F(t,2), 'flow');
    
    % Render individual components in SPM
    colind = cind;
    for c = 1:size(map,2)
        colind.conn = cind.conn(c,:);
        for f = 1:numel(map{:,c})
            ax(f,2) = nexttile(T); hold on
            plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,i),2)), coords_AAL90, origin, sphereScale, [], map{c}{f}, colind, [], strength.summary, rdux);
            title(strcat("Contrast: ", strcont{c}));
            subtitle(strcat("Component ", num2str(f)));
        end
    end
end
clear c f colind m n map a col climits ncomp T s


%% Save figures
for t = 1:max(unique((ind)))
	saveas(F(t,1), fullfile(path{5}, dirName, strcat(fN{1},"_NBS_Threshold", string(join(strsplit(num2str(tstat(t)),'.'),'')), "_all")), 'png');
    saveas(F(t,2), fullfile(path{5}, dirName, strcat(fN{1},"_NBS_Threshold", string(join(strsplit(num2str(tstat(t)),'.'),'')), "_individual")), 'png');
end
% Save as MATLAB figure
savefig(F, fullfile(path{5}, dirName, strjoin({fN{1},'NBS'},'_')), 'compact');
