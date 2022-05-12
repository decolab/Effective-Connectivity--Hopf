function [F] = singlefig(cortex, nbs, memberships, h, ind, N, contrast, cind, fPan, fDim, fInds, ROI, origin, strcont)


%% Visualize each contrast separately

% Locate thresholds & contrasts with significant differences
thresh = ind(:,1);
contr = ind(:,2);
ind = unique(contr);
nm = cell(1,max(ind));

% Identify contrast of interest
for c = 1:length(ind)
    
    % Identify significant components
    cont = strsplit(strcont{ind(c)});
    cont = strjoin({cont{1}, cont{3}}, " v. ");		% Generate contrast label
    V = memberships(:,logical(h{:,cont}));
    
    % Generate node color map
    nCol = zeros(N.ROI, 3);
    nCol(zscore(mean(V,2))<0,:) = repmat(cind.node(1,:), sum(zscore(mean(V,2))<0),1);
    nCol(zscore(mean(V,2))>0,:) = repmat(cind.node(2,:), sum(zscore(mean(V,2))>0),1);
    
    % Isolate contrasts of interest
   	a = (sum(abs(contrast) == abs(contrast(ind(c),:)),2)) == N.conditions;
    scont = string(strcont(a));
    
    % Identify threshold of interest
    th = thresh(contr == ind(c));
    for t = 1:length(th)
        
        %% Combined image
        % Open figure & set axes
        F(c,t,1) = figure('Position', fDim(1,:));
        sgtitle("Significant Connections");

        % Set up adjacency matrix chart
        ax(1,1) = subplot(fPan(1), fPan(2), fInds{1}); hold on; grid on;	% adjancecy matrix
        set(ax(1,1), {'YLim','YTick','YTickLabel'}, {[0.6 N.ROI+0.5], 5:5:N.ROI, ROI{5:5:N.ROI,"Label"}});
        set(ax(1,1), {'XLim','XTick','XTickLabel','XTickLabelRotation'}, {[0.5 N.ROI+0.5], 5:5:N.ROI, ROI{5:5:N.ROI,"Label"}, -45});
        set(ax(1,1), {'Color', 'FontSize'},{'w', 6});
        pbaspect(ax(1,1), [1 1 1]);

        % Calculate scatter marker width in points
        currentunits = get(ax(1,1),'Units');
        set(ax(1,1), {'Units'}, {'Points'});
        axpos = get(ax(1,1),'Position');
        set(ax(1,1), 'Units', currentunits); hold on;
        markerWidth = 1/diff(xlim(ax(1,1)))*axpos(3)-1;
        
        % Isolate contrasts of interest
        map = nbs(th(t),a);
        s = scatter([],[]);
        nm{c,t} = zeros(N.ROI, N.conditions);

        % Extract significant connections
        for m = 1:numel(map)
            for n = 1:numel(map{m})
                [y, x] = find(map{m}{n});

                % Plot significant connections
                s(1,m) = scatter(ax(1,1), x, y, markerWidth(1)^2, cind.conn(m,:), 'filled', 's');

                % Label which nodes in which communities
                nm{c,t}(unique(y),m) = nm{t}(unique(y),m)+ones(length(unique(y)), 1);
            end
        end
        
        % Render overall network in SPM
        ax(2,1) = subplot(fPan(1), fPan(2), fInds{2}); hold on;
        plot_nodes_in_cortex(cortex, zscore(mean(V,2)), ROI{:,{'x','y','z'}}, origin, [], map, cind, [], []);
        
        % Plot legend in connectivity matrix
        legend(ax(1,1), s(~cellfun(@isempty, map)), scont(~cellfun(@isempty, map)), 'Location','southoutside', 'Color','w', 'FontSize',10);   %  'Orientation','horizontal',
        clear s lgnd s2
        
        
        %% Individual components in SPM format
        F(c,t,2) = figure('Position', fDim(2,:));
        T = tiledlayout(F(c,t,2), 'flow', 'TileSpacing','compact');

        % Render individual components
        for c2 = 1:size(map,2)
            for f = 1:numel(map{:,c2})
                % update color mapping
                ci = cind;
                ci.conn = cind.conn(c2,:);
                
                % Plot directed graphs per component
                ax(f,2) = nexttile(T); hold on
                plot_nodes_in_cortex(cortex, zscore(mean(V,2)), ROI{:,{'x','y','z'}}, origin, [], map{c2}{f}, ci, [], []);
            end
        end
        clear c2 f ci
    end
end

% Output an empty figure if no valid targets (prevents an error)
if isempty(ind)
    F = [];
end

clear c f colind m n map a col climits ncomp T s r cl G S c2 t a fPan fInds


end