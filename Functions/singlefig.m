function [F] = singlefig(cortex, nbs, memberships, tstat, ttype, h, ind, N, contrast, cind, fDim, fInds, coords_ROI, labels_ROI, origin, sphereScale, strcont, rdux)


%% Visualize each contrast separately

% Locate thresholds & contrasts with significant differences
thresh = ind(:,1);
contr = ind(:,2);
ind = unique(contr);
nm = cell(1,max(ind));

% Identify contrast of interest
for c = 1:length(ind)
    
    % Generate node color map
    nCol = zeros(N.ROI, 3);
    cont = strsplit(strcont{ind(c)});
    cont = strjoin({cont{1}, cont{3}}, " v. ");         % Generate contrast label
    nCol(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))<0,:) = repmat(cind.node(1,:), sum(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))<0),1);
    nCol(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))>0,:) = repmat(cind.node(2,:), sum(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))>0),1);
    
    % Isolate contrasts of interest
   	a = (sum(abs(contrast) == abs(contrast(ind(c),:)),2)) == N.conditions;
    scont = strcont(a);
    
    % Identify threshold of interest
    th = thresh(contr == ind(c));
    for t = 1:length(th)
        
        %% Combined image
        % Open figure & set axes
        F(c,t,1) = figure('Position', [0 0 1240 1024]);
        ax(3,1) = subplot(fDim(1), fDim(2), fInds{3}); hold on; grid on;	% adjancecy matrix
        ax(1,1) = subplot(fDim(1), fDim(2), fInds{1}); axis off             % graph format
        sgtitle(F(c,t,1), ['Threshold: t-statistic = ', num2str(tstat(th(t)))]);

        % Set up adjacency matrix chart
        set(ax(3,1), {'YLim','YTick','YTickLabel','FontSize'}, {[0.6 N.ROI+0.5], 1:N.ROI, labels_ROI, 5});
        set(ax(3,1), {'XLim','XTick','XTickLabel','XTickLabelRotation'}, {[0.5 N.ROI+0.5], 5:5:N.ROI, labels_ROI(5:5:N.ROI), -45});
        set(ax(3,1), {'Color', 'FontSize'},{'w', 6});
        title(ax(3,1), "Significant Connections", 'FontSize',10);
        pbaspect(ax(3,1), [1 1 1]);

        % Calculate scatter marker width in points
        currentunits = get(ax(3,1),'Units');
        set(ax(3,1), {'Units'}, {'Points'});
        axpos = get(ax(3,1),'Position');
        set(ax(3,1), 'Units', currentunits); hold on;
        markerWidth = 1/diff(xlim(ax(3,1)))*axpos(3)-1;
        
        % Isolate contrasts of interest
        map = nbs(th(t),a);
        s = scatter([],[]);
        nm{c,t} = zeros(N.ROI, N.conditions);

        % Extract significant connections
        for m = 1:numel(map)
            for n = 1:numel(map{m})
                [y, x] = find(map{m}{n});

                % Plot significant connections
                s(1,m) = scatter(ax(3,1), x, y, markerWidth(1)^2, cind.conn(m,:), 'filled', 's');

                % Label which nodes in which communities
                nm{c,t}(unique(y),m) = nm{t}(unique(y),m)+ones(length(unique(y)), 1);
            end
        end
        
        % Generate directed graphs for all components
        g = cell(1, size(map,2));
        G = plot([],[]);
        for c2 = 1:size(map,2)
            g{c2} = sparse(zeros(N.ROI));
            for f = 1:numel(map{:,c2})
                g{c2} = g{c2} + map{c2}{f};
            end
            g{c2} = digraph(g{c2}, labels_ROI);
            G(c2) = plot(ax(1,1), g{c2}, 'Layout','force', 'UseGravity',true, 'EdgeColor',cind.conn(c2,:), 'NodeColor',nCol); hold on
        end
        title(ax(1,1), cont);
        legend(ax(1,1), G, strcat("Contrast: ", scont));
        
        % Render overall network in SPM
        ax(2,1) = subplot(fDim(1), fDim(2), fInds{2}); hold on;
        plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2)), coords_ROI, origin, sphereScale, [], map, cind, [], [], rdux);
        title(ax(2,1), "All Components", 'FontSize',10);
        
        % Plot legend in connectivity matrix
        % legend(ax(3,1), s, scont, 'Location','southoutside', 'Color','w', 'FontSize',10);   %  'Orientation','horizontal',
        clear s lgnd s2
        
        
        %% Individual components in SPM format
        F(c,t,2) = figure('Position', [0 0 1280 1024]);
        T = tiledlayout(F(c,t,2), 'flow');

        % Render individual components
        for c2 = 1:size(map,2)
            for f = 1:numel(map{:,c2})
                % update color mapping
                ci = cind;
                ci.conn = cind.conn(c2,:);
                
                % Plot directed graphs per component
                ax(f,2) = nexttile(T); hold on
                
                % Plot directed graphs per component
                plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2)), coords_ROI, origin, sphereScale, [], map{c2}{f}, ci, [], [], rdux);
                title(strcat("Contrast: ", scont{c2}), 'FontSize',12);
                subtitle(strcat("Component ", num2str(f)), 'FontSize',8);
            end
        end
        clear c2 f ci
        
        
        %% Individual components in graph format
        F(c,t,3) = figure('Position', [0 0 1280 1024]);
        T = tiledlayout(F(c,t,3), 'flow');

        % Render individual components
        for c2 = 1:size(map,2)
            for f = 1:numel(map{:,c2})
                % Generate directed graphs for each component
                g = digraph(map{c2}{f}, labels_ROI);

                % Isolate nodes & links of interest
                [r,cl] = find(map{c2}{f});
                S = subgraph(g, labels_ROI(union(r,cl)));

                % Plot directed graphs per component
                ax(f,3) = nexttile(T); axis off; hold on
                plot(S, 'Layout','force', 'UseGravity',true, 'EdgeColor',cind.conn(c2,:), 'NodeColor',nCol(union(r,cl),:));
                title(strcat("Contrast: ", scont{c2}), 'FontSize',12);
                subtitle(strcat("Component ", num2str(f)), 'FontSize',8);
            end
        end
    end
end
clear c f colind m n map a col climits ncomp T s r cl G S c2 t a fDim fInds


end