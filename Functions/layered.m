function [F] = layered(cortex, nbs, memberships, ttype, h, storarray, N, condName, contrast, cind, fDim, coords_ROI, labels_ROI, origin, sphereScale, strcont, rdux)

%% Visualization Settings

% Set main figure dimensions & indices
i = (fDim(2))*(0:fDim(1)-2);
fInds{1} = sort([i+1 i+2]);
fInds{2} = sort([i+3 i+4 i+5]);
fInds{3} = fInds{2}+3;
fInds{4} = fDim(1)*(fDim(2)-1:fDim(2));
fInds{4}(1) = fInds{4}(1)+1;


%% Visualize each contrast separately

% Locate thresholds & contrasts with significant differences
[~, ind] = find(storarray);
ind = unique(ind); ind = ind(logical(mod(ind,2)));

nm = cell(1,max(unique((ind))));

% Identify contrast of interest
for c = 1:length(ind)
    
    % Generate node color map
    nCol = zeros(N.ROI, 3);
    cont = strjoin(condName(logical(contrast(ind(c),:))), " v. ");          % Generate contrast label
    nCol(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))<0,:) = repmat(cind.node(1,:), sum(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))<0),1);
    nCol(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))>0,:) = repmat(cind.node(2,:), sum(zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2))>0),1);
    
    % Isolate contrasts of interest
   	a = (sum(abs(contrast) == abs(contrast(ind(c),:)),2)) == N.conditions;
    
    % Identify threshold of interest
    thresh = find(storarray(:,ind(c)));
    for t = 1:length(thresh)
        
        %% Combined image
        % Open figure & set axes
        F(ind(c),t,1) = figure('Position', [0 0 1240 1024]);
        ax(1,1) = subplot(fDim(1), fDim(2), fInds{1}); hold on;				% brain space network
        ax(2,1) = subplot(fDim(1), fDim(2), fInds{2}); hold on; grid on;	% adjancecy matrix
        ax(3,1) = subplot(fDim(1), fDim(2), fInds{3}); hold on;				% graph format
        ax(4,1) = subplot(fDim(1), fDim(2), fInds{4}); hold on;				% membership bar chart

        % Set up adjacency matrix chart
        set(ax(2,1), {'YLim','YTick','YTickLabel','FontSize'}, {[0.6 N.ROI+0.5], 1:N.ROI, labels_ROI, 5});
        set(ax(2,1), {'XLim','XTick','XTickLabel','XTickLabelRotation'}, {[0.5 N.ROI+0.5], 5:5:N.ROI, labels_ROI(5:5:N.ROI), -45});
        set(ax(2,1), {'Color', 'FontSize'},{'w', 6});
        title(ax(2,1), "Significant Connections", 'FontSize',12);
        pbaspect(ax(2,1), [1 1 1]);
        
        % Set up node membership chart
        set(ax(3,1),{'Color', 'FontSize'},{'w', 6});
        set(ax(3,1), {'XLim','XTick','XTickLabel','XTickLabelRotation'}, {[1 N.ROI], 1:N.ROI, labels_ROI(5:5:N.ROI), -45});
        yticklabels(ax(3,1), []);
        title(ax(3,1), "Node Memberships", 'FontSize',12);

        % Calculate scatter marker width in points
        currentunits = get(ax(2,1),'Units');
        set(ax(2,1), {'Units'}, {'Points'});
        axpos = get(ax(2,1),'Position');
        set(ax(2,1), 'Units', currentunits); hold on;
        markerWidth = 1/diff(xlim(ax(2,1)))*axpos(3)-1;
        
        % Isolate contrasts of interest
        map = nbs(thresh(t),a);
        s = scatter([],[]);
        nm{c,t} = zeros(N.ROI, N.conditions);

        % Extract significant connections
        for m = 1:numel(map)
            for n = 1:numel(map{m})
                [y, x] = find(map{m}{n});

                % Plot significant connections
                s(1,m) = scatter(ax(2,1), x, y, markerWidth(1)^2, cind.conn(m,:), 'filled', 's');

                % Label which nodes in which communities
                nm{c,t}(unique(y),m) = nm{t}(unique(y),m)+ones(length(unique(y)), 1);
            end
        end
        
        
        % Generate directed graphs for all components
        G = cell(1, size(map,2));
        for c2 = 1:size(map,2)
            G{c2} = sparse(zeros(N.ROI));
            for f = 1:numel(map{:,c2})
                G{c2} = G{c2} + map{c2}{f};
            end
            G{c2} = digraph(G{c2}, labels_ROI);
            plot(ax(3,1), G{c2}, 'Layout','force', 'UseGravity',true, 'EdgeColor',cind.conn(c2,:), 'NodeColor',nCol); hold on
        end
        legend(strcat("Contrast: ", strcont(a)));
        
        
%         % Plot node roles
%         b = bar(ax(4,1), 1:N.ROI, nm{t}');              % connection type
%         s2 = scatter(ax(3,1), 1:N.ROI, zeros(1,N.ROI), 'filled');	% node type
%         for c2 = 1:nnz(a)
%             b(c2).FaceColor = cind.conn(c2,:);
%         end
%         s2.CData = nCol; clear c2 s2
        
        % Render overall network in SPM
        axes(ax(1,1));
        plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2)), coords_ROI, origin, sphereScale, [], map, cind, strcont, [], rdux);
        sgtitle(F(ind(c),t,1), ['Threshold: t-statistic = ', num2str(tstat(thresh(t)))]);
        title(ax(1,1), "All Components", 'FontSize',12);
        
        % Plot legend in connectivity matrix
        legend(ax(2,1), s, strcont(a), 'Location','southoutside', 'Color','w', 'FontSize',10);   %  'Orientation','horizontal',
        clear s lgnd s2
        
        
        %% Individual components in SPM format
        F(ind(c),t,2) = figure('Position', [0 0 1280 1024]);
        T = tiledlayout(F(ind(c),t,2), 'flow');

        % Render individual components
        for c2 = 1:size(map,2)
            for f = 1:numel(map{:,c2})
                % Plot directed graphs per component
                ax(f,2) = nexttile(T); hold on
                
                % Plot directed graphs per component
                plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,cell2mat(h{ttype,cont})),2)), coords_ROI, origin, sphereScale, [], map{c2}{f}, cind, strcont, [], rdux);
                title(strcat("Contrast: ", strcont{c2}), 'FontSize',12);
                subtitle(strcat("Component ", num2str(f)), 'FontSize',8);
            end
        end
        
        
        %% Individual components in graph format
        F(ind(c),t,3) = figure('Position', [0 0 1280 1024]);
        T = tiledlayout(F(ind(c),t,3), 'flow');

        % Render individual components
        for c2 = 1:size(map,2)
            for f = 1:numel(map{:,c2})
                % Generate directed graphs for each component
                G = digraph(map{c2}{f}, labels_ROI);

                % Isolate nodes & links of interest
                [r,cl] = find(map{c2}{f});
                S = subgraph(G, labels_ROI(union(r,cl)));

                % Plot directed graphs per component
                ax(f,3) = nexttile(T); hold on
                plot(S, 'Layout','force', 'UseGravity',true, 'EdgeColor',cind.conn(c2,:), 'NodeColor',nCol(union(r,cl),:));
                title(strcat("Contrast: ", strcont{c2}), 'FontSize',12);
                subtitle(strcat("Component ", num2str(f)), 'FontSize',8);
            end
        end
    end
end
clear c f colind m n map a col climits ncomp T s r cl G S c2 t a fDim fInds


end