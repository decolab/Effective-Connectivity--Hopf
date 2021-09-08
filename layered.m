%% Set parameters

% Legends
strcont = {strjoin([labels.Properties.VariableNames(1), '<', labels.Properties.VariableNames(2)]), strjoin([labels.Properties.VariableNames(1), '>', labels.Properties.VariableNames(2)])};

% Set brain rendering parameters
cortex.file = fullfile(path{3},'OCD','Data','MNI152_T1_2mm_brain_mask.nii');	% file containing cortical atlas
cortex.color = [0.9 0.9 0.9];					% color for cortical rendering
cortex.transparency = 0.1;						% set to 1 for opaque cortex
cortex.val = 0.3;								% set isonormal line spacing
cortex.view = [-90 90];							% set camera angle
rdux = 0.7;						% proportion of surface faces to keep (<1)

% set color index
cind.node = [1 0 0; 0 0 1];
cind.conn = ['m'; 'c'];


%% locate significant components
j = h{:,'IC'};
if numel(j) > 1
	i = j{1};
	for k = 2:numel(j)
		i = union(i, j{k});
	end
	i = unique(i);
else
	i = j{:};
end
clear j


%% Extract individual networks

[ind(:,1), ind(:,2)] = find(storarray);
fN = strsplit(fileName, '_');

nm = cell(1,max(unique((ind))));
compROIs = nm;
for t = 1:max(unique((ind)))
    
    % Open figure
    F(t,1) = figure('Position', [0 0 1240 1024]);
    
    % Plot significant connections as binarized connectivity map
    ax(2,1) = subplot(4, 5, [3:5, 8:10, 13:15]); hold on; grid on
	set(ax(2,1), {'YTick','YTickLabel','FontSize'}, {1:N.ROI, label_ROI, 5});
    set(ax(2,1), {'XTick','XTickLabel','XTickLabelRotation'}, {5:5:N.ROI, label_ROI(5:5:N.ROI), -90});
    set(ax(2,1),{'Color', 'FontSize'},{'w', 6});
    xlim([0.5 N.ROI+0.5]); ylim([0.6 N.ROI+0.5]);
    title("Significant Connections", 'FontSize',12);
    pbaspect([1 1 1]);
    
    % Plot node memberships as bar chart
    ax(3,1) = subplot(4, 5, 16:20); hold on
    set(ax(3,1),{'Color', 'FontSize'},{'w', 6});
    title("Node Memberships", 'FontSize',12);
    
    
    % Calculate scatter Marker width in points
    currentunits = get(ax(2,1),'Units');
    set(ax(2,1), {'Units'}, {'Points'});
	axpos = get(ax(2,1),'Position');
	set(ax(2,1), 'Units', currentunits); hold on;
	markerWidth = 1/diff(xlim(ax(2,1)))*axpos(3)-1;
    
    % Compute mean distance map & get number of components per threshold
    map = nbs(ind(t,1),:); s = scatter([],[]); lgnd  = 0;
    nm{t} = cell(size(map));
    for c = 1:numel(map)
%         for n = 1:numel(map{c})
%             map{c}{n} = map{c}{n}.*mean(d, 3, 'omitnan')./10;	% this would be an excellent place to apply recursive programming
%             a(n) = sum(map{c}{n},'all');                        % scale transparency by strength
%         end
%         [~,~,a] = find(a);
%         a = a./max(a,[],'all','omitnan');

        % Extract significant connections
        nm{t}{c} = zeros(numel(label_ROI), numel(map{c}));
        for n = 1:numel(map{c})
            [y, x] = find(map{c}{n});
            col = cind.conn + (n-1)/(max(numel(map{c})-1,1)).*[0 0 1; 0 1 0];
            % col = [a(n) 1-a(n) 0; 0 1-a(n), a(n)];
            
            % Plot significant connections
            s(end+1,1) = scatter(ax(2,1), x, y, markerWidth(1)^2, cind.conn(c,:), 'filled', 's');
            lgnd(end+1,1) = c;
            
            % Label which nodes in which communities
            nm{t}{c}(unique(y),n) = ones(length(unique(y)), 1);
        end
        nm{t}{c} = (-nm{t}{c}).^(c+1);
    end
    delete(s(1,:)); s(1,:) = []; lgnd(1,:) = [];
    
    % Plot node roles
    nm{t} = cell2mat(nm{t});
    bar(ax(3,1), categorical(label_ROI), nm{t});
    
    % List nodes in each component
    compROIs{t} = repmat(label_ROI, [1 size(nm{t},2)]);
    compROIs{t}(~logical(nm{t})) = "-";
    
    % Render overall network in SPM
    ax(1,1) = subplot(4, 5, [1:2, 6:7, 11:12]); hold on;
    plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,i),2)), coords_ROI, origin, labels_ROI, sphereScale, [], map, cind, strcont, [], rdux);
    sgtitle(F(t,1), ['Threshold: t-statistic = ', num2str(tstat(ind(t,1)))]);
    title(ax(1,1), "All Components", 'FontSize',12);
    
    % Plot legend in connectivity matrix
    legend(ax(2,1), s(:,1), strcont(lgnd), 'Location','bestoutside', 'Color','w', 'FontSize',10);   %  'Orientation','horizontal',
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
            plot_nodes_in_cortex(cortex, zscore(mean(memberships(:,i),2)), coords_ROI, origin, labels_ROI, sphereScale, [], map{c}{f}, colind, [], [], rdux);
            title(strcat("Contrast: ", strcont{c}), 'FontSize',12);
            subtitle(strcat("Component ", num2str(f)), 'FontSize',8);
        end
    end
end
clear c f colind m n map a col climits ncomp T s


%% Save figures
% for t = 1:max(unique((ind)))
% 	saveas(F(t,1), fullfile(path{5}, dirName, strcat(fN{1},"_NBS_Threshold", string(join(strsplit(num2str(tstat(t)),'.'),'')), "_all")), 'png');
%     saveas(F(t,2), fullfile(path{5}, dirName, strcat(fN{1},"_NBS_Threshold", string(join(strsplit(num2str(tstat(t)),'.'),'')), "_individual")), 'png');
% end
% % Save as MATLAB figure
% savefig(F, fullfile(path{5}, dirName, strjoin({fN{1},'NBS'},'_')), 'compact');
