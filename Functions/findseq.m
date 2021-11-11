function [pair] = findseq(storarray, np)
%   FINDPAIR returns paired contrast indices in the variable pair.  If no
% contrast returns paired results, pair is empty.

% Find paired results
pair = cell(size(storarray,2), 1);
for k = 1:np:size(storarray,2)
    [r,c] = find(storarray(:,k:k+(np-1)));
    if numel(unique(c)) == np && numel(unique(c)) > 1
        co = [k:k+(np-1)]';
        c = co(c);
        [~, ia, ~] = unique(r);
        id = setdiff(1:size(r,1), ia, 'stable');
        id = ismember(r,r(id), 'rows');
        pair{k} = horzcat(r(id), c(id));
    elseif numel(unique(c)) == np
        co = [k:k+(np-1)]';
        pair{k} = horzcat(r, co(c));
    end
end
clear r c ia id k co

% Remove non-paired results
pair = pair(~cellfun(@isempty, pair));

% Convert to array
pair = cell2mat(pair);

end

