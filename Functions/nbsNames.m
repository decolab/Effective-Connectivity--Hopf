function [name] = nbsNames(nbs, ROI)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here


name = nbs;
for c = 1:size(nbs,2)
	for t = 1:size(nbs,1)
        for l = 1:numel(nbs{t,c})
            [row,col] = find(nbs{t,c}{l});
            name{t,c}{l} = horzcat(ROI(col), ROI(row));
            name{t,c}{l} = array2table(name{t,c}{l}, 'VariableNames',{'Origin','Terminus'});
        end
	end
end

end

