function [inlets] = auto_inlet(vessel_tree)
%UNTITLED Summary of this function goes here
%   Detailed explanation goes here
term_list = false(length(vessel_tree),1);
for i = 1:length(vessel_tree)
    row = vessel_tree(i,:);
    if row(2) == 0
        child_count = sum(row(1)==vessel_tree(:,7));
        if child_count == 0
            term_list(row(1)) = true;
        end
    end
end
term_indx = find(term_list);
% Find Inlets: 
% The inlets are located at the bottom of the tree. 
[~, low_ind] = mink(vessel_tree(term_list,5),3);

%figure
%drawVessels(vessel_tree)
diam = zeros(3,6);

for i = 1:3
    diam(i,:) = vessel_tree(term_indx(low_ind(i)),1:6);
    %diam(i,4) = vessel_tree(term_indx(low_ind(i)),6);
    %scatter3(point_vec(2),point_vec(1),point_vec(3));
end

diam_sort = sortrows(diam,6);

inlets = diam_sort(:,1)';

end

