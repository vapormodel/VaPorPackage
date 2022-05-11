function [vessel_tree] = autofit(vessel_tree,target_volume)
vessel_tree = vesselAltRotate(vessel_tree,'x',90);

% Create the Mask:
y_mask = squeeze(any(target_volume,1));
y_lower = find(squeeze(any(y_mask,2)),1);
y_upper = find(squeeze(any(y_mask,2)), 1, 'last');
z_mask = squeeze(any(target_volume,2));
z_lower = find(squeeze(any(z_mask,1)),1);
z_upper = find(squeeze(any(z_mask,1)), 1, 'last');
x_mask = squeeze(any(target_volume,3));
x_lower = find(squeeze(any(x_mask,2)),1);
x_upper = find(squeeze(any(x_mask,2)), 1, 'last');
%planecut('x',1,z_mask)

% Find Min X in tree:
tree_x_min = 20e3;
for i = 1:length(vessel_tree)
    this_row = vessel_tree(i,:);
    %if this_row(2) ~= 0
        if this_row(3) < tree_x_min
            tree_x_min = this_row(3);
        end
    %end
end
tree_x_max = -20e3;
for i = 1:length(vessel_tree)
    this_row = vessel_tree(i,:);
    %if this_row(2) ~= 0
        if this_row(3) > tree_x_max
            tree_x_max = this_row(3);
        end
    %end
end
vessel_tree = vesselAltScale(vessel_tree,(x_upper-x_lower)/(tree_x_max-tree_x_min),'x');

tree_y_min = 20e3;
for i = 1:length(vessel_tree)
    this_row = vessel_tree(i,:);
    %if this_row(2) ~= 0
        if this_row(4) < tree_y_min
            tree_y_min = this_row(4);
        end
    %end
end
tree_y_max = -20e3;
for i = 1:length(vessel_tree)
    this_row = vessel_tree(i,:);
    %if this_row(2) ~= 0
        if this_row(4) > tree_y_max
            tree_y_max = this_row(4);
        end
    %end
end
vessel_tree = vesselAltScale(vessel_tree,(y_upper-y_lower)/(tree_y_max-tree_y_min),'y');

tree_z_min = 20e3;
for i = 1:length(vessel_tree)
    this_row = vessel_tree(i,:);
    %if this_row(2) ~= 0
        if this_row(5) < tree_z_min
            tree_z_min = this_row(5);
        end
    %end
end
tree_z_max = -20e3;
for i = 1:length(vessel_tree)
    this_row = vessel_tree(i,:);
    %if this_row(2) ~= 0
        if this_row(5) > tree_z_max
            tree_z_max = this_row(5);
        end
    %end
end
vessel_tree = vesselAltScale(vessel_tree,(z_upper-z_lower)/(tree_z_max-tree_z_min),'z');

tree_x_centre = (tree_x_min+tree_x_max)/2;
tree_y_centre = (tree_y_min+tree_y_max)/2;
tree_z_centre = (tree_z_min+tree_z_max)/2;
dom_x_centre = (x_lower+x_upper)/2;
dom_y_centre = (y_lower+y_upper)/2;
dom_z_centre = (z_lower+z_upper)/2;

vessel_tree = vesselAltTranslate(vessel_tree,[dom_x_centre-tree_x_centre,dom_y_centre-tree_y_centre,dom_z_centre-tree_z_centre]);

end