function [vesselTree, bins] = connectPoint(vesselTree, domain, point, minVertex, otherVertex, type, bins, attach_fam)
%CONNECTPOINT Connects the new point to the tree
%   Detailed explanation goes here
if nargin ~= 8
    attach_fam = 1;
end
persistent numberVessels

%bins = bins(1+(attach_fam-1)*numel(domain):(attach_fam)*numel(domain),:);

if isempty(numberVessels)
    numberVessels = max(vesselTree(:,1));
end

if type==0
    vesselTree(numberVessels+1,1) = numberVessels+1;
    
    if vesselTree(minVertex,2) > 9
        vesselTree(numberVessels+1,2) = vesselTree(minVertex,2);
    else
        vesselTree(numberVessels+1,2) = vesselTree(minVertex,2)+10;
    end
    vesselTree(numberVessels+1,3:5) = point;
    vesselTree(numberVessels+1,7) = minVertex;
    % Adds new point to vessel tree and adds connection to MinVertex.
    
    numberVessels = numberVessels+1;
    % Increasing the size of VesselTree.
    
    Intersections = cubeIntersect(vesselTree(numberVessels,3:5),vesselTree(vesselTree(numberVessels,7),3:5));
    % Finds intersections of new segment with voxels for bucket list.
    
    for Loop = 1:size(Intersections,1)
        I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
        J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
        K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
        %bins{I,J,K}(end+1) = numberVessels; % Adds point to bucket.
        %bins{I,J,K} = unique(bins{I,J,K}); % Deletes duplicates in bucket.
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) + 1;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = numberVessels;
    end
    % Update buckets. Adds this new connection to the corresponding
    % buckets.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Connecting New point to OtherVertex %%%%%%
elseif type==1
    
    % Adding new point to new vessel tree and adding connection to
    % OtherVertex.
    vesselTree(numberVessels+1,1) = numberVessels+1;
    if vesselTree(otherVertex,2) > 9
        vesselTree(numberVessels+1,2) = vesselTree(otherVertex,2);
    else
        vesselTree(numberVessels+1,2) = vesselTree(otherVertex,2)+10;
    end
    vesselTree(numberVessels+1,3:5) = point;
    vesselTree(numberVessels+1,7) = otherVertex;
    numberVessels = numberVessels+1; % Increasing the size of VesselTree.
    
    Intersections = cubeIntersect(vesselTree(numberVessels,3:5),vesselTree(vesselTree(numberVessels,7),3:5));
    for Loop = 1:size(Intersections,1)
        I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
        J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
        K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) + 1;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = numberVessels;
    end
    % Update buckets. Adds this new connection to the corresponding
    % buckets.
    
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%% Connecting New Point to Line Segment %%%%%
else
    
    InsertPoint = min(minVertex,otherVertex); % Point above insertion (in terms of going down list).
    OtherPoint = max(minVertex,otherVertex); % Point on other end of segment.
    % Redefining the connection points in terms of where they lie in
    % the vessel tree list. Same points just might be referred to the
    % other way around.
    
    vesselTree(numberVessels+1,:) = [numberVessels+1,0,vesselTree(minVertex,3:5)+type*(vesselTree(otherVertex,3:5)-vesselTree(minVertex,3:5)),0,0];
    % Creating a new point part way along the line segment where the
    % intersection lies on the line segment.
    
    if vesselTree(InsertPoint,7) == OtherPoint % Working out which way around the segment should be split in the vessel tree.
        
        % Delete the data from the bucket
        Intersections = cubeIntersect(vesselTree(InsertPoint,3:5),vesselTree(vesselTree(InsertPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
%             binBool = true(size(bins{I,J,K}));
%             binBool(bins{I,J,K}==InsertPoint) = false;
%             bins{I,J,K} = bins{I,J,K}(binBool);
        subset_bins = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),2:bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1));
        idx_remove = find(subset_bins==InsertPoint);
        new_subset_bins = [subset_bins(1:idx_remove-1),subset_bins(idx_remove+1:end), 0];
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),2:bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = new_subset_bins;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) - 1;
            %bins{I,J,K} = sort(bins{I,J,K});
            %bins{I,J,K} = bins{I,J,K}([true;diff(bins{I,J,K}(:))>0]);
        end
        % There will be two new line segments to add to the bucket so
        % the original line that is split has to be removed.
        
        % Add one part back to the bucket.
        vesselTree(InsertPoint,7) = numberVessels+1;
        % update bucket
        Intersections = cubeIntersect(vesselTree(InsertPoint,3:5),vesselTree(vesselTree(InsertPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
            %bins{I,J,K}(end+1) = InsertPoint;
            %bins{I,J,K} = unique(bins{I,J,K});
            bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) + 1;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = InsertPoint;
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
        % Add other part back to the bucket.
        vesselTree(numberVessels+1,7) = OtherPoint;
        if  vesselTree(OtherPoint,2) > 9
            vesselTree(numberVessels+1,2) = vesselTree(OtherPoint,2);
        else
            vesselTree(numberVessels+1,2) = vesselTree(OtherPoint,2)+10;
        end
        Intersections = cubeIntersect(vesselTree(numberVessels+1,3:5),vesselTree(vesselTree(numberVessels+1,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
            %bins{I,J,K}(end+1) = numberVessels+1;
            %bins{I,J,K} = unique(bins{I,J,K});
            bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) + 1;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = numberVessels+1;
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
    elseif vesselTree(OtherPoint,7) == InsertPoint % Exactly the same as above but referring to the points the other way around
        
        % Delete the data from the bucket.
        Intersections = cubeIntersect(vesselTree(OtherPoint,3:5),vesselTree(vesselTree(OtherPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
            %binBool = true(size(bins{I,J,K}));
            %binBool(bins{I,J,K}==OtherPoint) = false;
            %bins{I,J,K} = bins{I,J,K}(binBool);
            %bins{I,J,K} = unique(bins{I,J,K});
            subset_bins = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),2:bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1));
        idx_remove = find(subset_bins==OtherPoint);
        new_subset_bins = [subset_bins(1:idx_remove-1),subset_bins(idx_remove+1:end), 0];
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),2:bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = new_subset_bins;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) - 1;
        end
        % There will be two new line segments to add to the bucket so
        % the original line that is split has to be removed.
        
        % Add one part back to the bucket.
        vesselTree(OtherPoint,7) = numberVessels+1;
        Intersections = cubeIntersect(vesselTree(OtherPoint,3:5),vesselTree(vesselTree(OtherPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
            %bins{I,J,K}(end+1) = OtherPoint;
            %bins{I,J,K} = unique(bins{I,J,K});
            bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) + 1;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = OtherPoint;
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
        % Add other part back to the bucket.
        vesselTree(numberVessels+1,7) = InsertPoint;
        if vesselTree(InsertPoint,2) > 9
            vesselTree(numberVessels+1,2) = vesselTree(InsertPoint,2);
        else
            vesselTree(numberVessels+1,2) = vesselTree(InsertPoint,2)+10;
        end
        Intersections = cubeIntersect(vesselTree(numberVessels+1,3:5),vesselTree(vesselTree(numberVessels+1,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
            %bins{I,J,K}(end+1) = numberVessels+1;
            %bins{I,J,K} = unique(bins{I,J,K});
            bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) + 1;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = numberVessels+1;
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
    end
    % Creating new line segment from the originally generated point.
    
    numberVessels = numberVessels+1; % Increase vessel tree size by 1.
    
    vesselTree(numberVessels+1,1) = numberVessels+1;
    if vesselTree(numberVessels,2) > 9
        vesselTree(numberVessels+1,2) = vesselTree(numberVessels,2);
    else
        vesselTree(numberVessels+1,2) = vesselTree(numberVessels,2)+10;
    end
    vesselTree(numberVessels+1,3:5) = point;
    vesselTree(numberVessels+1,7) = numberVessels;
    numberVessels = numberVessels+1; % Increase vessel tree size by 1 (two points overall have been created).
    
    Intersections = cubeIntersect(vesselTree(numberVessels,3:5),vesselTree(vesselTree(numberVessels,7),3:5));
    for Loop = 1:size(Intersections,1)
        I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
        J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
        K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
        %bins{I,J,K}(end+1) = numberVessels;
        %bins{I,J,K} = unique(bins{I,J,K});
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1) + 1;
        bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(attach_fam-1)*numel(domain),1)) = numberVessels;
    end
    % Update buckets. Adds this new connection to the corresponding
    % buckets.
    
end
%bins(1+(attach_fam-1)*numel(domain):(attach_fam)*numel(domain),:) = bins;
end

