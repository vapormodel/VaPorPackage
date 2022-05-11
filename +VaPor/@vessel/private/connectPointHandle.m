function [bins] = connectPointHandle(obj, point, minVertex, otherVertex, type, bins)
%CONNECTPOINT Connects the new point to the tree
%   Detailed explanation goes here

persistent numberVessels

if isempty(numberVessels)
    numberVessels = max(obj.tree(:,1));
end

if type==0
    obj.tree(numberVessels+1,1) = numberVessels+1;
    
    if obj.tree(minVertex,2) > 9
        obj.tree(numberVessels+1,2) = obj.tree(minVertex,2);
    else
        obj.tree(numberVessels+1,2) = obj.tree(minVertex,2)+10;
    end
    obj.tree(numberVessels+1,3:5) = point;
    obj.tree(numberVessels+1,7) = minVertex;
    % Adds new point to vessel tree and adds connection to MinVertex.
    
    numberVessels = numberVessels+1;
    % Increasing the size of VesselTree.
    
    Intersections = cubeIntersect(obj.tree(numberVessels,3:5),obj.tree(obj.tree(numberVessels,7),3:5));
    % Finds intersections of new segment with voxels for bucket list.
    
    for Loop = 1:size(Intersections,1)
        I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
        J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
        K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
        bins{I,J,K}(end+1) = numberVessels; % Adds point to bucket.
        bins{I,J,K} = unique(bins{I,J,K}); % Deletes duplicates in bucket.
    end
    % Update buckets. Adds this new connection to the corresponding
    % buckets.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Connecting New point to OtherVertex %%%%%%
elseif type==1
    
    % Adding new point to new vessel tree and adding connection to
    % OtherVertex.
    obj.tree(numberVessels+1,1) = numberVessels+1;
    if obj.tree(otherVertex,2) > 9
        obj.tree(numberVessels+1,2) = obj.tree(otherVertex,2);
    else
        obj.tree(numberVessels+1,2) = obj.tree(otherVertex,2)+10;
    end
    obj.tree(numberVessels+1,3:5) = point;
    obj.tree(numberVessels+1,7) = otherVertex;
    numberVessels = numberVessels+1; % Increasing the size of VesselTree.
    
    Intersections = cubeIntersect(obj.tree(numberVessels,3:5),obj.tree(obj.tree(numberVessels,7),3:5));
    for Loop = 1:size(Intersections,1)
        I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
        J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
        K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
        bins{I,J,K}(end+1) = numberVessels;
        bins{I,J,K} = unique(bins{I,J,K});
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
    
    obj.tree(numberVessels+1,:) = [numberVessels+1,0,obj.tree(minVertex,3:5)+type*(obj.tree(otherVertex,3:5)-obj.tree(minVertex,3:5)),0,0];
    % Creating a new point part way along the line segment where the
    % intersection lies on the line segment.
    
    if obj.tree(InsertPoint,7) == OtherPoint % Working out which way around the segment should be split in the vessel tree.
        
        % Delete the data from the bucket
        Intersections = cubeIntersect(obj.tree(InsertPoint,3:5),obj.tree(obj.tree(InsertPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
            J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
            K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
            binBool = true(size(bins{I,J,K}));
            binBool(bins{I,J,K}==InsertPoint) = false;
            bins{I,J,K} = bins{I,J,K}(binBool);
            %bins{I,J,K} = sort(bins{I,J,K});
            %bins{I,J,K} = bins{I,J,K}([true;diff(bins{I,J,K}(:))>0]);
        end
        % There will be two new line segments to add to the bucket so
        % the original line that is split has to be removed.
        
        % Add one part back to the bucket.
        obj.tree(InsertPoint,7) = numberVessels+1;
        % update bucket
        Intersections = cubeIntersect(obj.tree(InsertPoint,3:5),obj.tree(obj.tree(InsertPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
            J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
            K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
            bins{I,J,K}(end+1) = InsertPoint;
            bins{I,J,K} = unique(bins{I,J,K});
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
        % Add other part back to the bucket.
        obj.tree(numberVessels+1,7) = OtherPoint;
        if  obj.tree(OtherPoint,2) > 9
            obj.tree(numberVessels+1,2) = obj.tree(OtherPoint,2);
        else
            obj.tree(numberVessels+1,2) = obj.tree(OtherPoint,2)+10;
        end
        Intersections = cubeIntersect(obj.tree(numberVessels+1,3:5),obj.tree(obj.tree(numberVessels+1,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
            J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
            K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
            bins{I,J,K}(end+1) = numberVessels+1;
            bins{I,J,K} = unique(bins{I,J,K});
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
    elseif obj.tree(OtherPoint,7) == InsertPoint % Exactly the same as above but referring to the points the other way around
        
        % Delete the data from the bucket.
        Intersections = cubeIntersect(obj.tree(OtherPoint,3:5),obj.tree(obj.tree(OtherPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
            J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
            K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
            binBool = true(size(bins{I,J,K}));
            binBool(bins{I,J,K}==OtherPoint) = false;
            bins{I,J,K} = bins{I,J,K}(binBool);
            %bins{I,J,K} = unique(bins{I,J,K});
        end
        % There will be two new line segments to add to the bucket so
        % the original line that is split has to be removed.
        
        % Add one part back to the bucket.
        obj.tree(OtherPoint,7) = numberVessels+1;
        Intersections = cubeIntersect(obj.tree(OtherPoint,3:5),obj.tree(obj.tree(OtherPoint,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
            J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
            K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
            bins{I,J,K}(end+1) = OtherPoint;
            %bins{I,J,K} = unique(bins{I,J,K});
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
        % Add other part back to the bucket.
        obj.tree(numberVessels+1,7) = InsertPoint;
        if obj.tree(InsertPoint,2) > 9
            obj.tree(numberVessels+1,2) = obj.tree(InsertPoint,2);
        else
            obj.tree(numberVessels+1,2) = obj.tree(InsertPoint,2)+10;
        end
        Intersections = cubeIntersect(obj.tree(numberVessels+1,3:5),obj.tree(obj.tree(numberVessels+1,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
            J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
            K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
            bins{I,J,K}(end+1) = numberVessels+1;
            %bins{I,J,K} = unique(bins{I,J,K});
        end
        % Update buckets. Adding part of the split line segment back
        % into the corresponding bucket.
        
    end
    % Creating new line segment from the originally generated point.
    
    numberVessels = numberVessels+1; % Increase vessel tree size by 1.
    
    obj.tree(numberVessels+1,1) = numberVessels+1;
    obj.tree(numberVessels+1,3:5) = point;
    obj.tree(numberVessels+1,7) = numberVessels;
    numberVessels = numberVessels+1; % Increase vessel tree size by 1 (two points overall have been created).
    
    Intersections = cubeIntersect(obj.tree(numberVessels,3:5),obj.tree(obj.tree(numberVessels,7),3:5));
    for Loop = 1:size(Intersections,1)
        I = max([1,min([size(obj.domainLimits,1),Intersections(Loop,1)])]);
        J = max([1,min([size(obj.domainLimits,2),Intersections(Loop,2)])]);
        K = max([1,min([size(obj.domainLimits,3),Intersections(Loop,3)])]);
        bins{I,J,K}(end+1) = numberVessels;
        %bins{I,J,K} = unique(bins{I,J,K});
    end
    % Update buckets. Adds this new connection to the corresponding
    % buckets.
    
end
end

