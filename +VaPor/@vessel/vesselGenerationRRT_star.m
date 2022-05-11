function [VesselOut,ProbMapDensityCumSum] = vesselGenerationRRT_star(VesselIn,ProbMap,Domain,NoIterations,WeightFactor,Epsilon)


%%%%%% Preallocate Space for Vessel Tree %%%%%%%%
NoPoints = size(VesselIn,1); % Number of points comprising the current vessel tree.
VesselOut = zeros(NoPoints+2*NoIterations,7); % Preallocate maximum possible size of new vessel tree (current size + 2*no_iter).
VesselOut(1:NoPoints,:) = VesselIn; % Insert original tree into new tree to start.
% Initialising the new vessel tree.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Generate Probability Map %%%%%%%%%%%%%%%%%
ProbMap(~Domain) = 0; % Everything outside the desired domain has zero probablilty.
ProbMap(isnan(ProbMap)) = 0; % Remove any NaNs from the data and set them to zero.
% Initialising the vessel generation probablity map

ProbMapLogical = logical(ProbMap);
ProbMapIndex = zeros(numel(ProbMapLogical(ProbMapLogical)),3);
[ProbMapIndex(:,1),ProbMapIndex(:,2),ProbMapIndex(:,3)] = ind2sub(size(ProbMapLogical),find(ProbMapLogical));
% Indexing points for ease of retrieval.

ProbMapDensity = (ProbMap(ProbMapLogical).^WeightFactor)/(sum(ProbMap(ProbMapLogical).^WeightFactor));
ProbMapDensityCumSum = cumsum(ProbMapDensity);
% Creating cummulitive sum of probablity map.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Create Buckets for Searching %%%%%%%%%%%%%
Bucket = cell(size(Domain));
Bucket{max([1,min([size(Domain,1),floor(VesselOut(1,3))])]),max([1,min([size(Domain,2),floor(VesselOut(1,4))])]),max([1,min([size(Domain,3),floor(VesselOut(1,5))])])} = 1;
for M = 2:NoPoints
    Intersections = cubeIntersect(VesselOut(M,3:5),VesselOut(VesselOut(M,7),3:5));
    for Loop = 1:size(Intersections,1)
        I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
        J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
        K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
        Bucket{I,J,K}(end+1) = M;
        %Bucket{I,J,K} = unique(Bucket{I,J,K});
        if ~isempty(Bucket{I,J,K})
            Bucket{I,J,K} = sort(Bucket{I,J,K});
            Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
        end
    end
end
% Initialising the Buckets for quicker searching. It would be far to
% inefficient to search every segment to determine the closest segment.
% Instead each voxel in the domain acts a 'bucket' which contains a list of
% all segments that intersect the voxel. By only searching within the
% closest buckets first, the search time is dramatically reduced. Once the
% closest buckets are exhausted (contain no segments), the search is
% gradually expanded to include neighbouring buckets.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Create Checks for Progress Update %%%%%%%%
Check25 = false; Check50 = false; Check75 = false;
% Values used for the percent completion checking.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% Pick the voxels to use in advance:
voxelsToGeneratedNodes = randsample(length(ProbMapDensity), NoIterations, true, ProbMapDensity)';

N = 0;
%%%%%% Generate Vessel Tree %%%%%%%%%%%%%%%%%%%%%
%for N = 1:NoIterations % For the number of specified iterations.
for thisValue = voxelsToGeneratedNodes
    
    N = N + 1;
    
    %%%%%% Generate New Point %%%%%%%%%%%%%%%%%%%%%%%
    %New = ProbMapIndex(find(ProbMapDensityCumSum>=rand,1,'first'),:); % Find volume to generate point inside.
    New = ProbMapIndex(thisValue,:);
    New = New + [rand-0.5,rand-0.5,rand-0.5]; % Randomly distibute point around volume.
    % Randomly generating new point based on weight function.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Serch Nearby Buckets for Segments %%%%%%%%
    BucketList = []; % Initialise list of closest segments.
    Range = 1; % Distance from original bucket to search from.
    Check = 0; % Variable for breaking loop.
    while Check == 0
        X_Lim = [max(1,floor(New(1))-Range),min(size(Domain,1),floor(New(1))+Range)];
        Y_Lim = [max(1,floor(New(2))-Range),min(size(Domain,2),floor(New(2))+Range)];
        Z_Lim = [max(1,floor(New(3))-Range),min(size(Domain,3),floor(New(3))+Range)];
        % Setting limits for bucket search.
        %         for I = X_Lim(1):X_Lim(2)
        %             for J = Y_Lim(1):Y_Lim(2)
        %                 for K = Z_Lim(1):Z_Lim(2)
        %                     if ~isempty(Bucket{I,J,K})
        %                         BucketList = [BucketList, Bucket{I,J,K}]; % Accumulate any segments found in a list.
        %                     end
        %                     % Searching through every bucket for possible vessel
        %                     % connections.
        %                 end
        %             end
        %         end
        searchArea = Bucket(X_Lim(1):X_Lim(2), Y_Lim(1):Y_Lim(2), Z_Lim(1):Z_Lim(2));
        searchArea = searchArea(~cellfun('isempty',searchArea));
        listLength = sum(cellfun('length',searchArea),1);
        if listLength ~= 0
            BucketList = zeros(listLength,1);
            bucketPos = 1;
            for counter = 1:length(searchArea)
                numInsert = numel(searchArea{counter});
                BucketList(bucketPos:bucketPos+numInsert-1) = searchArea{counter};
                bucketPos = bucketPos + numInsert;
            end
        else
            BucketList = [];
        end
        
        
        if ~isempty(BucketList)
            Check = 1; % If list is not empty, break the search.
        else
            Range = Range + 1; % Increase search range if list is empty.
        end
        % Checking if any vessels have been found.
        
        if Range > max(size(Domain))
            BucketList = Vessel(:,1)';
            Check = 1;
            % If no connections are found within the domain. Search the
            % whole tree (as it lies fully outside the domain). Very
            % inefficient if a large amount of vessels lie outside of the
            % domain for generation.
        end
        
    end
    %BucketList = unique(BucketList); % Removing duplicate connections found.
    
    BucketList = sort(BucketList);
    BucketList = BucketList([true;diff(BucketList(:))>0]);
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Find Closest Segment in Results %%%%%%%%%%
    Dist = zeros(size(BucketList,1),2);
    for M2 = 1:size(BucketList,1)
        M = BucketList(M2);
        if M == 1
            Dist(M2,:) = [norm(New-VesselOut(1,3:5)),0];
        else
            vesselOutVector = VesselOut(M,:);
            param1 = vesselOutVector(3:5);
            param2 = VesselOut(vesselOutVector(7),3:5);
            [D,T] = distPointToLineSeg(param1,param2,New);
            Dist(M2,:) = [D,T];
        end
    end
    % Find the distances between the point and all segments listed in the
    % BucketList.
    
    %MinVertex = find(Dist(:,1)==min(Dist(:,1)),1);
    [~, MinVertex] = min(Dist(:,1));
    % Find the segment with shortest distance. In case of equal distance
    % picks earliest occurance vertex.
    
    T = Dist(MinVertex,2);
    % How far the intersection is along the line segment. If T = 0 or T = 1
    % then the shortest distance is one of the two endpoints of the
    % segment.
    
    MinVertex = BucketList(MinVertex);
    % The node number of the vessel tree that defines the segment to which
    % the new point will be attached.
    
    OtherVertex = VesselOut(MinVertex,7);
    % The point at the other end of the line segment defined by MinVertex.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    
    %%%%%% Restricting Size of Segment %%%%%%%%%%%%%%
    if T==0
        PointOfConnection = VesselOut(MinVertex,3:5);
    elseif T == 1
        PointOfConnection = VesselOut(VesselOut(MinVertex,7),3:5);
    else
        PointOfConnection = VesselOut(MinVertex,3:5)+T*(VesselOut(VesselOut(MinVertex,7),3:5)-VesselOut(MinVertex,3:5));
    end
    
    if norm(New - PointOfConnection) > Epsilon
        New = PointOfConnection + Epsilon/norm(New - PointOfConnection)*(New - PointOfConnection);
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Connecting New point to MinVertex %%%%%%%%
    if T==0
        
        VesselOut(NoPoints+1,1) = NoPoints+1;
        %VesselOut(NoPoints+1,2) = VesselOut(MinVertex,2)+10; % Assign the region of the new vertex to the region of MinVertex.
        if VesselOut(MinVertex,2) > 9
            VesselOut(NoPoints+1,2) = VesselOut(MinVertex,2);
        else
            VesselOut(NoPoints+1,2) = VesselOut(MinVertex,2)+10;
        end
        VesselOut(NoPoints+1,3:5) = New;
        VesselOut(NoPoints+1,7) = MinVertex;
        % Adds new point to vessel tree and adds connection to MinVertex.
        
        NoPoints = NoPoints+1;
        % Increasing the size of VesselTree.
        
        Intersections = cubeIntersect(VesselOut(NoPoints,3:5),VesselOut(VesselOut(NoPoints,7),3:5));
        % Finds intersections of new segment with voxels for bucket list.
        addList = zeros(size(Intersections,1),3);
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
            addList(Loop, :) = [I, J, K];
            Bucket{I,J,K}(end+1) = NoPoints; % Adds point to bucket.
            %Bucket{I,J,K} = unique(Bucket{I,J,K}); % Deletes duplicates in bucket.
            %Bucket{I,J,K} = sort(Bucket{I,J,K});
            %Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
        end
        %BucketRes = cellfun(@removeDuplicates,Bucket, 'UniformOutput', false);
        % Update buckets. Adds this new connection to the corresponding
        % buckets.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
       % Bucket = BucketRes;
       %Bucket = removeDuplicates(Bucket, addList);
        
        
        %%%%%% Connecting New point to OtherVertex %%%%%%
    elseif T==1
        
        % Adding new point to new vessel tree and adding connection to
        % OtherVertex.
        VesselOut(NoPoints+1,1) = NoPoints+1;
        %VesselOut(NoPoints+1,2) = VesselOut(OtherVertex,2)+10; % Assign the region of the new vertex to the the region of OtherVertex
        if VesselOut(OtherVertex,2) > 9
            VesselOut(NoPoints+1,2) = VesselOut(OtherVertex,2);
        else
            VesselOut(NoPoints+1,2) = VesselOut(OtherVertex,2)+10;
        end
        VesselOut(NoPoints+1,3:5) = New;
        VesselOut(NoPoints+1,7) = OtherVertex;
        NoPoints = NoPoints+1; % Increasing the size of VesselTree.
        
        Intersections = cubeIntersect(VesselOut(NoPoints,3:5),VesselOut(VesselOut(NoPoints,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
            Bucket{I,J,K}(end+1) = NoPoints;
            %Bucket{I,J,K} = unique(Bucket{I,J,K});
            %Bucket{I,J,K} = sort(Bucket{I,J,K});
            %Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
        end
        % Update buckets. Adds this new connection to the corresponding
        % buckets.
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%% Connecting New Point to Line Segment %%%%%
    else
        
        InsertPoint = min(MinVertex,OtherVertex); % Point above insertion (in terms of going down list).
        OtherPoint = max(MinVertex,OtherVertex); % Point on other end of segment.
        % Redefining the connection points in terms of where they lie in
        % the vessel tree list. Same points just might be referred to the
        % other way around.
        
        VesselOut(NoPoints+1,:) = [NoPoints+1,0,VesselOut(MinVertex,3:5)+T*(VesselOut(OtherVertex,3:5)-VesselOut(MinVertex,3:5)),0,0];
        % Creating a new point part way along the line segment where the
        % intersection lies on the line segment.
        
        if VesselOut(InsertPoint,7) == OtherPoint % Working out which way around the segment should be split in the vessel tree.
            
            % Delete the data from the bucket
            Intersections = cubeIntersect(VesselOut(InsertPoint,3:5),VesselOut(VesselOut(InsertPoint,7),3:5));
            for Loop = 1:size(Intersections,1)
                I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
                J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
                K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
                BucketBool = true(size(Bucket{I,J,K}));
                BucketBool(Bucket{I,J,K}==InsertPoint) = false;
                Bucket{I,J,K} = Bucket{I,J,K}(BucketBool);
                %Bucket{I,J,K} = unique(Bucket{I,J,K});
%                 if ~isempty(Bucket{I,J,K})
%                 Bucket{I,J,K} = sort(Bucket{I,J,K});
%                 Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
%                 end
            end
            % There will be two new line segments to add to the bucket so
            % the original line that is split has to be removed.
            
            % Add one part back to the bucket.
            VesselOut(InsertPoint,7) = NoPoints+1;
            % update bucket
            Intersections = cubeIntersect(VesselOut(InsertPoint,3:5),VesselOut(VesselOut(InsertPoint,7),3:5));
            for Loop = 1:size(Intersections,1)
                I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
                J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
                K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
                Bucket{I,J,K}(end+1) = InsertPoint;
                %Bucket{I,J,K} = unique(Bucket{I,J,K});
%                 Bucket{I,J,K} = sort(Bucket{I,J,K});
%                 Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
            end
            % Update buckets. Adding part of the split line segment back
            % into the corresponding bucket.
            
            % Add other part back to the bucket.
            VesselOut(NoPoints+1,7) = OtherPoint;
            %VesselOut(NoPoints+1,2) = VesselOut(OtherPoint,2); % Assign new point to region of OtherPoint.
            if  VesselOut(OtherPoint,2) > 9
                VesselOut(NoPoints+1,2) = VesselOut(OtherPoint,2);
            else
                VesselOut(NoPoints+1,2) = VesselOut(OtherPoint,2)+10;
            end
            Intersections = cubeIntersect(VesselOut(NoPoints+1,3:5),VesselOut(VesselOut(NoPoints+1,7),3:5));
            for Loop = 1:size(Intersections,1)
                I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
                J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
                K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
                Bucket{I,J,K}(end+1) = NoPoints+1;
                %Bucket{I,J,K} = unique(Bucket{I,J,K});
%                 Bucket{I,J,K} = sort(Bucket{I,J,K});
%                 Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
            end
            % Update buckets. Adding part of the split line segment back
            % into the corresponding bucket.
            
        elseif VesselOut(OtherPoint,7) == InsertPoint % Exactly the same as above but referring to the points the other way around
            
            % Delete the data from the bucket.
            Intersections = cubeIntersect(VesselOut(OtherPoint,3:5),VesselOut(VesselOut(OtherPoint,7),3:5));
            for Loop = 1:size(Intersections,1)
                I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
                J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
                K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
                BucketBool = true(size(Bucket{I,J,K}));
                BucketBool(Bucket{I,J,K}==OtherPoint) = false;
                Bucket{I,J,K} = Bucket{I,J,K}(BucketBool);
                %Bucket{I,J,K} = unique(Bucket{I,J,K});
%                 if ~isempty(Bucket{I,J,K})
%                     Bucket{I,J,K} = sort(Bucket{I,J,K});
%                     Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
%                 end
            end
            % There will be two new line segments to add to the bucket so
            % the original line that is split has to be removed.
            
            % Add one part back to the bucket.
            VesselOut(OtherPoint,7) = NoPoints+1;
            Intersections = cubeIntersect(VesselOut(OtherPoint,3:5),VesselOut(VesselOut(OtherPoint,7),3:5));
            for Loop = 1:size(Intersections,1)
                I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
                J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
                K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
                Bucket{I,J,K}(end+1) = OtherPoint;
                %Bucket{I,J,K} = unique(Bucket{I,J,K});
%                 Bucket{I,J,K} = sort(Bucket{I,J,K});
%                 Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
            end
            % Update buckets. Adding part of the split line segment back
            % into the corresponding bucket.
            
            % Add other part back to the bucket.
            VesselOut(NoPoints+1,7) = InsertPoint;
            %VesselOut(NoPoints+1,2) = VesselOut(InsertPoint,2); % Assign new point to the region of InsertPoint.
            if VesselOut(InsertPoint,2) > 9
                VesselOut(NoPoints+1,2) = VesselOut(InsertPoint,2);
            else
                VesselOut(NoPoints+1,2) = VesselOut(InsertPoint,2)+10;
            end
            Intersections = cubeIntersect(VesselOut(NoPoints+1,3:5),VesselOut(VesselOut(NoPoints+1,7),3:5));
            for Loop = 1:size(Intersections,1)
                I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
                J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
                K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
                Bucket{I,J,K}(end+1) = NoPoints+1;
                %Bucket{I,J,K} = unique(Bucket{I,J,K});
%                 Bucket{I,J,K} = sort(Bucket{I,J,K});
%                 Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
            end
            % Update buckets. Adding part of the split line segment back
            % into the corresponding bucket.
            
        end
        % Creating new line segment from the originally generated point.
        
        NoPoints = NoPoints+1; % Increase vessel tree size by 1.
        
        VesselOut(NoPoints+1,1) = NoPoints+1;
        %VesselOut(NoPoints+1,2) = VesselOut(NoPoints,2);
        if VesselOut(NoPoints,2) > 9
            VesselOut(NoPoints+1,2) = VesselOut(NoPoints,2);
        else
            VesselOut(NoPoints+1,2) = VesselOut(NoPoints,2)+10;
        end
        VesselOut(NoPoints+1,3:5) = New;
        VesselOut(NoPoints+1,7) = NoPoints;
        NoPoints = NoPoints+1; % Increase vessel tree size by 1 (two points overall have been created).
        
        Intersections = cubeIntersect(VesselOut(NoPoints,3:5),VesselOut(VesselOut(NoPoints,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(Domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(Domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(Domain,3),Intersections(Loop,3)])]);
            Bucket{I,J,K}(end+1) = NoPoints;
            %Bucket{I,J,K} = unique(Bucket{I,J,K});
%             Bucket{I,J,K} = sort(Bucket{I,J,K});
%             Bucket{I,J,K} = Bucket{I,J,K}([true;diff(Bucket{I,J,K}(:))>0]);
        end
        % Update buckets. Adds this new connection to the corresponding
        % buckets.
        
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    %%%%%%%% Display Progress %%%%%%%%%%%%%%%%%%%%%%%
    if ~Check25 && N/NoIterations >= 0.25
        disp('25% Complete'); Check25 = true;
    elseif Check25 && ~Check50 && N/NoIterations >= 0.5
        disp('50% Complete'); Check50 = true;
    elseif Check25 && Check50 && ~Check75 && N/NoIterations >= 0.75
        disp('75% Complete'); Check75 = true;
    elseif Check25 && Check50 && Check75 && N==NoIterations
        disp('100% Complete');
    end
    % Perform percent completion checking and output.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Finalise Output %%%%%%%%%%%%%%%%%%%%%%%%%%
V_logical = logical(VesselOut(:,1));
VesselOut = VesselOut(V_logical,:);
% This deletes any unused rows from the initialisation of VesselOut
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end

function res = removeDuplicates(res,modifiedBuckets)
% if ~isempty(input)
%     res = sort(input);
%     res = res([true;diff(res(:))>0]);
% else
%     res = [];
% end

%modifiedBuckets = unique(modifiedBuckets, 'rows');


% Extract the size of input:
for counter=1:length(modifiedBuckets)
    ii = modifiedBuckets(counter,1);
    jj = modifiedBuckets(counter,2);
    kk = modifiedBuckets(counter,3);
            %if ~isempty(res{ii,jj,kk})
                res{ii,jj,kk} = sort(res{ii,jj,kk});
                res{ii,jj,kk} = res{ii,jj,kk}([true;diff(res{ii,jj,kk}(:))>0]);
            %end
end
end

