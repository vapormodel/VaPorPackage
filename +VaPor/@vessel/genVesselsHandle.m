function genVesselsHandle(obj)
%GENVESSELS Generates vessels DRYly
%   This function generates vessels based on the principle of not repeating
%   yourself. 
%   Author:         Luke Fulford - luke.fulford@ed.ac.uk
%   Last Modified:  29/01/2020

clear connectPointHandle

%% Allocate a matrix to hold the generated vessels: 
treeSave = obj.tree;
obj.growTree(obj.generateIterations);

%% Create the probability maps:
[probIndex, probDensity, probCumSum] =...
    genProb(obj.domainWeighting, obj.domainLimits, obj.weightFactor);

%% Bin the vessels into their voxels: 
bins = genBins(treeSave, obj.domainLimits);
% By 'binning' the vessels into their voxels we can more effectivley find
% the vessels that are local to a given point.

%% Choose the voxels that we want to visit in advance: 
voxelsToRecieve = ...
    randsample(length(probDensity), obj.generateIterations, true, probDensity)';
% Using rand sample is more efficient than selecting vessels one-by-one.
% However it does mean that the probability of a voxel getting a new
% termination is independant of the number of terminations it already has.

%% Loop through the new locations and add vessels:
progress = 0;
for thisVoxel = voxelsToRecieve
    progress = progress + 1; % Increment the loop-tracker.
    
    %% Pick a random point somewhere around this voxel:
    point = probIndex(thisVoxel,:) + [rand-0.5, rand-0.5, rand-0.5];
    
    %% Search for vessels that are 'close'
    closeVessels = binSearch(bins, point, obj.domainLimits, 1, obj.tree);
    
    %% Calculate how close the found points are. 
    distances = calcDists(point,closeVessels,obj.tree);
    
    %% Pick the point to use.
    [~, minVertex] = min(distances(:,1));
    % How this point is selected is one of the features that will differ
    % between generation methods - for RRT & RRT* it is always based on the
    % closest point. For some of the more exotic methods this will be
    % chosen based on some function. 
    
    closestType = distances(minVertex, 2);
    % This indicates which part of the line segment is closest - 0 or 1
    % indicates that it is nearest a node, whilst another value indicates
    % that we need to add a node. 
    
    minVertex = closeVessels(minVertex);
    % We need to unwrap the value of minVertex to recover the actual node
    % number. 
    
    otherVertex = obj.tree(minVertex, 7);
    % Find the parent of minVertex. 
    
    %% If we're using RRT* then the length of the segment is limited:
    if obj.useRRTstar
        if closestType==0
            PointOfConnection = obj.tree(minVertex,3:5);
        elseif closestType == 1
            PointOfConnection = obj.tree(obj.tree(minVertex,7),3:5);
        else
            PointOfConnection = obj.tree(minVertex,3:5)+closestType*(obj.tree(obj.tree(minVertex,7),3:5)-obj.tree(minVertex,3:5));
        end
        
        if norm(point - PointOfConnection) > obj.RRTstarEpsilon
            point = PointOfConnection + obj.RRTstarEpsilon/norm(point - PointOfConnection)*(point - PointOfConnection);
        end
    end
    
    %% Connect the point: 
    [bins] = connectPointHandle(obj, point, minVertex, otherVertex, closestType, bins);
    
end
%% Cleanup the tree: 
obj.tree = obj.tree(logical(obj.tree(:,1)),:);
end

