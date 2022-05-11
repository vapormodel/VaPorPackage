function [newTree,probCumSum] = genVessels(existingTree,probMap,domain,...
    numIterations,weightFactor,allowZero,joinRandom,confineJoin, regionGen, RRTstarLimit)
%GENVESSELS Generates vessels DRYly
%   This function generates vessels based on the principle of not repeating
%   yourself.
%   Author:         Luke Fulford - luke.fulford@ed.ac.uk
%   Last Modified:  03/02/2020

clear connectPoint

%% Housekeeping
% Here are admin steps.

% If we are using join random we need to search additional bins.
if joinRandom
    startDist = 4;
else
    startDist = 1;
end

% Parameters to calculate the number of nodes considered:
init_num = 20;
div_num = 40;

% For confine join we need to load some files:
if confineJoin || regionGen
    load('terrMask - simp expand.mat', 'regionMask');
    maskKey = [3 2 6 4 5 7];
    reverseKey = [0, 2, 1, 4, 5, 3, 6];
end
%TODO: move this to the constructor.


%% Allocate a matrix to hold the generated vessels:
newTree = mallocTree(existingTree, numIterations);

%% Create the probability maps:
[probIndex, probDensity, probCumSum] =...
    genProb(probMap, domain, weightFactor);

%% Bin the vessels into their voxels:
if ~confineJoin && ~regionGen
    bins = genBins(existingTree, domain, allowZero);
else
    bins = genBins(existingTree, domain, allowZero, confineJoin, reverseKey);
end
% By 'binning' the vessels into their voxels we can more effectivley find
% the vessels that are local to a given point.
reverse_mat = zeros(size(domain));
reverse_mat(:) = 1:numel(domain);
%% Choose the voxels that we want to visit in advance:
voxelsToRecieve = ...
    randsample(length(probDensity), numIterations, true, probDensity)';
% Using rand sample is more efficient than selecting vessels one-by-one.
% However it does mean that the probability of a voxel getting a new
% termination is independant of the number of terminations it already has.
%% Loop through the new locations and add vessels:
progress = 0;
for thisVoxel = voxelsToRecieve
    progress = progress + 1; % Increment the loop-tracker.
    if mod(progress, 1000) == 0
        fprintf('%d of %d (%.2f%%)\n', progress, numIterations, progress*100/numIterations);
    end
    %% Pick a random point somewhere around this voxel:
    point = probIndex(thisVoxel,:) + [rand-0.5, rand-0.5, rand-0.5];
    
    %% If joining is confined, choose a tree now:
    if confineJoin
        %attach_fam = randsample(6, 1, true, regionMask(probIndex(thisVoxel, 1), probIndex(thisVoxel, 2), probIndex(thisVoxel, 3), :));
        
        % Rather than using an individual voxel, use a combination of 2
        % voxels either side of the chosen one - i.e. a cube 5x5x5 with the
        % selected voxel in the center.
        expansion_amount = 2;
        lower_x = max(probIndex(thisVoxel,1) - expansion_amount, 1);
        upper_x = probIndex(thisVoxel,1) + expansion_amount;
        lower_y = max(probIndex(thisVoxel,2) - expansion_amount, 1);
        upper_y = probIndex(thisVoxel,2) + expansion_amount;
        lower_z = max(probIndex(thisVoxel,3) - expansion_amount, 1);
        upper_z = probIndex(thisVoxel,3) + expansion_amount;
        
        % Verify that the new limits are in the brain matter
        
        prob_mat = regionMask(lower_x:upper_x, lower_y:upper_y, lower_z:upper_z, :);
        prob_vec = squeeze(sum(sum(sum(prob_mat))))./sum(prob_mat(:));
        attach_fam = randsample(6, 1, true, prob_vec);
    elseif regionGen
        attach_fam = randsample(6, 1, true, [0.4, 0.05, 0.05, 0.4, 0.05, 0.05]);
    end
    %% Search for vessels that are 'close'
    if confineJoin
        %         closeVessels = binSearch(bins(:,:,:,attach_fam), point, domain, startDist, newTree);
        %closeVessels = binSearch(bins(1+(attach_fam-1)*numel(domain):(attach_fam)*numel(domain),:), point, domain, startDist, newTree);
        closeVessels = binSearch(bins, point, domain, startDist, newTree, attach_fam, reverse_mat);
    else
        closeVessels = binSearch(bins, point, domain, startDist, newTree, 1, reverse_mat);
    end
    
    %% Calculate how close the found points are.
    distances = calcDists(point,closeVessels,newTree);
    %% Check if we are using join random:
    if joinRandom
        %% Calculate the number of nodes to consider:
        num_nodes =...
            max(init_num, ceil(div_num*length(distances(:,1))/100));
        
        %% Find num_nodes shortest distances:
        [Dist_Short, Dist_Short_Idx] = mink(distances(:,1),num_nodes);
        
        %% Randomly select a lucky node:
        minVertex = randsample(Dist_Short_Idx, 1, true, (1./distances(Dist_Short_Idx).^1)./(sum(1./Dist_Short.^1)));
    else
        %% Pick the point to use.
        [~, minVertex] = min(distances(:,1));
        % How this point is selected is one of the features that will differ
        % between generation methods - for RRT & RRT* it is always based on the
        % closest point. For some of the more exotic methods this will be
        % chosen based on some function.
    end
    
    closestType = distances(minVertex, 2);
    % This indicates which part of the line segment is closest - 0 or 1
    % indicates that it is nearest a node, whilst another value indicates
    % that we need to add a node.
    
    minVertex = closeVessels(minVertex);
    % We need to unwrap the value of minVertex to recover the actual node
    % number.
    
    otherVertex = newTree(minVertex, 7);
    % Find the parent of minVertex.
    
    %% Only bifurications are permissable.
    % So if a node already has more than 2 children we'll need to join in
    % the middle.
    if otherVertex ~= -1
        if closestType == 0
            if sum(newTree(:,7)==minVertex) > 1
                closestType = 0.5;
            end
        elseif closestType == 1
            if sum(newTree(:,7)==newTree(minVertex,7)) > 1
                closestType = 0.5;
            end
        end
    end
    
    %% If we're using RRT* then the length of the segment is limited:
    if nargin == 10
        if closestType==0
            PointOfConnection = newTree(minVertex,3:5);
        elseif closestType == 1
            PointOfConnection = newTree(newTree(minVertex,7),3:5);
        else
            PointOfConnection = newTree(minVertex,3:5)+closestType*(newTree(newTree(minVertex,7),3:5)-newTree(minVertex,3:5));
        end
        
        if norm(point - PointOfConnection) > RRTstarLimit
            point = PointOfConnection + RRTstarLimit/norm(point - PointOfConnection)*(point - PointOfConnection);
        end
    end
    
    %% Connect the point:
    if confineJoin
        %squeezed_bins = bins(:,:,:,attach_fam);
        %squeezed_bins = reshape(bins(attach_fam,:,:,:), size(bins(attach_fam,:,:,:), [2:4]));
        %         [newTree, bins(:,:,:,attach_fam)] = connectPoint(newTree, domain, point, minVertex, otherVertex, closestType, bins(:,:,:,attach_fam));
        %[newTree, bins(1+(attach_fam-1)*numel(domain):(attach_fam)*numel(domain),:)] = connectPoint(newTree, domain, point, minVertex, otherVertex, closestType, bins(1+(attach_fam-1)*numel(domain):(attach_fam)*numel(domain),:));
    [newTree, bins] = connectPoint(newTree, domain, point, minVertex, otherVertex, closestType, bins, attach_fam);
    else
        [newTree, bins] = connectPoint(newTree, domain, point, minVertex, otherVertex, closestType, bins, 1);
    end
    
end
%% Cleanup the tree:
newTree = newTree(logical(newTree(:,1)),:);
end

