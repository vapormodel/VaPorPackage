function [MdotVessel,Vessel2Volume,Volume2Vessel,MdotVoxels] = vesselDomainExchangeFunc(Vessel,L,DomainFull,DomainSource,BloodFlow,BranchTerminationOption)


%%%%%% Create Inter-Domain Cells %%%%%%%%%%%%%%%%
Vessel2Volume = cell(size(Vessel,1),1);
Volume2Vessel = cell(size(DomainFull));
% Create cells that store information about inter-domain transfer. For each
% voxel, Volume2Vessel contains information about intersecting vessel
% segments for each voxel and how much transfer occurs between them.
% Vessel2Volume contains information about intersecting voxels for each
% vessel segment and how much transfer occurs between them.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Initial Guess For Inter-Domain Transfer %%
MdotVessel = BloodFlow * L/sum(L);
% Inter-Domain transfer based on proportional length of segment. This is
% subject to alteration further on if certain segments do not have voxel
% intersections.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Find Intersections for Segments %%%%%%%%%%
for iv = 2:size(Vessel,1)
    if Vessel(iv,1) == 0
        continue
    end
    Intersections = cubeIntersect(Vessel(iv,3:5),Vessel(Vessel(iv,7),3:5));
    % Find intersections between vessels and cubes.
    
    if ~isempty(Intersections)
        IntersectionsBool = true(size(Intersections,1),1);
        for iw = 1:size(Intersections,1)
            I = Intersections(iw,1); J = Intersections(iw,2); K = Intersections(iw,3);
            if I>=1 && I<=size(DomainFull,1) && J>=1 && J<=size(DomainFull,2) && K>=1 && K<=size(DomainFull,3)
                if ~DomainFull(Intersections(iw,1),Intersections(iw,2),Intersections(iw,3))
                    IntersectionsBool(iw) = false;
                end
            else
                IntersectionsBool(iw) = false;
            end
        end
        Intersections = Intersections(IntersectionsBool,:);
    end
    % Deleting intersections not within domain.
    
    
    if isempty(Intersections) && (round(Vessel(iv,3))>=1 && round(Vessel(iv,3))<=size(DomainFull,1))...
            && (round(Vessel(iv,4))>=1 && round(Vessel(iv,4))<=size(DomainFull,2))...
            && (round(Vessel(iv,5))>=1 && round(Vessel(iv,5))<=size(DomainFull,3))
        if DomainFull(round(Vessel(iv,3)),round(Vessel(iv,4)),round(Vessel(iv,5)))
            Intersections = [round(Vessel(iv,3)),round(Vessel(iv,4)),round(Vessel(iv,5))];
        end
    end
    % Check if vessel is within a cube (but doesn't intersect the boundries
    % of that cube).
    
    if ~isempty(Intersections)
        Intersections(:,end+1) = 0;
    end
    % Pre-allocate source term to any intersections found.
    
    Vessel2Volume{iv} = Intersections;
    % Update cell.
end
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Branch Termination List %%%%%%%%%%%%%%%%%%
BranchTerminationList = ones(size(Vessel,1),1);
if BranchTerminationOption
    for n=2:numel(BranchTerminationList)
        if Vessel(n,1) == 0
            continue
        end
        BranchTerminationList(Vessel(n,7)) = 0;
    end
end
% Checks whether branches have connections. If they do then they are
% not considered terminations. Only branches with terminations are allowed
% mass transfer. If the Option_BranchTerminationsOnly = [false] then all 
% segments are considered 'branch terminations' for the purposes of the
% solver.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Voxel Inter-Domain Transfer %%%%%%%%%%%%%%
MdotVoxels = zeros(size(DomainFull));
% Preallocate array for inter-domain mass transfer in voxels.

for iv = 2:size(Vessel,1) % For every line segment.
    if BranchTerminationList(iv) ~= 0
        if ~isempty(Vessel2Volume{iv}) % If there is an intersection.
            for  iw=1:size(Vessel2Volume{iv},1)
                I = Vessel2Volume{iv}(iw,1);
                J = Vessel2Volume{iv}(iw,2);
                K = Vessel2Volume{iv}(iw,3); % Get location of intersection voxel.
                if DomainSource(I,J,K) % If intersection voxel is in domain that allows inter-domain transfer.
                    Vessel2Volume{iv}(iw,4) = MdotVessel(iv)/numel(Vessel2Volume{iv}(iw,:)); % Update Vessel2Volume to include this specific transfer.
                    MdotVoxels(I,J,K) = MdotVoxels(I,J,K) + MdotVessel(iv)/numel(Vessel2Volume{iv}(iw,:)); % Update total transfer to voxel.
                end
            end
        end
    end
end
% This allocates the transfers from each segement to the corresponding
% voxels. If there are any segments that do not transfer then the result
% will be unbalanced. This is then corrected in the subsequent section.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Adjust Inter-Domain Transfer %%%%%%%%%%%%%
if sum(sum(sum(MdotVoxels))) == 0
    error('Error: No intersections between vessel tree and domain found.')
    % Check if there is any transfer to any voxel. Otherwise domains are
    % not connected and program cannot proceed.
else
    FracSource = BloodFlow/sum(sum(sum(MdotVoxels))); % Amount MdotVessels need to be adjusted to balance the inter-domain transfer.
    MdotVoxels = MdotVoxels * FracSource; % Adjusting the voxel inter-domain transfer to match overall BloodFlow.
    for iv = 1:size(Vessel,1) % For every line segment.
        if ~isempty(Vessel2Volume{iv}) % If there is an intersection.
            Vessel2Volume{iv}(:,4) = Vessel2Volume{iv}(:,4)*FracSource; % Adjusting the segment inter-domain transfer so that overall mass transfer is conserved.
            MdotVessel(iv) = sum(Vessel2Volume{iv}(:,4)); % Overall segment transfer is sum of flow to each voxel.
            for iw = 1:size(Vessel2Volume{iv},1) % For every voxel that this segment intersects.
                ijk = Vessel2Volume{iv}(iw,1:3);
                Volume2Vessel{ijk(1),ijk(2),ijk(3)}(end+1,:) = [iv,Vessel2Volume{iv}(iw,4)]; % Update stored value for inter-domain transfer in Volume2Vessel.
            end
        else
            MdotVessel(iv) = 0; % If there are no intersections segment inter-domain flow is zero.
        end
    end
end
% As some vessels do not transfer blood to their corresponding
% intersections, the overall transfer becomes unbalanced. This part fixes
% the imbalance by setting transfer in appropriate segments to zero and
% proportionally increasing inter-domain flow in the segments that do
% transfer flow.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

end

function Intersections = cubeIntersect(Pt1,Pt2)
% Checks whether a line defined by two points intersects cubes.

Intersections = zeros(50,3);
% Initialise output.
interCounter = 1;

Dvec = Pt2-Pt1;
% Distance between two points.

% Creating Bounding Box for Vessel
Imin = floor(min([Pt1(1),Pt2(1)]));
Imax = ceil(max([Pt1(1),Pt2(1)]));
Jmin = floor(min([Pt1(2),Pt2(2)]));
Jmax = ceil(max([Pt1(2),Pt2(2)]));
Kmin = floor(min([Pt1(3),Pt2(3)]));
Kmax = ceil(max([Pt1(3),Pt2(3)]));
% Limits the search of intersections to only those within the x, y
% and z limits of the two points.

for I2=Imin:Imax % Bounding X region of cylinder.
    for J2=Jmin:Jmax % Bounding Y region of cylinder.
        for K2=Kmin:Kmax % Bounding Z region of cylinder.
            
            boundbox_min = [I2-0.5,J2-0.5,K2-0.5];
            boundbox_max = [I2+0.5,J2+0.5,K2+0.5];
            % Values of cubes are assumed to be at centre so cube
            % walls are plus and minus 0.5 away.
            
            Intersect = false;
            % Initialise Intersect.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_min(1)-Pt1(1))/Dvec(1);
                if t2>0 && t2<1
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses lower x boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_min(2)-Pt1(2))/Dvec(2);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses lower y boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_min(3)-Pt1(3))/Dvec(3);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses lower z boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_max(1)-Pt1(1))/Dvec(1);
                if t2>0 && t2<1
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses upper x boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_max(2)-Pt1(2))/Dvec(2);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Z_Int = Pt1(3) + t2*Dvec(3);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Z_Int >= boundbox_min(3) && Z_Int <= boundbox_max(3)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses upper y boundary.
            
            if ~Intersect % Skips if already found an intersection.
                t2 = (boundbox_max(3)-Pt1(3))/Dvec(3);
                if t2>0 && t2<1
                    X_Int = Pt1(1) + t2*Dvec(1);
                    Y_Int = Pt1(2) + t2*Dvec(2);
                    if X_Int >= boundbox_min(1) && X_Int <= boundbox_max(1) && Y_Int >= boundbox_min(2) && Y_Int <= boundbox_max(2)
                        Intersect = true;
                    end
                end
            end
            % Check if line crosses upper z boundary.
            
            if Intersect
                Intersections(interCounter,:) = [I2,J2,K2];
                % If the vessel intersected any of the boundries
                % then it intersects the voxel and it is added to
                % the list of intersections.
                interCounter = interCounter+1;
            end
        end
    end
end

%IntersectionsMask = logical(Intersections);
%Intersections = Intersections(IntersectionsMask(1,:),:);
Intersections = Intersections(1:interCounter-1,:);
end