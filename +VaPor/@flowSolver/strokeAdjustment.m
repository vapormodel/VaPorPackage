function obj = strokeAdjustment(obj, reduction_frac)
disp('%%%%%%% Adjusting Flowrates from Stroke %%%%%%%%%%%') % Display
disp('Creating Network Trees')
tic

% Save the Old (pre-stroke) perfusion before alterations: 
obj.preStrokePerfusion = obj.measuredPerfusion;

numArt = size(obj.arteries.tree,1);
numVein = size(obj.veins.tree,1);
nxnynz = numel(obj.grey_white(obj.grey_white));
row_convert = find(obj.grey_white); % goes from row number back to obj.grey_white
row_convert2 = (1:nxnynz) + numArt;
GM_WM_convert = zeros(size(obj.grey_white));
GM_WM_convert(row_convert) = row_convert2;

error_amt = 1e-15;

Tree = cell(numArt+nxnynz+numVein,1);
% data stored in [location, flowrate to that location]
% for arteries and veins [location, flowrate to that location, vessel segment flow alteration refers to]
% for 3D flow [location, flowrate to that location, number 1 to 6 to reflect the directional flow]


% find downstream flows and weights for arteries
for art_loop = 1:numArt
    
    if obj.FdotArt(art_loop,1) <= 0
        % probablility of travelling down artery length
        if abs(obj.FdotArt(art_loop,1)) > 0
            Tree{art_loop}(end+1,:) = [obj.arteries.tree(art_loop,7),abs(obj.FdotArt(art_loop,1)),art_loop];
        end
        
        % probablility of transferring domain on artery length
        if ~isempty(obj.Vessel2VolumeArt{art_loop})
            for n = 1:size(obj.Vessel2VolumeArt{art_loop},1)
                if obj.Vessel2VolumeArt{art_loop}(n,4) ~= 0
                    Tree{art_loop}(end+1,:) = [GM_WM_convert(obj.Vessel2VolumeArt{art_loop}(n,1),obj.Vessel2VolumeArt{art_loop}(n,2),obj.Vessel2VolumeArt{art_loop}(n,3))...
                        ,obj.Vessel2VolumeArt{art_loop}(n,4),art_loop];
                end
            end
        end
        
    elseif obj.FdotArt(art_loop,2) >= 0
        
        % probablility of travelling down artery length
        if abs(obj.FdotArt(art_loop,2)) > 0 && abs(obj.FdotArt(art_loop,2)) > error_amt
            Tree{obj.arteries.tree(art_loop,7)}(end+1,:) = [art_loop,abs(obj.FdotArt(art_loop,2)),art_loop];
        end
        
        % probablility of transferring domain on artery length
        if ~isempty(obj.Vessel2VolumeArt{art_loop})
            for n = 1:size(obj.Vessel2VolumeArt{art_loop},1)
                if obj.Vessel2VolumeArt{art_loop}(n,4) ~= 0
                    Tree{obj.arteries.tree(art_loop,7)}(end+1,:) = [GM_WM_convert(obj.Vessel2VolumeArt{art_loop}(n,1),obj.Vessel2VolumeArt{art_loop}(n,2),obj.Vessel2VolumeArt{art_loop}(n,3))...
                        ,obj.Vessel2VolumeArt{art_loop}(n,4),art_loop];
                end
            end
        end
        
    end
    
end

% find downstream flows and weights for domain
for domain_loop = 1:nxnynz
    down_temp = [];
    
    [I,J,K] = ind2sub(size(obj.grey_white),row_convert(domain_loop));
    
    Ax = obj.voxelSize^2; Ay = obj.voxelSize^2; Az = obj.voxelSize^2;
    if obj.porous.U(I,J,K)<0, down_temp(end+1,:) = [GM_WM_convert(I-1,J,K),Ax*obj.bloodDensity*abs(obj.porous.U(I,J,K)),1]; end
    if obj.porous.U(I+1,J,K)>0, down_temp(end+1,:) = [GM_WM_convert(I+1,J,K),Ax*obj.bloodDensity*abs(obj.porous.U(I+1,J,K)),2]; end
    if obj.porous.V(I,J,K)<0, down_temp(end+1,:) = [GM_WM_convert(I,J-1,K),Ay*obj.bloodDensity*abs(obj.porous.V(I,J,K)),3]; end
    if obj.porous.V(I,J+1,K)>0, down_temp(end+1,:) = [GM_WM_convert(I,J+1,K),Ay*obj.bloodDensity*abs(obj.porous.V(I,J+1,K)),4]; end
    if obj.porous.W(I,J,K)<0, down_temp(end+1,:) = [GM_WM_convert(I,J,K-1),Az*obj.bloodDensity*abs(obj.porous.W(I,J,K)),5]; end
    if obj.porous.W(I,J,K+1)>0, down_temp(end+1,:) = [GM_WM_convert(I,J,K+1),Az*obj.bloodDensity*abs(obj.porous.W(I,J,K+1)),6]; end
    
    if ~isempty(obj.Volume2VesselVein{I,J,K})
        for n = 1:size(obj.Volume2VesselVein{I,J,K},1)
            feed = obj.Volume2VesselVein{I,J,K}(n,1);
            if obj.SplitVein(feed,1)~=0
                down_temp(end+1,:) = [feed+numArt+nxnynz,obj.SplitVein(n,2)*abs(obj.Volume2VesselVein{I,J,K}(n,2)),obj.Volume2VesselVein{I,J,K}(n,1)];
                down_temp(end+1,:) = [obj.veins.tree(feed,7)+numArt+nxnynz,(1-obj.SplitVein(n,2))*abs(obj.Volume2VesselVein{I,J,K}(n,2)),obj.Volume2VesselVein{I,J,K}(n,1)];
            else
                if obj.FdotVein(feed,2) <= 0
                    feed = obj.veins.tree(feed,7);
                end
                down_temp(end+1,:) = [feed+numArt+nxnynz,abs(obj.Volume2VesselVein{I,J,K}(n,2)),obj.Volume2VesselVein{I,J,K}(n,1)];
                
            end
        end
    end
    
    if ~isempty(down_temp)
        Tree{domain_loop+numArt} = [Tree{domain_loop+numArt};down_temp];
    end
    
end

% find downstream flows and weights for veins
for vein_loop = 1:numVein
    
    if obj.FdotVein(vein_loop,2) < 0
        Tree{vein_loop+numArt+nxnynz}(end+1,:) = [obj.veins.tree(vein_loop,7)+numArt+nxnynz,abs(obj.FdotVein(vein_loop,2)),vein_loop];
    elseif obj.FdotVein(vein_loop,1) > 0
        Tree{obj.veins.tree(vein_loop,7)+numArt+nxnynz}(end+1,:) = [vein_loop+numArt+nxnynz,abs(obj.FdotVein(vein_loop,1)),vein_loop];
    end
    
end
toc

% Do the backwards flow for stroke point in arterial domain only
tic
Tree2 = cell(numArt,1);

for art_loop = 2:numArt
    
    if obj.FdotArt(art_loop,1) < 0 && abs(obj.FdotArt(art_loop,1)) > error_amt
        %             Tree2{art_loop}(end+1,:) = [obj.arteries.tree(art_loop,7),abs(obj.FdotArt(art_loop,1)),art_loop];
        Tree2{obj.arteries.tree(art_loop,7)}(end+1,:) = [art_loop,abs(obj.FdotArt(art_loop,1)),art_loop];
        
    elseif obj.FdotArt(art_loop,2) > 0 && abs(obj.FdotArt(art_loop,2)) > error_amt
        %             Tree2{obj.arteries.tree(art_loop,7)}(end+1,:) = [art_loop,abs(obj.FdotArt(art_loop,2)),art_loop];
        Tree2{art_loop}(end+1,:) = [obj.arteries.tree(art_loop,7),abs(obj.FdotArt(art_loop,2)),art_loop];
        
    end
    
end
toc


disp('Adjusting Arterial Flowrates (Downstream)')
tic
% start_position = 1;
% start_position = inlet(2,1); % arterial node of stroke occurence
% start_position = 825; % left middle cerebral artery
start_position = obj.strokeLocation;
Destinations = Tree{start_position};
%reduction_frac = 0.5;
stroke_source = zeros(size(obj.grey_white));
Destinations(:,2) = Destinations(:,2) * reduction_frac;
while ~isempty(Destinations)
    
    Destinations_new = [];
    
    for n = 1:size(Destinations,1)
        if Destinations(n,1) <= numArt
            
            % Reducing flow down that stream
            if obj.FdotArt(Destinations(n,3),1) <= 0 && obj.FdotArt(Destinations(n,3),2) <= 0
                obj.FdotArt(Destinations(n,3),:) = obj.FdotArt(Destinations(n,3),:) + Destinations(n,2); % Upstream and Downstream
            else
                obj.FdotArt(Destinations(n,3),:) = obj.FdotArt(Destinations(n,3),:) - Destinations(n,2); % Upstream and Downstream
            end
            
            if ~isempty(Tree{Destinations(n,1)})
                Destination_alter = [Tree{Destinations(n,1)}(:,1),...
                    (Destinations(n,2)*Tree{Destinations(n,1)}(:,2)/sum(Tree{Destinations(n,1)}(:,2))),...
                    Tree{Destinations(n,1)}(:,3)];
                Destinations_new = [Destinations_new;Destination_alter];
            end
        else
            
            % Reducing flow down that stream
            if obj.FdotArt(Destinations(n,3),1) <= 0 && obj.FdotArt(Destinations(n,3),2) <= 0
                obj.FdotArt(Destinations(n,3),2) = obj.FdotArt(Destinations(n,3),2) + Destinations(n,2); % Upstream only
            else
                obj.FdotArt(Destinations(n,3),1) = obj.FdotArt(Destinations(n,3),1) - Destinations(n,2); % Upstream only
            end
            
            [I,J,K] = ind2sub(size(obj.grey_white),row_convert(Destinations(n,1)-numArt));
            
            index = find(obj.Volume2VesselArt{I,J,K}(:,1)==Destinations(n,3));
            obj.Volume2VesselArt{I,J,K}(index,2) = obj.Volume2VesselArt{I,J,K}(index,2) - Destinations(n,2);
            index = find(obj.Vessel2VolumeArt{Destinations(n,3)}(:,1) == I & obj.Vessel2VolumeArt{Destinations(n,3)}(:,2) == J & obj.Vessel2VolumeArt{Destinations(n,3)}(:,3) == K);
            obj.Vessel2VolumeArt{Destinations(n,3)}(index,4) = obj.Vessel2VolumeArt{Destinations(n,3)}(index,4) - Destinations(n,2);
            
            obj.MdotVoxelsArt(I,J,K) = obj.MdotVoxelsArt(I,J,K) - Destinations(n,2);
            
            stroke_source(I,J,K) = stroke_source(I,J,K) + Destinations(n,2);
        end
    end
    Destinations = Destinations_new;
    %     size(Destinations,1)
end

toc
disp('Adjusting Porous Flowrates (Downstream)')
tic
stroke_source_veins = zeros(size(obj.veins.tree,1),1);

stroke_source(stroke_source<error_amt) = 0;
stroke_source2 = zeros(size(stroke_source));
while sum(sum(sum(stroke_source))) > 0
    for I = 1:size(stroke_source,1)
        for J = 1:size(stroke_source,2)
            for K = 1:size(stroke_source,3)
                if stroke_source(I,J,K)
                    Destinations = Tree{GM_WM_convert(I,J,K)};
                    outflow = sum(Destinations(:,2));
                    frac = Destinations(:,2)/outflow;
                    for n = 1:size(Destinations,1)
                        if  Destinations(n,1) <= (numArt + nxnynz)
                            
                            [i2,j2,k2] = ind2sub(size(obj.grey_white),row_convert(Destinations(n,1)-numArt));
                            
                            if Destinations(n,3) == 1, obj.porous.U(I,J,K) = obj.porous.U(I,J,K) + 1/obj.bloodDensity*1/obj.voxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 2, obj.porous.U(I+1,J,K) = obj.porous.U(I+1,J,K) - 1/obj.bloodDensity*1/obj.voxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 3, obj.porous.V(I,J,K) = obj.porous.V(I,J,K) + 1/obj.bloodDensity*1/obj.voxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 4, obj.porous.V(I,J+1,K) = obj.porous.V(I,J+1,K) - 1/obj.bloodDensity*1/obj.voxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 5, obj.porous.W(I,J,K) = obj.porous.W(I,J,K) + 1/obj.bloodDensity*1/obj.voxelSize^2*frac(n)*stroke_source(I,J,K); end
                            if Destinations(n,3) == 6, obj.porous.W(I,J,K+1) = obj.porous.W(I,J,K+1) - 1/obj.bloodDensity*1/obj.voxelSize^2*frac(n)*stroke_source(I,J,K); end
                            
                            
                            stroke_source2(i2,j2,k2) = stroke_source2(i2,j2,k2) + frac(n)*stroke_source(I,J,K);
                        else
                            
                            if obj.FdotVein(Destinations(n,3),2) > 0
                                obj.FdotVein(Destinations(n,3),2) = obj.FdotVein(Destinations(n,3),2) - frac(n)*stroke_source(I,J,K);
                            else
                                obj.FdotVein(Destinations(n,3),1) = obj.FdotVein(Destinations(n,3),1) + frac(n)*stroke_source(I,J,K);
                            end
                            
                            
                            index = find(obj.Volume2VesselVein{I,J,K}(:,1)==Destinations(n,3));
                            obj.Volume2VesselVein{I,J,K}(index,2) = obj.Volume2VesselVein{I,J,K}(index,2) + frac(n)*stroke_source(I,J,K);
                            if ~isempty(obj.Vessel2VolumeVein{Destinations(n,3)})
                                index = find(obj.Vessel2VolumeVein{Destinations(n,3)}(:,1) == I & obj.Vessel2VolumeVein{Destinations(n,3)}(:,2) == J & obj.Vessel2VolumeVein{Destinations(n,3)}(:,3) == K);
                                obj.Vessel2VolumeVein{Destinations(n,3)}(index,4) = obj.Vessel2VolumeVein{Destinations(n,3)}(index,4) + frac(n)*stroke_source(I,J,K);
                            end
                            
                            
                            stroke_source_veins(Destinations(n,1)-numArt-nxnynz) = stroke_source_veins(Destinations(n,1)-numArt-nxnynz) + frac(n)*stroke_source(I,J,K);
                        end
                    end
                end
            end
        end
    end
    stroke_source = stroke_source2;
    stroke_source(stroke_source<error_amt) = 0;
    %     sum(sum(sum(stroke_source)))
    stroke_source2 = zeros(size(stroke_source));
end
toc

disp('Adjusting Venous Flowrates (Downstream)')
tic
stroke_source_veins(stroke_source_veins<error_amt) = 0;
stroke_source_veins2 = zeros(size(stroke_source_veins));
while sum(sum(sum(stroke_source_veins))) > 0
    for iv = 1:size(stroke_source_veins,1)
        if stroke_source_veins(iv)
            if ~ismember(iv,obj.veins.termPoints)
                Destinations = Tree{iv+nxnynz+numArt};
                if ~isempty(Destinations)
                    outflow = sum(Destinations(:,2));
                    frac = Destinations(:,2)/outflow;
                    for n = 1:size(Destinations,1)
                        
                        if obj.FdotVein(Destinations(n,3),1) <= 0 && obj.FdotVein(Destinations(n,3),2) <= 0
                            obj.FdotVein(Destinations(n,3),:) = obj.FdotVein(Destinations(n,3),:) + frac(n)*stroke_source_veins(iv); % Upstream and Downstream
                        else
                            obj.FdotVein(Destinations(n,3),:) = obj.FdotVein(Destinations(n,3),:) - frac(n)*stroke_source_veins(iv); % Upstream and Downstream
                        end
                        if Destinations(n,1)-numArt-nxnynz > 0
                            stroke_source_veins2(Destinations(n,1)-numArt-nxnynz) = stroke_source_veins2(Destinations(n,1)-numArt-nxnynz) + frac(n)*stroke_source_veins(iv);
                        end
                    end
                end
            end
        end
    end
    stroke_source_veins = stroke_source_veins2;
    stroke_source_veins(stroke_source_veins<error_amt) = 0;
    %     sum(stroke_source_veins)
    stroke_source_veins2 = zeros(size(stroke_source_veins));
end
toc


disp('Adjusting Arterial Flowrates (Upstream)')
tic
Destinations = Tree2{start_position};
Destinations(:,2) = Destinations(:,2) * reduction_frac;
while ~isempty(Destinations)
    Destinations_new = [];
    
    for n = 1:size(Destinations,1)
        
        % Reducing flow down that stream
        if -obj.FdotArt(Destinations(n,3),1) <= 0 && -obj.FdotArt(Destinations(n,3),2) <= 0
            obj.FdotArt(Destinations(n,3),:) = obj.FdotArt(Destinations(n,3),:) - Destinations(n,2); % Upstream and Downstream
        else
            obj.FdotArt(Destinations(n,3),:) = obj.FdotArt(Destinations(n,3),:) + Destinations(n,2); % Upstream and Downstream
        end
        
        if ~isempty(Tree2{Destinations(n,1)})
            Destination_alter = [Tree2{Destinations(n,1)}(:,1),...
                (Destinations(n,2)*Tree2{Destinations(n,1)}(:,2)/sum(Tree2{Destinations(n,1)}(:,2))),...
                Tree2{Destinations(n,1)}(:,3)];
            Destinations_new = [Destinations_new;Destination_alter];
        end
    end
    Destinations = Destinations_new;
    %     size(Destinations,1)
end

toc

% redo massFlowrate
measurePerfusion
obj.measuredPerfusion = MeasuredPerfusion;

% Update the term point values
obj.arteries.termFlows = abs(obj.FdotArt(obj.arteries.termPoints,1))';

disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%')
