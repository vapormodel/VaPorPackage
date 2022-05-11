function newStrokeAdjustmentSimple(obj, removed_flow)
%NEWSTROKEADJUSTMENT Adjusts stroke flows.
%   Detailed explanation goes here

if isempty(obj.preStrokePerfusion)
    obj.preStrokePerfusion = obj.measuredPerfusion;
end

% Only find the terminations if we really need to - since seeking them
% takes a long while.
if isempty(obj.strokeTerminations)
    obj.strokeTerminations = findTerms(obj.strokeLocation, obj.arteries.tree);
end

strokeMask = false(size(obj.grey_white));
deductions = zeros(size(obj.grey_white));
for i = 1:length(obj.strokeTerminations)
    dest = obj.Vessel2VolumeArt{obj.strokeTerminations(i)};
    for j = 1:size(dest,1)
        strokeMask(dest(j,1),dest(j,2),dest(j,3)) = true;
        deductions(dest(j,1),dest(j,2),dest(j,3)) = deductions(dest(j,1),dest(j,2),dest(j,3)) + dest(j,4);
    end
end
obj.strokeMask = strokeMask;
obj.strokeMask(~obj.grey_white) = false;
obj.deductions = deductions;
clear strokeMask deductions

vein_flow_normal = sum(obj.MdotVoxelsVein(obj.strokeMask));

scaled_for_flow = 1 + ((removed_flow*sum(obj.deductions,'all'))./vein_flow_normal);

%adjustedflow = distance * (venous_flow_to_assign/total_weighting);

obj.MdotVoxelsVein(obj.strokeMask) = scaled_for_flow * obj.MdotVoxelsVein(obj.strokeMask);

obj.MdotVoxelsArt = obj.MdotVoxelsArt - removed_flow*obj.deductions;
%obj.MdotVoxelsVein = obj.MdotVoxelsVein + obj.deductions;

obj.MdotVoxelsOverall = obj.MdotVoxelsArt + obj.MdotVoxelsVein;

obj.solve3D;

% Calculate the reduction in the artierial tree:
obj.MdotVesselArt(obj.strokeTerminations) = 0;
obj.bloodFlow = sum(obj.MdotVesselArt);

% Volume2VesselArt and Vessel2VolumeArt need recalculating.
Volume2VesselArtTemp = obj.Volume2VesselArt;
Vessel2VolumeArtTemp = cell(size(obj.Vessel2VolumeArt));
MdotVesselArt = zeros(size(obj.MdotVesselArt));
for I = 1:size(Volume2VesselArtTemp, 1)
    for J = 1:size(Volume2VesselArtTemp, 2)
        for K = 1:size(Volume2VesselArtTemp, 3)
            if ~isempty(Volume2VesselArtTemp{I,J,K})
                if obj.MdotVoxelsArt(I,J,K) == 0
                    divisor = 1;
                else
                    divisor = sum(Volume2VesselArtTemp{I,J,K}(:,2));
                end
                Volume2VesselArtTemp{I,J,K}(:,2) = obj.MdotVoxelsArt(I,J,K) * Volume2VesselArtTemp{I,J,K}(:,2)/divisor;
                for idx = 1:size(Volume2VesselArtTemp{I,J,K},1)
                    Vessel2VolumeArtTemp{Volume2VesselArtTemp{I,J,K}(idx,1)} = [Vessel2VolumeArtTemp{Volume2VesselArtTemp{I,J,K}(idx,1)}; I, J, K, Volume2VesselArtTemp{I,J,K}(idx,2)];
                    if ~isnan(Volume2VesselArtTemp{I,J,K}(idx,2))
                        MdotVesselArt(Volume2VesselArtTemp{I,J,K}(idx,1)) = MdotVesselArt(Volume2VesselArtTemp{I,J,K}(idx,1)) + Volume2VesselArtTemp{I,J,K}(idx,2);
                    end
                end
            end
        end
    end
end

% Calculate the reduction in the venous tree:
Volume2VesselVeinTemp = obj.Volume2VesselVein;
Vessel2VolumeVeinTemp = cell(size(obj.Vessel2VolumeVein));
MdotVesselVein = zeros(size(obj.MdotVesselVein));
% Step over the cells:
for I = 1:size(Volume2VesselVeinTemp, 1)
    for J = 1:size(Volume2VesselVeinTemp, 2)
        for K = 1:size(Volume2VesselVeinTemp, 3)
            if ~isempty(Volume2VesselVeinTemp{I,J,K})
                if obj.MdotVoxelsVein(I,J,K) == 0 || sum(Volume2VesselVeinTemp{I,J,K}(:,2)) == 0
                    divisor = 1;
                else
                    divisor = sum(Volume2VesselVeinTemp{I,J,K}(:,2));
                end
                Volume2VesselVeinTemp{I,J,K}(:,2) = obj.MdotVoxelsVein(I,J,K) * Volume2VesselVeinTemp{I,J,K}(:,2)/divisor;
                if isnan(sum(Volume2VesselVeinTemp{I,J,K}(:,2)))
                    fprintf('Nan Found');
                end
                for idx = 1:size(Volume2VesselVeinTemp{I,J,K},1)
                    Vessel2VolumeVeinTemp{Volume2VesselVeinTemp{I,J,K}(idx,1)} = [Vessel2VolumeVeinTemp{Volume2VesselVeinTemp{I,J,K}(idx,1)}; I, J, K, Volume2VesselVeinTemp{I,J,K}(idx,2)];
                    if ~isnan(Volume2VesselVeinTemp{I,J,K}(idx,2))
                        MdotVesselVein(Volume2VesselVeinTemp{I,J,K}(idx,1)) = MdotVesselVein(Volume2VesselVeinTemp{I,J,K}(idx,1)) + Volume2VesselVeinTemp{I,J,K}(idx,2);
                    end
                end
            end
        end
    end
end
obj.Volume2VesselVein = Volume2VesselVeinTemp;
obj.Vessel2VolumeVein = Vessel2VolumeVeinTemp;
obj.Volume2VesselArt = Volume2VesselArtTemp;
obj.Vessel2VolumeArt = Vessel2VolumeArtTemp;
obj.MdotVesselVein = MdotVesselVein;
obj.solve1D

% FIN
end

