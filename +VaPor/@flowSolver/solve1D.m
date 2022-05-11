function [obj] = solve1D(obj)
%SOLVE1D Summary of this function goes here
%   Detailed explanation goes here

% Adjust the termial flows.
obj.arteries.termFlows = obj.arteries.termFlows./sum(obj.arteries.termFlows) * obj.bloodFlow;
obj.veins.termFlows = obj.veins.termFlows/sum(obj.veins.termFlows) * obj.bloodFlow;

% Call the direct solver for the arterial tree:
[obj.FdotArt, obj.SplitArt] = DirectSolve1D(...
    obj.arteries.tree,obj.MdotVesselArt,obj.arteries.termPoints,obj.arteries.termFlows,[],[]);

% Call direct solver for the venous tree:
[obj.FdotVein, obj.SplitVein] = DirectSolve1D(...
    obj.veins.tree,obj.MdotVesselVein,[],[],obj.veins.termPoints,-obj.veins.termFlows);

VesselDiameterAdjust = 1;
VesselDiameterToFlowAdjust = 1;

% Resolve for the new diameters and then update the geometry.
Dart = VesselDiameterAdjust*0.0332*abs(VesselDiameterToFlowAdjust*obj.FdotArt(:,2)).^0.3703; % Diameter Flowrate correlation.
Dart(Dart<=1e-5) = 1e-5; % Set minimum diameter as 10 micron.
% Artery diameters based on flowrates. Any diameter less than 10 micron
% is then set to 10 micron. This avoids any zero volume line segments at
% branch terminations.

Dvein = VesselDiameterAdjust*0.0332*abs(VesselDiameterToFlowAdjust*obj.FdotVein(:,2)).^0.3703; % Diameter Flowrate correlation.
Dvein(Dvein<=1e-5) = 1e-5; % Set minimum diameter as 10 micron.
% Venous diameters based on flowrates. Any diameter less than 10 micron
% is then set to 10 micron. This avoids any zero volume line segments at
% branch terminations.

obj.arteries.tree(:,6) = Dart;
obj.veins.tree(:,6) = Dvein;

[obj.geometry.art.L, obj.geometry.art.AL, obj.geometry.art.V, obj.geometry.art.Davg, obj.geometry.art.Aavg] = vesselGeometry(obj.arteries.tree,obj.voxelSize);
[obj.geometry.vein.L, obj.geometry.vein.AL, obj.geometry.vein.V, obj.geometry.vein.Davg, obj.geometry.vein.Aavg] = vesselGeometry(obj.veins.tree,obj.voxelSize);

end

