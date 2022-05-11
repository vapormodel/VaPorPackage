function obj = vesselDomainExchange(obj)
%VESSELDOMAINEXCHANGE Summary of this function goes here
%   Detailed explanation goes here

[L1, ~, Vol1, ~, ~] = vesselGeometry(obj.arteries.tree,obj.voxelSize);

obj.geometry.art.L = L1;
obj.geometry.art.V = Vol1;

[MdotVessel,Vessel2Volume,Volume2Vessel,MdotVoxels] = ...
    vesselDomainExchangeFunc(obj.arteries.tree,L1,obj.tissue,obj.grey_white,obj.bloodFlow,true);

obj.MdotVesselArt = MdotVessel;
obj.Vessel2VolumeArt = Vessel2Volume;
obj.Volume2VesselArt = Volume2Vessel;
obj.MdotVoxelsArt = MdotVoxels;

[L1, ~, Vol1, ~, ~] = vesselGeometry(obj.veins.tree,obj.voxelSize);

obj.geometry.vein.L = L1;
obj.geometry.vein.V = Vol1;

[MdotVessel,Vessel2Volume,Volume2Vessel,MdotVoxels] = ...
    vesselDomainExchangeFunc(obj.veins.tree,L1,obj.tissue,obj.grey_white,-obj.bloodFlow,true);

obj.MdotVesselVein = MdotVessel;
obj.Vessel2VolumeVein = Vessel2Volume;
obj.Volume2VesselVein = Volume2Vessel;
obj.MdotVoxelsVein = MdotVoxels;

if obj.options.counterCurrentFlow
    obj.MdotVoxelsOverall = MdotVoxelsArt-(Perfusion*VoxelSize^3);
    obj.MdotVoxelsOverallB = (Perfusion*VoxelSize^3)+Mdot2;
else
    obj.MdotVoxelsOverall = obj.MdotVoxelsArt+obj.MdotVoxelsVein;
end

end

