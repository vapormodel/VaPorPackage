function [obj] = loadVesselSWC(obj,pathway, postProcess)
%LOADVESSELS loads vessels from a SWC file and applies basic cleanup.

obj.tree = load(pathway);
if ~isempty(postProcess)
    obj.tree = postProcess(obj, obj.tree);
else
    obj.tree = autofit(obj.tree, obj.domainLimits);
end

end

