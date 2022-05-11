function exportPy(flows, filename)
%EXPORTPY Exports variables for Python
%   Detailed explanation goes here
Vessel1 = flows.arteries.tree;
Vessel2 = flows.veins.tree;
SplitArt = flows.SplitArt;
SplitVein = flows.SplitVein;
FdotArt = flows.FdotArt;
FdotVein = flows.FdotVein;
MdotArt = flows.MdotVesselArt;
MdotVein = flows.MdotVesselVein;
Mdot1 = flows.MdotVoxelsArt;
Mdot2 = flows.MdotVoxelsVein;
Vol1 = flows.geometry.art.V;
Vol2 = flows.geometry.vein.V;
L1 = flows.geometry.art.L;
L2 = flows.geometry.vein.L;
U = flows.porous.U;
V = flows.porous.V;
W = flows.porous.W;
GM_WM = flows.grey_white;
DomTot = flows.tissue;
InletPoints = flows.arteries.termPoints;
OutletPoints = flows.veins.termPoints;
InletFlows = flows.arteries.termFlows;
OutletFlows = flows.veins.termFlows;
Porosity = flows.porosity;
VoxelSize = flows.voxelSize;
MassFlow = flows.massFlow;
Rho_b = flows.bloodDensity;
Rho = flows.density;

save(filename,'-v7.3','-regexp','^(?!(flows)$).');
end

