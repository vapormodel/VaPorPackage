function [obj] = solveit(obj, T_amb, T_blood)
%SOLVE Solves the temperatures
%   Solves the temperatures using BiCG.

if isempty(obj.T_Solve)
% Calculate the Graetz number and then apply relation for the Nusselt
% Number.
Gz1 = obj.flow_obj.bloodDensity*obj.bloodCp*(max(abs(obj.flow_obj.FdotArt),[],2)/obj.flow_obj.bloodDensity./(obj.flow_obj.geometry.art.Aavg)).*obj.flow_obj.arteries.tree(:,6).^2/obj.bloodKc./obj.flow_obj.geometry.art.L;
Gz2 = obj.flow_obj.bloodDensity*obj.bloodCp*(max(abs(obj.flow_obj.FdotVein),[],2)/obj.flow_obj.bloodDensity./(obj.flow_obj.geometry.vein.Aavg)).*obj.flow_obj.veins.tree(:,6).^2/obj.bloodKc./obj.flow_obj.geometry.vein.L;
obj.Nu.art = obj.Nu.base+0.155*exp(1.58*log10(Gz1));
obj.Nu.vein = obj.Nu.base+0.155*exp(1.58*log10(Gz2));

% Capped as the equation is only valid for Gz<1000
obj.Nu.art(obj.Nu.art>20) = 20;
obj.Nu.vein(obj.Nu.vein>20) = 20;

% Calculate Vessel Inter-Domain Heat Transfer Terms

% Arterial Segments
% obj.Beta.art = zeros(size(obj.flow_obj.arteries.tree,1),1);
% for n = 2:size(obj.flow_obj.arteries.tree,1)
%     if size(obj.flow_obj.Vessel2VolumeArt{n},1)~=0
%         obj.Beta.art(n) = (obj.Nu.art(n)*obj.bloodKc./obj.flow_obj.geometry.art.Davg(n)).*obj.flow_obj.geometry.art.AL(n) / size(obj.flow_obj.Vessel2VolumeArt{n},1);
%     else
%         obj.Beta.art(n) = 0;
%     end
% end
obj.Beta.art = beta_calc(obj.flow_obj.arteries.tree,obj.flow_obj.Vessel2VolumeArt,obj.Nu.art, obj.bloodKc,obj.flow_obj.geometry.art.Davg,obj.flow_obj.geometry.art.AL);

% Venous Segments
% obj.Beta.vein = zeros(size(obj.flow_obj.veins.tree,1),1);
% for n = 2:size(obj.flow_obj.veins.tree,1)
%     if size(obj.flow_obj.Vessel2VolumeVein{n},1)~=0
%         obj.Beta.vein(n) = (obj.Nu.vein(n)*obj.bloodKc./obj.flow_obj.geometry.vein.Davg(n)).*obj.flow_obj.geometry.vein.AL(n) / size(obj.flow_obj.Vessel2VolumeVein{n},1);
%     else
%         obj.Beta.vein(n)=0;
%     end
% end
obj.Beta.vein = beta_calc(obj.flow_obj.veins.tree,obj.flow_obj.Vessel2VolumeVein,obj.Nu.vein, obj.bloodKc,obj.flow_obj.geometry.vein.Davg,obj.flow_obj.geometry.vein.AL);
% Setting up inter-domain heat transfer terms for each venous line
% segment.

%matrix_builder(obj, T_amb, h_surf, T_blood, [T_blood, T_blood, T_blood]);
matrix_builder(obj);
%decompMat(obj);
[obj.luFac.L,obj.luFac.U] = ilu(obj.T_Solve);
end

bc_builder(obj, T_amb, T_blood, [T_blood, T_blood, T_blood]);

%T_N = obj.T_Solve\obj.D;
%T_N = divMat(obj);
[T_N,flag, relres, numIt] = bicg(obj.T_Solve, obj.D, 1e-10, 20000, obj.luFac.L, obj.luFac.U, 35.*ones(size(obj.D, 1), 1));
if flag ~= 0
                warning('VaPor:IterativeMethodUnconverged', 'The BICG solver has not converged to the specified tolerance.');
                fprintf('The residual was %.2e\n', relres);
            end
fprintf('Took %d iterations.\n', numIt);
NumDomTot = numel(obj.flow_obj.tissue(obj.flow_obj.tissue));
RowConvert = find(obj.flow_obj.tissue);  % goes from row number to I,J,K
NumDomRows = NumDomTot + obj.RowAdjustment; % Total rows assigned to voxels.
VesselRow = size(obj.flow_obj.arteries.tree,1);

Tt = zeros(size(obj.flow_obj.tissue));
Tt(RowConvert) = T_N(1:NumDomTot);

if obj.Beta.porous == Inf
    Tb = Tt;
else
    Tb = T_N(NumDomTot+1:NumDomRows);
end

T_Art = T_N(NumDomRows+1:NumDomRows+VesselRow);
T_Vein = T_N(NumDomRows+1+VesselRow:end);

Tt(~obj.flow_obj.tissue) = NaN;
Tb(~obj.flow_obj.tissue) = NaN;

% Export the results:
obj.temperatures.raw = T_N;
obj.temperatures.tissue = Tt;
obj.temperatures.art = T_Art;
obj.temperatures.vein = T_Vein;



end

