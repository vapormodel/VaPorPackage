function [obj] = solve_trans(obj, T_amb, T_blood, time_to_sim, time_step, init_cond)
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
    matrix_builder_CN(obj, time_step);
    decompMat(obj);
    obj.save_timestep = time_step;
end

bc_builder(obj, T_amb, T_blood, [T_blood, T_blood, T_blood]);


if length(init_cond) == 1
    T_NOld = init_cond*ones(size(obj.D));
else
    % We are given the full inital condition - sanity check size!
    if length(init_cond) == length(obj.D)
        T_NOld = init_cond;
    else
        error('Dimention mismatch on inital condition!');
    end
end

NumDomTot = numel(obj.flow_obj.tissue(obj.flow_obj.tissue));
RowConvert = find(obj.flow_obj.tissue);  % goes from row number to I,J,K
NumDomRows = NumDomTot + obj.RowAdjustment; % Total rows assigned to voxels.
VesselRow = size(obj.flow_obj.arteries.tree,1);

number_of_time_steps = ceil(time_to_sim/time_step);
time_to_sim = number_of_time_steps*time_step;
SaveTime = time_step;

T_TransientStore = zeros(size(obj.D,1),round(number_of_time_steps*time_step/SaveTime)+1);
T_TransientStore(:,1) = T_NOld;

T_TransientAdjust = zeros(size(obj.D));

Vol = obj.flow_obj.voxelSize^3;

T_TransientAdjust(1:NumDomTot) = Vol*(1-obj.flow_obj.porosity(RowConvert)).*obj.flow_obj.density(RowConvert).*obj.Cp(RowConvert)/time_step;

if obj.Beta.porous == Inf
    T_TransientAdjust(1:NumDomTot) = T_TransientAdjust(1:NumDomTot) +...
        Vol*obj.flow_obj.porosity(RowConvert)*obj.flow_obj.bloodDensity*obj.bloodCp/time_step;
else
    T_TransientAdjust(NumDomTot+1:NumDomRows) = Vol*obj.flow_obj.porosity(RowConvert)*obj.flow_obj.bloodDensity.*obj.bloodCp/time_step;
end

T_TransientAdjust(NumDomRows+1:NumDomRows+VesselRow) = obj.flow_obj.geometry.art.V*obj.flow_obj.bloodDensity*obj.bloodCp/time_step;

T_TransientAdjust(NumDomRows+VesselRow+1:end) = obj.flow_obj.geometry.vein.V*obj.flow_obj.bloodDensity*obj.bloodCp/time_step;


T_TransientAdjust(NumDomRows + obj.flow_obj.arteries.termPoints) = 0;


%T_N = obj.T_Solve\obj.D;

D_Backup = obj.D;
for TimeNo = 1:number_of_time_steps
    
    disp(['time_step: ' num2str(TimeNo) '/' num2str(number_of_time_steps)])
    %             D = D - T_TransientAdjust.*T_NOld;
    T_TransientAdjust2 = obj.T_Solve*T_NOld;
    T_TransientAdjust2(NumDomRows + obj.flow_obj.arteries.termPoints) = 0;
    obj.D = D_Backup - 2*(T_TransientAdjust.*T_NOld) - T_TransientAdjust2;
    
    %disp(['Solving Matrix Inversion size: ' num2str(size(T_Solve,1)) ' by ' num2str(size(T_Solve,2))]) % Display.
    tic % Start timing.
    
    T_N = divMat(obj); % Solve the linear system.
    
    toc % Finish timing.
    %disp('Matrix Inversion Completed') % Display.
    
    %         D = D + T_TransientAdjust.*T_NOld;
    T_NOld = T_N;
    
    if rem(TimeNo*time_step,SaveTime) == 0
        T_TransientStore(:,round(TimeNo*time_step/SaveTime+1)) = T_N;
    end
    
    
    
    
end

obj.Transient_Temp.tissue = T_TransientStore(1:NumDomTot,:);
obj.Transient_Temp.art = T_TransientStore(NumDomRows+1:NumDomRows+VesselRow,:);
obj.Transient_Temp.vein = T_TransientStore(NumDomRows+1+VesselRow:end,:);



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
obj.temperatures.tissue = Tt;
obj.temperatures.art = T_Art;
obj.temperatures.vein = T_Vein;



end

