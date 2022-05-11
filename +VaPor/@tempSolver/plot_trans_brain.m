function plot_trans_brain(obj, timestep, pos, climits , full_head)
%PLOTBRAIN Summary of this function goes here
%   Detailed explanation goes here
RowConvert = find(obj.flow_obj.tissue);
Tt = zeros(size(obj.flow_obj.tissue));
Tt(RowConvert) = obj.Transient_Temp.tissue(:,timestep);
if nargin < 5
    Tt(~obj.flow_obj.grey_white)=nan;
else
    Tt(~obj.flow_obj.tissue)=nan;
end
planecut('z', pos, Tt);
if nargin > 3
    if ~isempty(climits)
        caxis(climits);
    end
end
end