function plotbrain(obj, pos)
%PLOTBRAIN Summary of this function goes here
%   Detailed explanation goes here
to_plot = obj.temperatures.tissue;
to_plot(~obj.flow_obj.grey_white)=nan;
planecut('z', pos, to_plot);
end

