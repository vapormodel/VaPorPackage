function [beta] = beta_calc(tree,vessel2vol,Nu, Kc,Davg,AL)
%BETA_CALC Summary of this function goes here
%   Detailed explanation goes here
beta = zeros(size(tree,1),1);
for n = 2:size(tree,1)
    if size(vessel2vol{n},1)~=0
        beta(n) = (Nu(n)*Kc./Davg(n)).*AL(n) / size(vessel2vol{n},1);
    else
        beta(n) = 0;
    end
end
end

