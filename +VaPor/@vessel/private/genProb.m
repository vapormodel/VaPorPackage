function [idx,density,cumSum] = genProb(probMap, domain, weightFactor)
%GENPROB Summary of this function goes here
%   Detailed explanation goes here

probMap(~domain) = 0; % Everything outside the desired domain has zero probablilty.
probMap(isnan(probMap)) = 0; % Remove any NaNs from the data and set them to zero.
% Initialising the vessel generation probablity map

boolMap = logical(probMap);
idx = zeros(numel(boolMap(boolMap)),3);
[idx(:,1),idx(:,2),idx(:,3)] = ind2sub(size(boolMap),find(boolMap));
% Indexing points for ease of retrieval.

density = (probMap(boolMap).^weightFactor)/(sum(probMap(boolMap).^weightFactor));
cumSum = cumsum(density);
% Creating cummulitive sum of probablity map.
end

