function [distances] = calcDists(point,bins,vesselTree)
%CALCDISTS Summary of this function goes here
%   Detailed explanation goes here
distances = zeros(size(bins,1),2);
param1s = zeros(size(bins,1),3);
param2s = zeros(size(bins,1),3);
param2list = zeros(size(bins,1),1);
paramExist = false(size(bins,1),3);

for M2 = 1:size(bins,1)
    M = bins(M2);
    if M == 1
        distances(M2,:) = [norm(point-vesselTree(1,3:5)),0];
    else
        vesselOutVector = vesselTree(M,:);
        param1s(M2,:) = vesselOutVector(3:5);
        param2list(M2) = vesselOutVector(7);
        paramExist(M2,:) = [true, true, true];
    end
end
param2s(paramExist) = vesselTree(param2list(paramExist(:,1)), 3:5);

distances(paramExist(:,1),:) = distPointToLineSegMat(param1s(paramExist(:, 1),:),param2s(paramExist(:,1),:),point);


%[D,T] = distPointToLineSeg(param1,param2,point);
%        distances(M2,:) = [D,T];
end


% for M2 = 1:size(bins,1)
%     M = bins(M2);
%     if M == 1
%         distances(M2,:) = [norm(point-vesselTree(1,3:5)),0];
%     else
%         vesselOutVector = vesselTree(M,:);
%         param1 = vesselOutVector(3:5);
%         param2 = vesselTree(vesselOutVector(7),3:5);
%         [D,T] = distPointToLineSeg(param1,param2,point);
%         distances(M2,:) = [D,T];
%     end
% end
% end

