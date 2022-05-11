function [bucketList] = binSearch(bins,point,domain,initRange,vesselTree, attach_fam, reverse_mat)
%BINSEARCH Searches bins for vessels close to the desired point.
%   Detailed explanation goes here


adjustor_val = (attach_fam-1)*numel(domain);


bucketList = []; % Initialise list of closest segments.
range = initRange; % Distance from original bucket to search from.
check = 0; % Variable for breaking loop.
while check == 0
    x_Lim = [max(1,floor(point(1))-range),min(size(domain,1),floor(point(1))+range)];
    y_Lim = [max(1,floor(point(2))-range),min(size(domain,2),floor(point(2))+range)];
    z_Lim = [max(1,floor(point(3))-range),min(size(domain,3),floor(point(3))+range)];
    % Expand the limits:
    %x_Lim_exp = x_Lim(1):1:x_Lim(2);
    %y_Lim_exp = y_Lim(1):1:y_Lim(2);
    %z_Lim_exp = z_Lim(1):1:z_Lim(2);
    %searchArea = bins(x_Lim(1):x_Lim(2), y_Lim(1):y_Lim(2), z_Lim(1):z_Lim(2));
    %searchArea = searchArea(~cellfun('isempty',searchArea));
    %listLength = sum(cellfun('length',searchArea),1);
    %ind_to_check = sub2ind(size(domain),x_Lim_exp,y_Lim_exp,z_Lim_exp);
    %bucketList = [];
    bucketList = zeros(1,10000);
    bucket_pos = 1;
%     for marker = 1:length(ind_to_check)
%         if bins(ind_to_check(marker),1) ~= 1
%             bucketList = [bucketList, bins(ind_to_check(marker),2:bins(ind_to_check(marker),1))];
%         end
%     end
for x_val = x_Lim(1):1:x_Lim(2)
    for y_val = y_Lim(1):1:y_Lim(2)
        for z_val = z_Lim(1):1:z_Lim(2)
            ind_to_check = reverse_mat(x_val,y_val,z_val);
            if bins(ind_to_check+adjustor_val,1) ~= 1
                %bucketList = [bucketList, bins(ind_to_check+adjustor_val,2:bins(ind_to_check+adjustor_val,1))];
                num_to_ins = bins(ind_to_check+adjustor_val,1)-1;
                bucketList(bucket_pos:bucket_pos+num_to_ins-1) = bins(ind_to_check+adjustor_val,2:bins(ind_to_check+adjustor_val,1));
                bucket_pos = num_to_ins + bucket_pos;
            end
        end
    end
end
if bucket_pos == 1
    bucketList = [];
else
    bucketList = bucketList(1:bucket_pos-1);
end
%     if listLength ~= 0
%         bucketList = zeros(listLength,1);
%         bucketPos = 1;
%         for counter = 1:length(searchArea)
%             numInsert = numel(searchArea{counter});
%             bucketList(bucketPos:bucketPos+numInsert-1) = searchArea{counter};
%             bucketPos = bucketPos + numInsert;
%         end
%        %bucketList = horzcat(searchArea{:})';
%     else
%         bucketList = [];
%     end
    
    
    if ~isempty(bucketList)
        check = 1; % If list is not empty, break the search.
    else
        range = range + 1; % Increase search range if list is empty.
    end
    % Checking if any vessels have been found.
    
    if range > max(size(domain))
        bucketList = vesselTree(:,1)';
        check = 1;
        % If no connections are found within the domain. Search the
        % whole tree (as it lies fully outside the domain). Very
        % inefficient if a large amount of vessels lie outside of the
        % domain for generation.
    end
    
end

bucketList = sort(bucketList);
bucketList = bucketList([true;diff(bucketList(:))>0]);
bucketList = bucketList(bucketList(:)>(max(bucketList)-100e3));
bucketList = bucketList';
end

