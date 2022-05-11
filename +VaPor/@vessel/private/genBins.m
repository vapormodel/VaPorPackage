function [bins] = genBins(vessels, domain, allowZero, confineJoin, reverseKey)
%GENBINS Summary of this function goes here
%   Detailed explanation goes here


if nargin < 4 || ~confineJoin
%     bins = cell(size(domain));
%     bins{max([1,min([size(domain,1),floor(vessels(1,3))])]),max([1,min([size(domain,2),floor(vessels(1,4))])]),max([1,min([size(domain,3),floor(vessels(1,5))])])} = 1;
%     for M = 2:size(vessels,1)
%         if vessels(M,2) ~= 0 || allowZero
%             intersections = cubeIntersect(vessels(M,3:5),vessels(vessels(M,7),3:5));
%             for Loop = 1:size(intersections,1)
%                 I = max([1,min([size(domain,1),intersections(Loop,1)])]);
%                 J = max([1,min([size(domain,2),intersections(Loop,2)])]);
%                 K = max([1,min([size(domain,3),intersections(Loop,3)])]);
%                 bins{I,J,K}(end+1) = M;
%                 if ~isempty(bins{I,J,K})
%                     bins{I,J,K} = sort(bins{I,J,K});
%                     bins{I,J,K} = bins{I,J,K}([true;diff(bins{I,J,K}(:))>0]);
%                 end
%             end
%         end
%     end
bins = zeros(numel(domain)*6,200);
    bins(:,1) = 1;
    for M = 2:size(vessels,1)
        Intersections = cubeIntersect(vessels(M,3:5),vessels(vessels(M,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
            if vessels(M, 2) ~= 0 && vessels(M, 2) ~= 10
                if vessels(M,2) < 10
                    idx_val = vessels(M,2);
                else
                    idx_val = vessels(M,2) - 10;
                end
                
                % bins{I,J,K,Family}(end+1) = M;
                % Update insert pos:
                bins(sub2ind(size(domain),I,J,K),1) = bins(sub2ind(size(domain),I,J,K),1) + 1;
                bins(sub2ind(size(domain),I,J,K),bins(sub2ind(size(domain),I,J,K),1)) = M;
%                 if ~isempty(bins{I,J,K,Family})
%                     bins{I,J,K,Family} = sort(bins{I,J,K,Family});
%                     bins{I,J,K,Family} = bins{I,J,K,Family}([true;diff(bins{I,J,K,Family}(:))>0]);
%                 end
                 if bins(sub2ind(size(domain),I,J,K),1) > 2
                    bins(sub2ind(size(domain),I,J,K),2:bins(sub2ind(size(domain),I,J,K))) = sort(bins(sub2ind(size(domain),I,J,K),2:bins(sub2ind(size(domain),I,J,K),1)));
                    %bins{I,J,K,Family} = bins{I,J,K,Family}([true;diff(bins{I,J,K,Family}(:))>0]);
                end
            end
        end
    end
else
%     bins = cell([size(domain),6]);
%     bins{max([1,min([size(domain,1),floor(vessels(1,3))])]),max([1,min([size(domain,2),floor(vessels(1,4))])]),max([1,min([size(domain,3),floor(vessels(1,5))])]),1} = 1;
%     bins{max([1,min([size(domain,1),floor(vessels(1,3))])]),max([1,min([size(domain,2),floor(vessels(1,4))])]),max([1,min([size(domain,3),floor(vessels(1,5))])]),2} = 1;
%     bins{max([1,min([size(domain,1),floor(vessels(1,3))])]),max([1,min([size(domain,2),floor(vessels(1,4))])]),max([1,min([size(domain,3),floor(vessels(1,5))])]),3} = 1;
%     bins{max([1,min([size(domain,1),floor(vessels(1,3))])]),max([1,min([size(domain,2),floor(vessels(1,4))])]),max([1,min([size(domain,3),floor(vessels(1,5))])]),4} = 1;
%     bins{max([1,min([size(domain,1),floor(vessels(1,3))])]),max([1,min([size(domain,2),floor(vessels(1,4))])]),max([1,min([size(domain,3),floor(vessels(1,5))])]),5} = 1;
%     bins{max([1,min([size(domain,1),floor(vessels(1,3))])]),max([1,min([size(domain,2),floor(vessels(1,4))])]),max([1,min([size(domain,3),floor(vessels(1,5))])]),6} = 1;
    %bins = spalloc(numel(domain)*6,4000,10e6);
    bins = zeros(numel(domain)*6,400);
    bins(:,1) = 1;
    for M = 2:size(vessels,1)
        Intersections = cubeIntersect(vessels(M,3:5),vessels(vessels(M,7),3:5));
        for Loop = 1:size(Intersections,1)
            I = max([1,min([size(domain,1),Intersections(Loop,1)])]);
            J = max([1,min([size(domain,2),Intersections(Loop,2)])]);
            K = max([1,min([size(domain,3),Intersections(Loop,3)])]);
            if vessels(M, 2) ~= 0 && vessels(M, 2) ~= 10
                if vessels(M,2) < 10
                    idx_val = vessels(M,2);
                else
                    idx_val = vessels(M,2) - 10;
                end
                Family = reverseKey(idx_val);
                % bins{I,J,K,Family}(end+1) = M;
                % Update insert pos:
                bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),1) = bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),1) + 1;
                bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),1)) = M;
%                 if ~isempty(bins{I,J,K,Family})
%                     bins{I,J,K,Family} = sort(bins{I,J,K,Family});
%                     bins{I,J,K,Family} = bins{I,J,K,Family}([true;diff(bins{I,J,K,Family}(:))>0]);
%                 end
                 if bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),1) > 2
                    bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),2:bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),1)) = sort(bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),2:bins(sub2ind(size(domain),I,J,K)+(Family-1)*numel(domain),1)));
                    %bins{I,J,K,Family} = bins{I,J,K,Family}([true;diff(bins{I,J,K,Family}(:))>0]);
                end
            end
        end
    end
end

