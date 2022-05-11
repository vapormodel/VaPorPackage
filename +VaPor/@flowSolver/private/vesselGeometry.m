function [L, AL, Vol, Davg, Aavg] = vesselGeometry(Vessel,h)

% Vessel(:,3:5) = Vessel(:,3:5) - 0.5; % Vessels need to be adjusted by half a node because they are defined geometrically by nodes


% Preallocate data vectors
D=zeros(size(Vessel,1),1);
Davg=zeros(size(Vessel,1),1);
L=zeros(size(Vessel,1),1);
AL=zeros(size(Vessel,1),1);
Vol=zeros(size(Vessel,1),1);
Aavg=zeros(size(Vessel,1),1);

% First Node
D(1) = Vessel(1,6);
Davg(1) = Vessel(1,6);

% Every Other Node
for n=2:size(Vessel,1)
    D(n) = Vessel(n,6);
    
    curr = [Vessel(n,3),Vessel(n,4),Vessel(n,5)]; %  Endpoint on Vessel
    if curr(1) ~= 0
        prev = [Vessel(Vessel(n,7),3),Vessel(Vessel(n,7),4),Vessel(Vessel(n,7),5)];
        L(n) = h*norm(curr-prev); % Length of Vessel Segment
        AL(n) = sqrt(((D(n) - Vessel(Vessel(n,7),6))/2)^2+L(n)^2)*pi*(D(n) + Vessel(Vessel(n,7),6))/2;
        Davg(n)=(D(n)+Vessel(Vessel(n,7),6))/2;
        Aavg(n) = (1/12)*pi*(D(n)^2 + Vessel(Vessel(n,7),6)^2 + D(n)*Vessel(Vessel(n,7),6));
        Vol(n) = L(n)*Aavg(n);
    end
end

end