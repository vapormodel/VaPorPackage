MassFlow = zeros(size(obj.grey_white));
MassFlow1 = zeros(size(obj.grey_white));
MassFlow2 = zeros(size(obj.grey_white));
MassFlowCheck = zeros(size(obj.grey_white));
U = obj.porous.U;
V = obj.porous.V;
W = obj.porous.W;
for I=1:size(obj.grey_white,1)
    for J=1:size(obj.grey_white,2)
        for K=1:size(obj.grey_white,3)
            if obj.grey_white(I,J,K)
                    MassFlow1(I,J,K) = MassFlow1(I,J,K) + abs(obj.MdotVoxelsArt(I,J,K));

%                     if U(I,J,K)>0, MassFlow1(I,J,K) = MassFlow1(I,J,K) + abs(U(I,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                     if U(I+1,J,K)<0, MassFlow1(I,J,K) = MassFlow1(I,J,K) + abs(U(I+1,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                     if V(I,J,K)>0, MassFlow1(I,J,K) = MassFlow1(I,J,K) + abs(V(I,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                     if V(I,J+1,K)<0, MassFlow1(I,J,K) = MassFlow1(I,J,K) + abs(V(I,J+1,K))*obj.voxelSize^2*obj.bloodDensity; end
%                     if W(I,J,K)>0, MassFlow1(I,J,K) = MassFlow1(I,J,K) + abs(W(I,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                     if W(I,J,K+1)<0, MassFlow1(I,J,K) = MassFlow1(I,J,K) + abs(W(I,J,K+1))*obj.voxelSize^2*obj.bloodDensity; end
%                     
%                     if obj.options.counterCurrentFlow
%                         MassFlow2(I,J,K) = MassFlow2(I,J,K) + (Perfusion(I,J,K)*obj.voxelSize^3);
%                         
%                         if U2(I,J,K)>0, MassFlow2(I,J,K) = MassFlow2(I,J,K) + abs(U2(I,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                         if U2(I+1,J,K)<0, MassFlow2(I,J,K) = MassFlow2(I,J,K) + abs(U2(I+1,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                         if V2(I,J,K)>0, MassFlow2(I,J,K) = MassFlow2(I,J,K) + abs(V2(I,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                         if V2(I,J+1,K)<0, MassFlow2(I,J,K) = MassFlow2(I,J,K) + abs(V2(I,J+1,K))*obj.voxelSize^2*obj.bloodDensity; end
%                         if W2(I,J,K)>0, MassFlow2(I,J,K) = MassFlow2(I,J,K) + abs(W2(I,J,K))*obj.voxelSize^2*obj.bloodDensity; end
%                         if W2(I,J,K+1)<0, MassFlow2(I,J,K) = MassFlow2(I,J,K) + abs(W2(I,J,K+1))*obj.voxelSize^2*obj.bloodDensity; end
%                         
%                          MassFlowCheck(I,J,K) = obj.voxelSize^2*obj.bloodDensity*(U(I,J,K)-U(I+1,J,K)+V(I,J,K)-V(I,J+1,K)+W(I,J,K)-W(I,J,K+1))+...
%                              obj.voxelSize^2*obj.bloodDensity*(U2(I,J,K)-U2(I+1,J,K)+V2(I,J,K)-V2(I,J+1,K)+W2(I,J,K)-W2(I,J,K+1))+...
%                              Mdot1(I,J,K) + Mdot2(I,J,K);
%                     else
%                          MassFlowCheck(I,J,K) = obj.voxelSize^2*obj.bloodDensity*(U(I,J,K)-U(I+1,J,K)+V(I,J,K)-V(I,J+1,K)+W(I,J,K)-W(I,J,K+1)) + obj.MdotVoxelsArt(I,J,K) + obj.MdotVoxelsVein(I,J,K);
%                     end
                       
                
            end
        end
    end
end

MassFlow = MassFlow1+MassFlow2;

MassFlow(~obj.grey_white) = NaN; % kg/s
MassFlow1(~obj.grey_white) = NaN; % kg/s
MassFlow2(~obj.grey_white) = NaN; % kg/s
MeasuredPerfusion = MassFlow/obj.voxelSize^3./obj.density/obj.bloodDensity*10^6/10*60; % ml/100g/min

