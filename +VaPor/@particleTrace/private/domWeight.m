for DomainLoop = 1:obj.NumGM_WM
    DownTemp = [];
    
    [I,J,K] = ind2sub(size(flows.grey_white),obj.RowConvert(DomainLoop));
    
    Ax = flows.voxelSize^2; Ay = flows.voxelSize^2; Az = flows.voxelSize^2;
    if flows.porous.U(I,J,K)<0, DownTemp(end+1,:) = [obj.GM_WM_Convert(I-1,J,K),Ax*flows.bloodDensity*abs(flows.porous.U(I,J,K)),-obj.MTT(I,J,K)]; end % MTT for domain already calculated in massFlowrate
    if flows.porous.U(I+1,J,K)>0, DownTemp(end+1,:) = [obj.GM_WM_Convert(I+1,J,K),Ax*flows.bloodDensity*abs(flows.porous.U(I+1,J,K)),-obj.MTT(I,J,K)]; end
    if flows.porous.V(I,J,K)<0, DownTemp(end+1,:) = [obj.GM_WM_Convert(I,J-1,K),Ay*flows.bloodDensity*abs(flows.porous.V(I,J,K)),-obj.MTT(I,J,K)]; end
    if flows.porous.V(I,J+1,K)>0, DownTemp(end+1,:) = [obj.GM_WM_Convert(I,J+1,K),Ay*flows.bloodDensity*abs(flows.porous.V(I,J+1,K)),-obj.MTT(I,J,K)]; end
    if flows.porous.W(I,J,K)<0, DownTemp(end+1,:) = [obj.GM_WM_Convert(I,J,K-1),Az*flows.bloodDensity*abs(flows.porous.W(I,J,K)),-obj.MTT(I,J,K)]; end
    if flows.porous.W(I,J,K+1)>0, DownTemp(end+1,:) = [obj.GM_WM_Convert(I,J,K+1),Az*flows.bloodDensity*abs(flows.porous.W(I,J,K+1)),-obj.MTT(I,J,K)]; end
    
    if ~isempty(flows.Volume2VesselVein{I,J,K})
        for N = 1:size(flows.Volume2VesselVein{I,J,K},1)
            DownNode = flows.Volume2VesselVein{I,J,K}(N,1);
            if flows.SplitVein(DownNode,1)~=0
                DownTemp(end+1,:) = [DownNode+obj.NumArt+obj.NumGM_WM,flows.SplitVein(N,2)*abs(flows.Volume2VesselVein{I,J,K}(N,2)),-obj.MTT(I,J,K)];
                DownTemp(end+1,:) = [flows.veins.tree(DownNode,7)+NumArt+NumGM_WM,(1-flows.SplitVein(N,2))*abs(flows.Volume2VesselVein{I,J,K}(N,2)),-obj.MTT(I,J,K)];
            else
                if flows.FdotVein(DownNode,2) <= 0
                    DownNode = flows.veins.tree(DownNode,7);
                end
                DownTemp(end+1,:) = [DownNode+obj.NumArt+obj.NumGM_WM,abs(flows.Volume2VesselVein{I,J,K}(N,2)),-obj.MTT(I,J,K)];
                
            end
        end
    end
    
    obj.Tree{DomainLoop} = DownTemp;
    
end