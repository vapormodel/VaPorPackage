for VeinLoop = 1:obj.NumVein
    
    if flows.FdotVein(VeinLoop,2) < 0 % && abs(FdotVein(VeinLoop,2)) > ErrorAmt
        obj.Tree{VeinLoop+obj.NumArt+obj.NumGM_WM}(end+1,:) = [flows.veins.tree(VeinLoop,7)+obj.NumArt+obj.NumGM_WM,abs(flows.FdotVein(VeinLoop,2)),-flows.geometry.vein.V(VeinLoop)/(abs(flows.FdotVein(VeinLoop,2))/flows.bloodDensity)];
    elseif flows.FdotVein(VeinLoop,1) > 0 % && abs(FdotVein(VeinLoop,1)) > ErrorAmt
        obj.Tree{flows.veins.tree(VeinLoop,7)+obj.NumArt+obj.NumGM_WM}(end+1,:) = [VeinLoop+obj.NumArt+obj.NumGM_WM,abs(flows.FdotVein(VeinLoop,1)),flows.geometry.vein.V(VeinLoop)/(abs(flows.FdotVein(VeinLoop,1))/flows.bloodDensity)];
    end
    
end