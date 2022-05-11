for ArtLoop = 1:obj.NumArt
    
    if flows.SplitArt(ArtLoop,1)~=0
        for N = 1:size(flows.Vessel2VolumeArt{ArtLoop},1)
            if flows.Vessel2VolumeArt{ArtLoop}(N,4) ~= 0
                obj.Tree{ArtLoop+obj.NumGM_WM}(end+1,:) = [obj.GM_WM_Convert(flows.Vessel2VolumeArt{ArtLoop}(N,1),flows.Vessel2VolumeArt{ArtLoop}(N,2),flows.Vessel2VolumeArt{ArtLoop}(N,3))...
                    ,SplitArt(ArtLoop,2)*flows.Vessel2VolumeArt{ArtLoop}(N,4),0];
                obj.Tree{flows.arteries.tree(ArtLoop,7)+obj.NumGM_WM}(end+1,:) = [obj.GM_WM_Convert(flows.Vessel2VolumeArt{ArtLoop}(N,1),flows.Vessel2VolumeArt{ArtLoop}(N,2),flows.Vessel2VolumeArt{ArtLoop}(N,3))...
                    ,(1-SplitArt(ArtLoop,2))*flows.Vessel2VolumeArt{ArtLoop}(N,4),0];
            end
        end
        
        
    else
        
        
        if flows.FdotArt(ArtLoop,1) <= 0
            % probablility of travelling down artery length
            if abs(flows.FdotArt(ArtLoop,1)) > 0
                obj.Tree{ArtLoop+obj.NumGM_WM}(end+1,:) = [flows.arteries.tree(ArtLoop,7)+obj.NumGM_WM,abs(flows.FdotArt(ArtLoop,1)),-flows.geometry.art.V(ArtLoop)/(abs(flows.FdotArt(ArtLoop,1))/flows.bloodDensity)];
            end
            
            % probablility of transferring domain on artery length
            if ~isempty(flows.Vessel2VolumeArt{ArtLoop})
                for N = 1:size(flows.Vessel2VolumeArt{ArtLoop},1)
                    if flows.Vessel2VolumeArt{ArtLoop}(N,4) ~= 0
                        obj.Tree{ArtLoop+obj.NumGM_WM}(end+1,:) = [obj.GM_WM_Convert(flows.Vessel2VolumeArt{ArtLoop}(N,1),flows.Vessel2VolumeArt{ArtLoop}(N,2),flows.Vessel2VolumeArt{ArtLoop}(N,3))...
                            ,flows.Vessel2VolumeArt{ArtLoop}(N,4),0];
                    end
                end
            end
            
        elseif flows.FdotArt(ArtLoop,2) >= 0
            
            % probablility of travelling down artery length
            if abs(flows.FdotArt(ArtLoop,2)) > 0
                obj.Tree{flows.arteries.tree(ArtLoop,7)+obj.NumGM_WM}(end+1,:) = [ArtLoop+obj.NumGM_WM,abs(flows.FdotArt(ArtLoop,2)),flows.geometry.art.V(ArtLoop)/(abs(flows.FdotArt(ArtLoop,2))/flows.bloodDensity)];
            end
            
            % probablility of transferring domain on artery length
            if ~isempty(flows.Vessel2VolumeArt{ArtLoop})
                for N = 1:size(flows.Vessel2VolumeArt{ArtLoop},1)
                    if flows.Vessel2VolumeArt{ArtLoop}(N,4) ~= 0
                        obj.Tree{flows.arteries.tree(ArtLoop,7)+obj.NumGM_WM}(end+1,:) = [obj.GM_WM_Convert(flows.Vessel2VolumeArt{ArtLoop}(N,1),flows.Vessel2VolumeArt{ArtLoop}(N,2),flows.Vessel2VolumeArt{ArtLoop}(N,3))...
                            ,flows.Vessel2VolumeArt{ArtLoop}(N,4),0];
                    end
                end
            end
            
        end
        
    end
    
end

