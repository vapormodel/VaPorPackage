classdef particleTrace < handle
    %PARTICLETRACE Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        MTT
        Tree
        T_TransientStore
        T_TransientStore_Save
        new_inlet
    end
    properties (Access = private)
        NumArt
        NumVein
        NumGM_WM
        NumDomTot
        RowConvert
        GM_WM_Convert
        flowobj
    end
    
    methods
        function obj = particleTrace(flows)
            %PARTICLETRACE Construct an instance of this class
            %   Detailed explanation goes here
            
            obj.flowobj = flows;
            
            obj.NumArt = size(flows.arteries.tree,1);
            obj.NumVein = size(flows.veins.tree,1);
            obj.NumGM_WM = numel(flows.grey_white(flows.grey_white));
            obj.NumDomTot = obj.NumGM_WM;
            obj.RowConvert = find(flows.grey_white);  % goes from row number to I,J,K
            obj.GM_WM_Convert = zeros(size(flows.grey_white));
            obj.GM_WM_Convert(obj.RowConvert) =  (1:numel(flows.grey_white(flows.grey_white)));
            obj.MTT = (flows.porosity*flows.voxelSize^3)./(flows.massFlow/flows.bloodDensity); % s
            
            obj.Tree = cell(obj.NumGM_WM+obj.NumArt+obj.NumVein,1);
            
            artWeight;
            domWeight;
            veinWeight;
            
            for N = 1:size(obj.Tree,1)
                if ~isempty(obj.Tree{N})
                    obj.Tree{N}(:,2) = obj.Tree{N}(:,2)/sum(obj.Tree{N}(:,2));
                end
            end
        end
        
        obj = traceParticles(obj, numParticles, parallel, overrideInlet);

    end
    
    methods(Static)
        exportPy(flows, filename)
    end
end

