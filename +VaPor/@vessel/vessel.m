classdef vessel < handle
    %VESSEL Holds vessel objects.
    %   This file is part of Craft. This is the vessel object used in
    %   VaPor.
    %   Author:         Luke Fulford - luke.fulford@ed.ac.uk
    %   Last Modified:  03/02/2020
    
    properties
        tree
        termPoints
        termFlows
        generateVessels
        generateIterations
        weightFactor
        useRRTstar
        RRTstarEpsilon
        allowZero
        joinRandom
        confineJoin
        regionGen
%     end
%     properties (Access = protected)
        basePath
        cleanFunc
        diameterBackup
        termPointBackup % Not sure this will be used.
        domainLimits
        domainWeighting
    end
    
    methods
        function obj = vessel(propStruct)
            %VESSEL Construct an instance of this class
            %  This is the constructor. Here the base tree and options
            %  struct are loaded and parsed.
            
            % Unpack the structure:
            obj.basePath = propStruct.vesselPath;
            obj.cleanFunc = propStruct.vesselCleanFunc;
            obj.termPoints = propStruct.vesselSpecTermPoints;
            obj.termFlows = propStruct.vesselSpecTermFlows;
            obj.generateVessels = propStruct.generateVessels;
            obj.generateIterations = propStruct.generateIterations;
            obj.weightFactor = propStruct.weightFactor;
            obj.useRRTstar = propStruct.useRRTstar;
            obj.RRTstarEpsilon = propStruct.RRTstarEpsilon;
            obj.domainLimits = propStruct.domainLimits;
            obj.domainWeighting = propStruct.domainWeighting;
            obj.allowZero = propStruct.allowZero;
            obj.joinRandom = propStruct.joinRandom;
            obj.confineJoin = propStruct.confineJoin;
            obj.regionGen = propStruct.regionGen;
            
            % From the base vessel tree path load the vessel:
            obj.loadVesselSWC(obj.basePath, obj.cleanFunc);
            
            if isempty(obj.termPoints)
                obj.termPoints = auto_inlet(obj.tree);
            end
            
            % Vessel Diameters are changed as part of the flowrate solving,
            % store a copy of the originals here, just in case we want to
            % check how they change later.
            obj.diameterBackup = obj.tree(:,6);
        end
        
        function draw(obj)
            %DRAW Draws the vessel tree.
            %   This method calls the drawVessels helper function to draw
            %   the vessel tree. As a conveniece it will also check that
            %   the number of vessels is not too excessive.
            
            % Check the number of vessels if this is greater than 15,000
            % warn the user & ask for confirmation.
            
            if length(obj.tree) > 10000
                % The tree is large, ask for confirmation:
                questString =...
                    sprintf(['Warning! This tree has %d nodes, this may'...
                    ' take a long time to draw! Do you want to continue'...
                    '? [y/N]: '], length(obj.tree));
                res = input(questString, 's');
                if isempty(res)
                    res = 'N';
                end
                switch res
                    case {'y', 'Y'}
                        continueDraw = 1;
                    otherwise
                        continueDraw = 0;
                end
            else
                continueDraw = 1;
            end
            if continueDraw
                drawVessels(obj.tree);
            end
        end
        
        function export(obj)
            %EXPORT exports the generated tree as a .swc file.
            %tempTree = obj.tree;
            %save('export.swc', 'tempTree', '-ascii', '-double');
            dlmwrite('export.swc',obj.tree,'precision','%.6f');
        end
        function save(obj, type)
            %SAVE saves the vessel class instance as a .mat file.
            switch type
                case 'a'
                    arteries = obj;
                    save('test.mat', 'arteries');
                case 'v'
                    veins = obj;
                    save('test.mat', 'veins');
                otherwise
                    save('test.mat', 'obj');
            end
        end
        function craft(obj)
            %CRAFT Expands vessel trees.
            %   Currently this function suports the normal cases of both
            %   RRT and RRT*.
            
            % Backup the terminal points:
            obj.termPointBackup = obj.tree(obj.termPoints, :);
            %             obj.tree = genVessels(obj.tree, obj.domainWeighting,...
            %                 obj.domainLimits, obj.generateIterations,...
            %                 obj.weightFactor, obj.allowZero, obj.joinRandom,...
            %                 obj.RRTstarEpsilon);

                if obj.useRRTstar
                    % Expand the tree using RRT*
                    obj.tree = obj.genVessels(obj.tree, obj.domainWeighting,...
                        obj.domainLimits, obj.generateIterations,...
                        obj.weightFactor, obj.allowZero, obj.joinRandom,...
                        obj.confineJoin, obj.regionGen, obj.RRTstarEpsilon);
                else
                    obj.tree = obj.genVessels(obj.tree, obj.domainWeighting,...
                        obj.domainLimits, obj.generateIterations,...
                        obj.weightFactor, obj.allowZero, obj.joinRandom,...
                        obj.confineJoin, obj.regionGen);
                end

        end
        function craftOld(obj)
            %CRAFT Expands vessel trees.
            %   Currently this function suports the normal cases of both
            %   RRT and RRT*.
            
            % Backup the terminal points:
            obj.termPointBackup = obj.tree(obj.termPoints, :);
            %obj.tree = genVessels(obj.tree, obj.domainWeighting, obj.domainLimits, obj.generateIterations, obj.weightFactor, obj.RRTstarEpsilon);
            if obj.useRRTstar
                % Expand the tree using RRT*
                obj.tree = vesselGenerationRRT_star(obj.tree,...
                    obj.domainWeighting,...
                    obj.domainLimits,...
                    obj.generateIterations,...
                    obj.weightFactor,...
                    obj.RRTstarEpsilon);
            else
                obj.tree = vesselGenerationRRT(obj.tree,...
                    obj.domainWeighting,...
                    obj.domainLimits,...
                    obj.generateIterations,...
                    obj.weightFactor);
            end
        end
        function growTree(obj, numIterations)
            obj.tree = mallocTree(obj.tree, numIterations);
        end
        
        function craftHandle(obj)
            genVesselsHandle(obj)
            
        end
        
        function add(obj, newRow)
            %ADD This method inserts a row into the tree - it can also
            %modify a row.
            
            
            
        end
    end
    methods(Static)
        [newTree,probCumSum] = genVessels(existingTree,probMap,domain,...
            numIterations,weightFactor,allowZero,joinRandom,confineJoin,regionGen,RRTstarLimit);
    end
end

