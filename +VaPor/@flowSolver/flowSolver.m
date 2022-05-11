classdef flowSolver < handle
    %FLOWSOLVER Summary of this class goes here
    %   Detailed explanation goes here
    
    properties
        arteries
        veins
        tissue
        grey_white
        bloodFlow
        voxelSize
        porosity
        tortuosity
        capDia
        density
        
        bloodDensity
        bloodVisc
        
        MdotVesselArt
        Vessel2VolumeArt
        Volume2VesselArt
        MdotVoxelsArt
        
        MdotVesselVein
        Vessel2VolumeVein
        Volume2VesselVein
        MdotVoxelsVein
        
        MdotVoxelsOverall
        MdotVoxelsOverallB
        
        FdotArt
        FdotVein
        SplitArt
        SplitVein
        
        strokeLocation
        strokeTerminations
        strokeMask
        deductions
        
        options
        
        measuredPerfusion
        preStrokePerfusion
        porous
        
        geometry
        
        massFlow
    end
    
    methods
        function obj = flowSolver(arteries, veins, domain, gm_wm, bloodFlow, voxelSize, porosity, tortuosity, rho_b, bloodVisc, capDia, density, strokeLocation)
            %FLOWSOLVER Construct an instance of this class
            %   Detailed explanation goes here
            obj.arteries = arteries;
            obj.veins = veins;
            obj.tissue = domain;
            obj.grey_white = gm_wm;
            obj.bloodFlow = bloodFlow;
            obj.voxelSize = voxelSize;
            obj.porosity = porosity;
            obj.tortuosity = tortuosity;
            obj.bloodDensity = rho_b;
            obj.bloodVisc = bloodVisc;
            obj.capDia = capDia;
            obj.density = density;
            obj.strokeLocation = strokeLocation;
            
            obj.options.counterCurrentFlow = false;
            
            obj.vesselDomainExchange;
        end
        
        function perfPlot(obj)
            %METHOD1 Summary of this method goes here
            %   Detailed explanation goes here
            planecut('z', 30, obj.measuredPerfusion);
        end
        
        newStrokeAdjustmentSimple(obj, removed_flow)
        
    end
end

