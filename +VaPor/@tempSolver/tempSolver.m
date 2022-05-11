classdef tempSolver < handle
    %UNTITLED Summary of this class goes here
    %   Detailed explanation goes here
    
    % Public Properties
    properties
        % To Move Elsewhere
        bloodCp
        bloodKc
        
        Kc
        Cp
        Q
        H
        
        Nu
        Beta
        
        %T_Solve
        D
        temperatures
        RowAdjustment
        
        % DECOMP:
        %luFac
        
        Transient_Temp
        save_timestep
    end
    properties (Transient = true)
        T_Solve
        luFac
    end
    % Private Properties
    properties (Access = private)
        flow_obj
        solve_method
    end
    
    methods
        function obj = tempSolver(flows, blood_cp, blood_kc, Nu, tissue_q, tissue_cp, tissue_kc, h)
            %TEMPSOLVER Construct an instance of this class
            %   This is the constructor method - a copy of this class is
            %   instantiated here.
            obj.flow_obj = flows;
            obj.solve_method = 'bicg';
            
            %These will need to be reconfigured:
            obj.bloodCp = blood_cp;
            obj.bloodKc = blood_kc;
            obj.Nu.base = Nu;
            obj.Beta.porous = Inf;
            obj.Kc = tissue_kc;
            obj.Cp = tissue_cp;
            obj.Q = tissue_q;
            obj.H = h;
%             if ~isempty(flows.preStrokePerfusion)
%                 reduction_matrix = 1-(flows.preStrokePerfusion-flows.measuredPerfusion)./flows.preStrokePerfusion;
%                 for i = 1:size(reduction_matrix, 1)
%                     for j = 1:size(reduction_matrix, 2)
%                         for k = 1:size(reduction_matrix, 3)
%                             if isnan(reduction_matrix(i,j,k)) && ~isnan(tissue_q(i,j,k))
%                                 reduction_matrix(i,j,k) = 1;
%                             end
%                         end
%                     end
%                 end
%                 obj.Q = obj.Q .* (0.2 + 0.8*reduction_matrix);
%             end
        end
        
        obj = solve(obj, T_amb, T_blood);
        
        obj = solveit(obj, T_amb, T_blood);
        
        plotbrain(obj, position);
        plot_trans_brain(obj, timestep, position,climits, full_head);
        
        function export_trans(obj)
            numberOfTimeSteps = size(obj.Transient_Temp.tissue,2);
            SaveTime = 1;
            trans3D = zeros(size(obj.flow_obj.tissue,1),size(obj.flow_obj.tissue,2),size(obj.flow_obj.tissue,3), numberOfTimeSteps);
            NumDomTot = numel(obj.flow_obj.tissue(obj.flow_obj.tissue));
            RowConvert = find(obj.flow_obj.tissue);
            for marker = 1:numberOfTimeSteps
                % Assemble the space matrix:
                thisData = obj.Transient_Temp.tissue(:,marker);
                Tt = zeros(size(obj.flow_obj.tissue));
                Tt(RowConvert) = thisData(1:NumDomTot);
                Tt(~obj.flow_obj.tissue) = NaN;
                Tt(~obj.flow_obj.grey_white) = NaN;
                
                trans3D(:,:,:,marker) = Tt;
            end
            
            nccreate('timeexport.nc','Transient', 'Format', 'netcdf4','DeflateLevel', 2,  'Dimensions',{'x', size(trans3D,1), 'y',size(trans3D,2),'z',size(trans3D,3),'time',size(trans3D,4)});
            ncwrite('timeexport.nc', 'Transient', trans3D);
            
            nccreate('timeexport.nc', 'time', 'Format', 'netcdf4', 'DeflateLevel', 2,  'Dimensions', {'time', numberOfTimeSteps});
            ncwrite('timeexport.nc', 'time', linspace(0,numberOfTimeSteps*SaveTime,numberOfTimeSteps));
            ncwriteatt('timeexport.nc','time','units','seconds since 2010-11-9 00:00:00 +0:00')
        end
        
        function [suggested_lim] = bounds_extractor(obj)
            % Locate the rows which correspond to the brain (rather than
            % all tissue)
            RowConvert = find(obj.flow_obj.tissue);
            running_mat = [];
            for row_idx = 1:length(RowConvert)
                if obj.flow_obj.grey_white(RowConvert(row_idx))
                    running_mat = [running_mat, obj.Transient_Temp.tissue(row_idx,:)];
                end
            end
            upper_lim = max(running_mat, [], 'all');
            lower_lim = min(running_mat, [], 'all');
            suggested_lim = [lower_lim, upper_lim];
        end
        
    end
    
    methods (Access = protected)
        obj = matrix_builder(obj);
    end
end

