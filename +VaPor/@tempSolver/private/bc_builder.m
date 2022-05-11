function bc_builder(obj, T_Out, BloodTemp, InletTemp)
%BC_BUILDER Builds boundary condition vector (D).

H_Out = obj.H;

%%%%%% Setting Up Reference Lists for obj.flow_obj.tissue %%%%
NumDomTot = numel(obj.flow_obj.tissue(obj.flow_obj.tissue));
RowConvert = find(obj.flow_obj.tissue);  % goes from row number to I,J,K
DomTotConvert = zeros(size(obj.flow_obj.tissue));
DomTotConvert(RowConvert) = 1:NumDomTot;  % goes from I,J,K to row number
% The sparse matrix to solve only needs those voxels that are within
% obj.flow_obj.tissue. These lists help convert sparse matrix row number to location in
% obj.flow_obj.tissue (given by I,J,K) and from locations in obj.flow_obj.tissue to row numbers in
% the sparse matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

Borders = bwperim(obj.flow_obj.tissue);


%%%%%% Check if 2 Phases in Porous Domain %%%%%%%
RowAdjustment = 0; % Value that adds rows into the matrix if required. (RowAdjustment = 0 adds no rows).
if obj.Beta.porous ~= Inf % If the inter-domain heat transfer is finite.
    RowAdjustment = NumDomTot; % Add rows to solve the tissue and blood phases separately.
end

obj.RowAdjustment = RowAdjustment;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Establishing Number of Rows in Domains %%%
NumDomRows = NumDomTot + RowAdjustment; % Total rows assigned to voxels.
VesselRow = size(obj.flow_obj.arteries.tree,1); % Total rows assigned to arteries.
% In the sparse matrix, the rows are divided up by domiain in the following
% order: [tissue -> blood -> arteries -> veins]. These values allow rapid
% identification of rows for arteries and veins.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%



D = zeros(NumDomRows+size(obj.flow_obj.arteries.tree,1)+size(obj.flow_obj.veins.tree,1),1);
% T_super will be a list of all interactions that will be fed into the
% sparse matrix stored as [row, col, value]. To create the list, it is fed
% by T_DataTemp on every iteration through the domain which is filled up
% by interactions. T_DataRowcount keeps track of where to add T_DataTemp
% at the end of T_Super.
% D contains all the values that are not variable dependent (in this case
% all the Mdot values). This will be the RHS of the equation in the linear
% solver (Mx = D where M is the sparse matrix, x is the list of variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling 3D Voxel Interactions %%%%%%%%%%
disp('Compiling Interactions within 3D Voxels') % Display.
tic % Start timing.

for I=1:size(obj.flow_obj.tissue,1)
    for J=1:size(obj.flow_obj.tissue,2)
        for K=1:size(obj.flow_obj.tissue,3)
            if obj.flow_obj.tissue(I,J,K)
                % For all voxels, if they are within obj.flow_obj.tissue.
                
                
                %%%%%% Establishing Row Values %%%%%%%%%%%%%%%%%%
                Row = DomTotConvert(I,J,K); BRow = Row + RowAdjustment;
                if I<size(obj.flow_obj.grey_white,1), RowUI = DomTotConvert(I+1,J,K); BRowUI = RowUI+RowAdjustment; end
                if I>1,             RowDI = DomTotConvert(I-1,J,K); BRowDI = RowDI+RowAdjustment; end
                if J<size(obj.flow_obj.grey_white,2), RowUJ = DomTotConvert(I,J+1,K); BRowUJ = RowUJ+RowAdjustment; end
                if J>1,             RowDJ = DomTotConvert(I,J-1,K); BRowDJ = RowDJ+RowAdjustment; end
                if K<size(obj.flow_obj.grey_white,3), RowUK = DomTotConvert(I,J,K+1); BRowUK = RowUK+RowAdjustment; end
                if K>1,             RowDK = DomTotConvert(I,J,K-1); BRowDK = RowDK+RowAdjustment; end
                % Establishes the row values in the matrix for the current
                % voxel and all neighbouring voxels. If there is only a
                % single phase being solved (because Beta23 = Inf) then all
                % BRow values equal Row values.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Preallocating Searchable Parameters %%%%%%
                Vol=obj.flow_obj.voxelSize^3;
                Ax=obj.flow_obj.voxelSize^2;
                Ay=obj.flow_obj.voxelSize^2;
                Az=obj.flow_obj.voxelSize^2;
                % Volume and areas of voxels.
                Option_CounterCurrentFlow = false;
                if Option_CounterCurrentFlow
                    Porosity_b = 2*Porosity(I,J,K);
                    Porosity_t = 1-2*Porosity_b;
                else
                    Porosity_b = obj.flow_obj.porosity(I,J,K);
                    Porosity_t = 1-Porosity_b;
                end
                
                Rho_t = obj.flow_obj.density(I,J,K);
                Cp_t = obj.Cp(I,J,K);
                Kc_t = obj.Kc(I,J,K);
                Q_t = obj.Q(I,J,K);
                % Establishing searchable parameters. This prevents
                % constantly having to look up values in matrix.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Boundary Interactions %%%%%%%%%%%%%%%%%%%%
                if Borders(I,J,K) % If voxel lies on domain boundary.
                    BorderMinX=0; BorderMaxX=0; BorderMinY=0; BorderMaxY=0; BorderMinZ=0; BorderMaxZ=0;
                    % Sets all boundary heat transfer to zero
                    
                    Option_LimitedCooling = false;
                    %%%%%% Check for Boundary Heat Transfer %%%%%
                    if ~Option_LimitedCooling || K >= LimitedCoolingHeight % Check if Option_LimitedCooling is disabled or both enabled and conditions are met.
                        if I==1              || ~obj.flow_obj.tissue(I-1,J,K), BorderMinX = 1; end
                        if I==size(obj.flow_obj.tissue,1) || ~obj.flow_obj.tissue(I+1,J,K), BorderMaxX = 1; end
                        if J==1              || ~obj.flow_obj.tissue(I,J-1,K), BorderMinY = 1; end
                        if J==size(obj.flow_obj.tissue,2) || ~obj.flow_obj.tissue(I,J+1,K), BorderMaxY = 1; end
                        if K==1              || ~obj.flow_obj.tissue(I,J,K-1), BorderMinZ = 1; end
                        if K==size(obj.flow_obj.tissue,3) || ~obj.flow_obj.tissue(I,J,K+1), BorderMaxZ = 1; end
                        % Check which faces are on bounderies for heat
                        % transfer to occur.
                    end
                    Option_AdiabaticBase = true;
                    if Option_AdiabaticBase && K == 1, BorderMinZ=0; end
                    % Set base of model to adiabatic.
                    
                    %if Option_AdiabaticArmEnds && K == find(any(any(GM_WM))==1,1,'first'), BorderMinZ=0; end
                    %if Option_AdiabaticArmEnds && K == find(any(any(GM_WM))==1,1,'last'), BorderMaxZ=0; end
                    % Set ends of arm model to adiabatic.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Establish Boundary Interactions %%%%%%%%%%
                    if H_Out == Inf % If boundary heat transfer is infinite, set voxel temperature equal to boundary temperature.
                        D(Row) = D(Row) + -(BorderMinX+BorderMaxX)*T_Out; % Tissue domain boundaries in the X direction.
                        D(Row) = D(Row) + -(BorderMinY+BorderMaxY)*T_Out; % Tissue domain boundaries in the Y direction.
                        D(Row) = D(Row) + -(BorderMinZ+BorderMaxZ)*T_Out; % Tissue domain boundaries in the Z direction.
                        D(BRow) = D(BRow) + -(BorderMinX+BorderMaxX)*T_Out; % Blood domain boundaries in the X direction.
                        D(BRow) = D(BRow) + -(BorderMinY+BorderMaxY)*T_Out; % Blood domain boundaries in the Y  direction.
                        D(BRow) = D(BRow) + -(BorderMinZ+BorderMaxZ)*T_Out; % Blood domain boundaries in the Z direction.
                        % Establish all interactions with boundaries if
                        % required.
                        
                    else % If boundary heat transfer is finite do heat transfer interactions between voxels and boundaries.
                        D(Row) = D(Row) + Porosity_t*(BorderMinX+BorderMaxX)*Ax*(-H_Out*T_Out); % Tissue domain boundaries in the X direction.
                        D(Row) = D(Row) + Porosity_t*(BorderMinY+BorderMaxY)*Ay*(-H_Out*T_Out); % Tissue domain boundaries in the Y direction.
                        D(Row) = D(Row) + Porosity_t*(BorderMinZ+BorderMaxZ)*Az*(-H_Out*T_Out); % Tissue domain boundaries in the Z direction.
                        D(BRow) = D(BRow) + Porosity_b*(BorderMinX+BorderMaxX)*Ax*(-H_Out*T_Out); % Blood domain boundaries in the X direction.
                        D(BRow) = D(BRow) + Porosity_b*(BorderMinY+BorderMaxY)*Ay*(-H_Out*T_Out); % Blood domain boundaries in the Y  direction.
                        D(BRow) = D(BRow) + Porosity_b*(BorderMinZ+BorderMaxZ)*Az*(-H_Out*T_Out); % Blood domain boundaries in the Z direction.
                        % Establish all interactions with boundaries if
                        % required.
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    
                    
                end
                
                %%%%%% Establish Metabolism %%%%%%%%%%%%%%%%%%%%%
                D(Row) = D(Row) - Vol*Q_t; % Tissue domain source for source.
                % Establishes the metabolic heat source term.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Establish Pennes Perfusion %%%%%%%%%%%%%%%
                if ~obj.flow_obj.grey_white(I,J,K)  % If outside of the brain domain.
                   D(Row) = D(Row) - Vol*obj.flow_obj.arteries.domainWeighting(I,J,K)*obj.bloodCp*BloodTemp; % Tissue domain interaction for source.
                end
                % Establishes the Pennes Perfusion source term for tissue
                % outside of the brain domain.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

            end
        end
    end
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Artery Segment Interactions %%%%
disp('Compiling Interactions within Arteries')
tic % Start timing.

for Node = 1:size(obj.flow_obj.Vessel2VolumeArt,1)
    
    %%%%%% Inlet Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(Node,obj.flow_obj.arteries.termPoints)
        Index = find(Node==obj.flow_obj.arteries.termPoints); % Find which inlet is being referred to.
        D(NumDomRows+Node) = InletTemp(Index); % Set node equal to InletTemp.
        % If node is an inlet, set temperature to corresponding InletTemp.
        
    end
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
obj.D = D;

end

