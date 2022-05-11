function matrix_builder_CN(obj, Timestep)
%MATRIX_BUILDER builds the sparse matrix for temperature solving.
%   Detailed explanation goes here

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


%%%%%% Initialising Storage Matrices %%%%%%%%%%%%
T_DataTemp = zeros(3e8,3);
T_DataTempPos = 0;
T_DataRowcount = 1;
T_Super = zeros(3e8,3);
%D = zeros(NumDomRows+size(obj.flow_obj.arteries.tree,1)+size(obj.flow_obj.veins.tree,1),1);
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
                    
                    %if Option_AdiabaticArmEnds && K == find(any(any(obj.flow_obj.grey_white))==1,1,'first'), BorderMinZ=0; end
                    %if Option_AdiabaticArmEnds && K == find(any(any(obj.flow_obj.grey_white))==1,1,'last'), BorderMaxZ=0; end
                    % Set ends of arm model to adiabatic.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Establish Boundary Interactions %%%%%%%%%%
                    if H_Out == Inf % If boundary heat transfer is infinite, set voxel temperature equal to boundary temperature.
                        %D(Row) = D(Row) + -(BorderMinX+BorderMaxX)*T_Out; % Tissue domain boundaries in the X direction.
                        %D(Row) = D(Row) + -(BorderMinY+BorderMaxY)*T_Out; % Tissue domain boundaries in the Y direction.
                        %D(Row) = D(Row) + -(BorderMinZ+BorderMaxZ)*T_Out; % Tissue domain boundaries in the Z direction.
                        %D(BRow) = D(BRow) + -(BorderMinX+BorderMaxX)*T_Out; % Blood domain boundaries in the X direction.
                        %D(BRow) = D(BRow) + -(BorderMinY+BorderMaxY)*T_Out; % Blood domain boundaries in the Y  direction.
                        %D(BRow) = D(BRow) + -(BorderMinZ+BorderMaxZ)*T_Out; % Blood domain boundaries in the Z direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(BorderMinX+BorderMaxX)]; % Tissue domain boundaries in the X direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(BorderMinY+BorderMaxY)]; % Tissue domain boundaries in the Y direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(BorderMinZ+BorderMaxZ)]; % Tissue domain boundaries in the Z direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(BorderMinX+BorderMaxX)]; % Blood domain boundaries in the X direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(BorderMinY+BorderMaxY)]; % Blood domain boundaries in the Y direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(BorderMinZ+BorderMaxZ)]; % Blood domain boundaries in the Z direction.
                        % Establish all interactions with boundaries if
                        % required.
                        
                    else % If boundary heat transfer is finite do heat transfer interactions between voxels and boundaries.
                        %D(Row) = D(Row) + Porosity_t*(BorderMinX+BorderMaxX)*Ax*(-H_Out*T_Out); % Tissue domain boundaries in the X direction.
                        %D(Row) = D(Row) + Porosity_t*(BorderMinY+BorderMaxY)*Ay*(-H_Out*T_Out); % Tissue domain boundaries in the Y direction.
                        %D(Row) = D(Row) + Porosity_t*(BorderMinZ+BorderMaxZ)*Az*(-H_Out*T_Out); % Tissue domain boundaries in the Z direction.
                        %D(BRow) = D(BRow) + Porosity_b*(BorderMinX+BorderMaxX)*Ax*(-H_Out*T_Out); % Blood domain boundaries in the X direction.
                        %D(BRow) = D(BRow) + Porosity_b*(BorderMinY+BorderMaxY)*Ay*(-H_Out*T_Out); % Blood domain boundaries in the Y  direction.
                        %D(BRow) = D(BRow) + Porosity_b*(BorderMinZ+BorderMaxZ)*Az*(-H_Out*T_Out); % Blood domain boundaries in the Z direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,Porosity_t*(BorderMinX+BorderMaxX)*Ax*-H_Out]; % Tissue domain boundaries in the X direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,Porosity_t*(BorderMinY+BorderMaxY)*Ay*-H_Out]; % Tissue domain boundaries in the Y direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,Porosity_t*(BorderMinZ+BorderMaxZ)*Az*-H_Out]; % Tissue domain boundaries in the Z direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,Porosity_b*(BorderMinX+BorderMaxX)*Ax*-H_Out]; % Blood domain boundaries in the X direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,Porosity_b*(BorderMinY+BorderMaxY)*Ay*-H_Out]; % Blood domain boundaries in the Y direction.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,Porosity_b*(BorderMinZ+BorderMaxZ)*Az*-H_Out]; % Blood domain boundaries in the Z direction.
                        % Establish all interactions with boundaries if
                        % required.
                    end
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Establish Conductive Interactions %%%%%%%%
                    if ~(I==1 || ~obj.flow_obj.tissue(I-1,J,K)) % If no boundary between current voxel and I-1.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowDI,(Porosity_t*Ax*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(Porosity_t*Ax*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDI,(Porosity_b*Ax*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(Porosity_b*Ax*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with I-1 if not a boundary.
                    
                    if ~(I==size(obj.flow_obj.tissue,1) || ~obj.flow_obj.tissue(I+1,J,K)) % If no boundary between current voxel and I+1.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowUI,(Porosity_t*Ax*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(Porosity_t*Ax*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUI,(Porosity_b*Ax*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(Porosity_b*Ax*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with I+1 if not a boundary.
                    
                    if ~(J==1 || ~obj.flow_obj.tissue(I,J-1,K)) % If no boundary between current voxel and J-1.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowDJ,(Porosity_t*Ay*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(Porosity_t*Ay*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDJ,(Porosity_b*Ay*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(Porosity_b*Ay*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with J-1 if not a boundary.
                    
                    if ~(J==size(obj.flow_obj.tissue,2) || ~obj.flow_obj.tissue(I,J+1,K)) % If no boundary between current voxel and J+1.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowUJ,(Porosity_t*Ay*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(Porosity_t*Ay*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUJ,(Porosity_b*Ay*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(Porosity_b*Ay*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with J+1 if not a boundary.
                    
                    if ~(K==1 || ~obj.flow_obj.tissue(I,J,K-1)) % If no boundary between current voxel and K-1.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowDK,(Porosity_t*Az*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(Porosity_t*Az*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDK,(Porosity_b*Az*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(Porosity_b*Az*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with K-1 if not a boundary.
                    
                    if ~(K==size(obj.flow_obj.tissue,3) || ~obj.flow_obj.tissue(I,J,K+1)) % If no boundary between current voxel and K+1.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowUK,(Porosity_t*Az*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(Porosity_t*Az*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                        
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUK,(Porosity_b*Az*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(Porosity_b*Az*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    end
                    % Conduction with K+1 if not a boundary.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                    %%%%%% Non-Boundary Interactions %%%%%%%%%%%%%%%%
                else % If voxel does not lie on domain boundary.
                    
                    
                    %%%%%% Establish Conductive Interactions %%%%%%%%
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowUI,(Porosity_t*Ax*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel (I+1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowDI,(Porosity_t*Ax*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel (I-1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(2*Porosity_t*Ax*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUI,(Porosity_b*Ax*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel (I+1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDI,(Porosity_b*Ax*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel (I-1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(2*Porosity_b*Ax*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    % Conduction with I-1 and I+1.
                    
                   T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowUJ,(Porosity_t*Ay*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel (J+1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowDJ,(Porosity_t*Ay*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel (J-1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(2*Porosity_t*Ay*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUJ,(Porosity_b*Ay*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel (J+1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDJ,(Porosity_b*Ay*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel (J-1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(2*Porosity_b*Ay*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    % Conduction with J-1 and J+1.
                    
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowUK,(Porosity_t*Az*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel (K+1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,RowDK,(Porosity_t*Az*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for neighbouring voxel (K-1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-(2*Porosity_t*Az*(Kc_t/obj.flow_obj.voxelSize))]; % Tissue domain interaction for current voxel.
                   T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUK,(Porosity_b*Az*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel (K+1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDK,(Porosity_b*Az*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for neighbouring voxel (K-1).
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-(2*Porosity_b*Az*(obj.bloodKc/obj.flow_obj.voxelSize))]; % Blood domain interaction for current voxel.
                    % Conduction with K-1 and K+1.
                    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                    
                    
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% 3D Advection Interactions %%%%%%%%%%%%%%%%
                if obj.flow_obj.porous.U(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDI,Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.U(I,J,K))]; end % Blood domain interaction for upstream voxel.
                if obj.flow_obj.porous.U(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.U(I,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from I-1.
                
                if obj.flow_obj.porous.U(I+1,J,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUI,Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.U(I+1,J,K))]; end % Blood domain interaction for upstream voxel .
                if obj.flow_obj.porous.U(I+1,J,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.U(I+1,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from I+1.
                
                if obj.flow_obj.porous.V(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDJ,Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.V(I,J,K))]; end % Blood domain interaction for upstream voxel.
                if obj.flow_obj.porous.V(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.V(I,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from J-1.
                
                if obj.flow_obj.porous.V(I,J+1,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUJ,Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.V(I,J+1,K))]; end % Blood domain interaction for upstream voxel.
                if obj.flow_obj.porous.V(I,J+1,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.V(I,J+1,K))]; end  % Blood domain interaction for current voxel.
                % Advection from J+1.
                
                if obj.flow_obj.porous.W(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDK,Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.W(I,J,K))]; end % Blood domain interaction for upstream voxel.
                if obj.flow_obj.porous.W(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.W(I,J,K))]; end  % Blood domain interaction for current voxel.
                % Advection from K-1.
                
                if obj.flow_obj.porous.W(I,J,K+1)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUK,Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.W(I,J,K+1))]; end % Blood domain interaction for upstream voxel.
                if obj.flow_obj.porous.W(I,J,K+1)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(obj.flow_obj.porous.W(I,J,K+1))]; end  % Blood domain interaction for current voxel.
                % Advection from K+1.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                Option_CounterCurrentFlow = false;
                Option_PseudoCounterCurrentFlow = false;
                %%%%%% 3D Advection Interactions (C-C) %%%%%%%%%%
                if Option_CounterCurrentFlow || Option_PseudoCounterCurrentFlow
                    if U2(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDI,Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(U2(I,J,K))]; end % Blood domain interaction for upstream voxel.
                    if U2(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(U2(I,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from I-1.
                    
                    if U2(I+1,J,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUI,Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(U2(I+1,J,K))]; end % Blood domain interaction for upstream voxel .
                    if U2(I+1,J,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ax*obj.flow_obj.bloodDensity*obj.bloodCp*abs(U2(I+1,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from I+1.
                    
                    if V2(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDJ,Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(V2(I,J,K))]; end % Blood domain interaction for upstream voxel.
                    if V2(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(V2(I,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from J-1.
                    
                    if V2(I,J+1,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUJ,Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(V2(I,J+1,K))]; end % Blood domain interaction for upstream voxel.
                    if V2(I,J+1,K)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Ay*obj.flow_obj.bloodDensity*obj.bloodCp*abs(V2(I,J+1,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from J+1.
                    
                    if W2(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowDK,Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(W2(I,J,K))]; end % Blood domain interaction for upstream voxel.
                    if W2(I,J,K)>0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(W2(I,J,K))]; end  % Blood domain interaction for current voxel.
                    % Advection from K-1.
                    
                    if W2(I,J,K+1)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRowUK,Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(W2(I,J,K+1))]; end % Blood domain interaction for upstream voxel.
                    if W2(I,J,K+1)<0, T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Az*obj.flow_obj.bloodDensity*obj.bloodCp*abs(W2(I,J,K+1))]; end  % Blood domain interaction for current voxel.
                    % Advection from K+1.
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Establish Metabolism %%%%%%%%%%%%%%%%%%%%%
                %D(Row) = D(Row) - Vol*Q_t; % Tissue domain source for source.
                % Establishes the metabolic heat source term.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Establish Pennes Perfusion %%%%%%%%%%%%%%%
                if ~obj.flow_obj.grey_white(I,J,K) % If outside of the brain domain.
                    %D(Row) = D(Row) - Vol*obj.flow_obj.arteries.domainWeighting(I,J,K)*obj.bloodCp*BloodTemp; % Tissue domain interaction for source.
                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-Vol*obj.flow_obj.arteries.domainWeighting(I,J,K)*obj.bloodCp]; % Tissue domain interaction for current voxel.
                end
                % Establishes the Pennes Perfusion source term for tissue
                % outside of the brain domain.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% 3D Inter-Domain Heat Transfer %%%%%%%%%%%%
                if ~isinf(obj.Beta.porous) % If the heat transfer coefficient is finite.
                    if Porosity_b ~= 0 % If there are two phases present.
                       T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-obj.Beta.porous]; % Tissue domain interaction for current voxel (tissue domain).
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,BRow,obj.Beta.porous]; % Tissue domain interaction for current voxel (blood domain).
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,Row,obj.Beta.porous]; % Blood domain interaction for current voxel (tissue domain).
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-obj.Beta.porous]; % Blood domain interaction for current voxel (blood domain).
                        % Establishes the inter-domain heat transfer
                        % between the tissue phase and the blood phase.
                        
                    else % If only tissue phase is present.
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,Row,2]; % Blood domain interaction for current voxel (tissue domain).
                        T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-2]; % Blood domain interaction for current voxel (blood domain).
                        % This sets the blood phase equal to the tissue
                        % phase. This avoids any voxels remaining undefined
                        % if the porosity is 0 (eg outside the brain
                        % domain).
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Transient Heat Transfer %%%%%%%%%%%%%%%%%%
                %                 if Option_TransientSolve
                                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-2*Vol*Porosity_t*Rho_t*Cp_t/Timestep];
                                    T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-2*Vol*Porosity_b*obj.flow_obj.bloodDensity*obj.bloodCp/Timestep];
                %                 end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Artery Voxel Interactions %%%%%%%%%%%%%%%%
                if ~isempty(obj.flow_obj.Volume2VesselArt{I,J,K}) % If this voxel has intersecting artery segments.
                    for mn = 1:size(obj.flow_obj.Volume2VesselArt{I,J,K},1) % For every intersecting artery segment.
                        
                        %%%%%% Preallocating Searchable Parameters %%%%%%
                        Node = obj.flow_obj.Volume2VesselArt{I,J,K}(mn,1); % Define arterial node.
                        NodeMdot = obj.flow_obj.Volume2VesselArt{I,J,K}(mn,2); % Define inter-domain mass transfer from node segment
                        Beta1_t = obj.Beta.art(Node); % Define inter-domain heat transfer coefficient from node segment.
                        % Establishing searchable parameters. This prevents
                        % constantly having to look up values in matrix.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        %%%%%% Split Condition  %%%%%%%%%%%%%%%%%%%%%%%%%
                        if obj.flow_obj.SplitArt(Node,1) == 1
                            % If flow is split along this segment, blood
                            % arrives from both nodes of the segment and no
                            % heat transfer occurs.
                            
                            NodeConn = obj.flow_obj.arteries.tree(Node,7); % Define connecting node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,NumDomRows+Node,obj.flow_obj.SplitArt(Node,2)*NodeMdot*obj.bloodCp]; % Artery domain interaction for node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,NumDomRows+NodeConn,(1-obj.flow_obj.SplitArt(Node,2))*NodeMdot*obj.bloodCp]; % Artery domain interaction for connecting node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-NodeMdot*obj.bloodCp]; % Blood domain interaction for current voxel.
                            % Inter-domain convection from mass transfer if
                            % there is a split on the segment. Here blood
                            % is taken porportionally from both nodes in
                            % the segment. No further heat transfer occurs
                            % in the segment.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            %%%%%% Non-Split Condition  %%%%%%%%%%%%%%%%%%%%%
                        else % If no split on segment.
                            
                            
                            %%%%%% Determine Direction of Flow  %%%%%%%%%%%%%
                            if obj.flow_obj.FdotArt(Node,2) < 0
                                NodeDown = obj.flow_obj.arteries.tree(Node,7); % Downstream node.
                                NodeUp = Node; % Upstream node.
                            else
                                NodeDown = Node; % Downstream node.
                                NodeUp = obj.flow_obj.arteries.tree(Node,7); % Upstream node.
                            end
                            % This determines which is the upstream and
                            % downstream nodes of the line segments. NodeUp
                            % is defined as the upstream and NodeDown is
                            % downstream.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Convection %%%%%%%%%%%%%%%%%%
                            % Blood source from leakage of arteries
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,NumDomRows+NodeUp,NodeMdot*obj.bloodCp]; % Blood domain interaction for upstream arterial node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-NodeMdot*obj.bloodCp]; % Blood domain interaction for current voxel.
                            % Inter-domain convection from mass transfer.
                            % This comes from the upstream node only.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Tissue) %%%%%%%%%
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,NumDomRows+NodeDown,Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for downstream node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,NumDomRows+NodeUp,Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for upstream node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-Porosity_t*Beta1_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,Row,Porosity_t*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,NumDomRows+NodeDown,-Porosity_t*Beta1_t*0.5]; % Downstream node interaction for downstream node.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,Row,Porosity_t*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,NumDomRows+NodeUp,-Porosity_t*Beta1_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % from tissue domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Blood) %%%%%%%%%%
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,NumDomRows+NodeDown,Porosity_b*Beta1_t*0.5]; % Blood domain interaction for downstream node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Porosity_b*Beta1_t*0.5]; % Blood domain interaction for upstream node.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,NumDomRows+NodeUp,Porosity_b*Beta1_t*0.5]; % Blood domain interaction for upstream node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Porosity_b*Beta1_t*0.5]; % Blood domain interaction for current voxel.
                            
                           T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,BRow,Porosity_b*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,NumDomRows+NodeDown,-Porosity_b*Beta1_t*0.5]; % Downstream node interaction for downstream node.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,BRow,Porosity_b*Beta1_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+NodeDown,NumDomRows+NodeUp,-Porosity_b*Beta1_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % from blood domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Vein Voxel Interactions %%%%%%%%%%%%%%%%%%
                if ~isempty(obj.flow_obj.Volume2VesselVein{I,J,K}) % If this voxel has intersecting vein segments.
                    for mn = 1:size(obj.flow_obj.Volume2VesselVein{I,J,K},1) % For every intersecting vein segment.
                        
                        %%%%%% Preallocating Searchable Parameters %%%%%%
                        Node = obj.flow_obj.Volume2VesselVein{I,J,K}(mn,1);
                        Beta2_t = obj.Beta.vein(Node);
                        % Establishing searchable parameters. This prevents
                        % constantly having to look up values in matrix.
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                        
                        %%%%%% Non-Split Condition  %%%%%%%%%%%%%%%%%%%%%
                        if obj.flow_obj.SplitVein(Node,1) ~= 1
                            % If flow is split along this segment, no heat
                            % transfer occurs.
                            
                            
                            %%%%%% Determine Direction of Flow  %%%%%%%%%%%%%
                            if obj.flow_obj.FdotVein(Node,2) < 0
                                NodeDown = obj.flow_obj.veins.tree(Node,7);
                                NodeUp = Node;
                            else
                                NodeDown = Node;
                                NodeUp = obj.flow_obj.veins.tree(Node,7);
                            end
                            % This determines which is the upstream and
                            % downstream nodes of the line segments. NodeUp
                            % is defined as the upstream and NodeDown is
                            % downstream.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Tissue) %%%%%%%%%
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,NumDomRows+VesselRow+NodeDown,Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for downstream node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,NumDomRows+VesselRow+NodeUp,Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for upstream node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [Row,Row,-Porosity_t*Beta2_t*0.5]; % Tissue domain interaction for current voxel.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,Row,Porosity_t*Beta2_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeDown,-Porosity_t*Beta2_t*0.5]; % Downstream node interaction for downstream node.
                            %
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,Row,Porosity_t*Beta2_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeUp,-Porosity_t*Beta2_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % from tissue domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                            %%%%%% Inter-Domain Conduction (Blood) %%%%%%%%%%
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,NumDomRows+VesselRow+NodeDown,Porosity_b*Beta2_t*0.5]; % Blood domain interaction for downstream node.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Porosity_b*Beta2_t*0.5]; % Blood domain interaction for upstream node.
                            
                           T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,NumDomRows+VesselRow+NodeUp,Porosity_b*Beta2_t*0.5]; % Blood domain interaction for upstream node.
                           T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [BRow,BRow,-Porosity_b*Beta2_t*0.5]; % Blood domain interaction for current voxel.
                            
                           T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,BRow,Porosity_b*Beta2_t*0.5];  % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeDown,-Porosity_b*Beta2_t*0.5]; % Downstream node interaction for downstream node.
                            
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,BRow,Porosity_b*Beta2_t*0.5]; % Downstream node interaction for current voxel.
                            T_DataTempPos = T_DataTempPos + 1; T_DataTemp(T_DataTempPos,:) = [NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeUp,-Porosity_b*Beta2_t*0.5]; % Downstream node interaction for upstream node.
                            % Inter-domain heat transfer to the segment
                            % Inter-domain heat transfer to the segment
                            % from blood domain. Transfers heat to the
                            % downstream node based on the average
                            % tempeartures of the upstream and downstream
                            % nodes.
                            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                            
                            
                        end
                        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                        
                    end
                end
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                

                
            end
        end
    end
    if T_DataTempPos > 3e8
        fprintf('Warning Preallocated Length has been exceeded');
    end
end
                %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
                T_DataTemp = T_DataTemp(1:T_DataTempPos, :);
                T_Super(T_DataRowcount:T_DataRowcount-1+size(T_DataTemp,1),:) = T_DataTemp; % Add T_DataTemp to T_Super.
                T_DataRowcount = T_DataRowcount + size(T_DataTemp,1); % Update T_DataRowcount.
                T_DataTemp = []; % Clear T_DataTemp
                % Compiling interactions into T_Super list.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Artery Segment Interactions %%%%
disp('Compiling Interactions within Arteries')
tic % Start timing.

for Node = 1:size(obj.flow_obj.Vessel2VolumeArt,1)
    
    %%%%%% Inlet Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    if ismember(Node,obj.flow_obj.arteries.termPoints)
        %Index = find(Node==obj.flow_obj.arteries.termPoints); % Find which inlet is being referred to.
        T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node,2]; % Set node equal to InletTemp.
        %D(NumDomRows+Node) = InletTemp(Index); % Set node equal to InletTemp.
        % If node is an inlet, set temperature to corresponding InletTemp.
        
        if ~obj.flow_obj.SplitArt(Node,1) % If current segment is not split.
            NodeConn = obj.flow_obj.arteries.tree(Node,7); % Establish connecting node.
            if obj.flow_obj.FdotArt(Node,1) < 0 % If flow is directed towards connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+Node,abs(obj.flow_obj.FdotArt(Node,1))*obj.bloodCp]; % Connecting node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+NodeConn,-abs(obj.flow_obj.FdotArt(Node,1))*obj.bloodCp]; % Connecting node interaction for connecting node.
            end
        end
        % Do advection flow for the downstream node of current segment.
        % Flow cannot be directed towards an inlet so no upstream advection
        % is performed.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%% Non-Inlet Nodes %%%%%%%%%%%%%%%%%%%%%%%%%%
    else
        
        %%%%%% Arterial Advection Interactions %%%%%%%%%%
        if Node~=1 % Node 1 does not contain a segment.
            if ~obj.flow_obj.SplitArt(Node,1) % If current segment is not split.
                NodeConn = obj.flow_obj.arteries.tree(Node,7); % Establish connecting node.
                if obj.flow_obj.FdotArt(Node,2) > 0 % If flow is directed towards current node.
                    T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+NodeConn,abs(obj.flow_obj.FdotArt(Node,2))*obj.bloodCp]; % Current node interaction for connecting node.
                    T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node,-abs(obj.flow_obj.FdotArt(Node,2))*obj.bloodCp]; % Current node interaction for current node.
                end
                if obj.flow_obj.FdotArt(Node,1) < 0 % If flow is directed towards connecting node.
                    T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+Node,abs(obj.flow_obj.FdotArt(Node,1))*obj.bloodCp]; % Connecting node interaction for current node.
                    T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+NodeConn,-abs(obj.flow_obj.FdotArt(Node,1))*obj.bloodCp]; % Connecting node interaction for connecting node.
                end
            end
        end
        % Do advection flow for the downstream node of current segment.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        %%%%%% Arterial Conduction Interactions %%%%%%%%%
        if Node~=1 % Node 1 does not contain a segment.
            NodeConn = obj.flow_obj.arteries.tree(Node,7); % Establish connecting node.
            if ~(NodeConn == 1 && ismember(1,obj.flow_obj.arteries.termPoints)) % Imbalance occurs if there is conduction to an inlet point.
                T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+NodeConn, obj.flow_obj.geometry.art.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.art.L(Node))]; % Current node interaction for connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node, -obj.flow_obj.geometry.art.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.art.L(Node))]; % Current node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+Node, obj.flow_obj.geometry.art.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.art.L(Node))]; % Connecting node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+NodeConn,NumDomRows+NodeConn, -obj.flow_obj.geometry.art.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.art.L(Node))]; % Connecting node interaction for current node.
                
            end
        end
        % Do conduction for the segment. Both current and connecting nodes
        % are updated here so that connections do not have to be searched
        % for in each node.
        % Conduction should not happen to an inlet point as it is defined
        % directly as the inlet temperature.
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        
        
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
        %%%%%% Transient Heat Transfer %%%%%%%%%%%%%%%%%%
        %     if Option_TransientSolve
                T_DataTemp = [T_DataTemp;NumDomRows+Node,NumDomRows+Node,-2*obj.flow_obj.geometry.art.V(Node)*obj.flow_obj.bloodDensity*obj.bloodCp/Timestep];
        %     end
        %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        
    end
    
    %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
    if size(T_DataTemp,1) ~= 0
        T_Super(T_DataRowcount:T_DataRowcount-1+size(T_DataTemp,1),:) = T_DataTemp; % Add T_DataTemp to T_Super.
        T_DataRowcount = T_DataRowcount + size(T_DataTemp,1); % Update T_DataRowcount.
        T_DataTemp = []; % Clear T_DataTemp.
        % Compiling interactions into T_Super list.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Venous Segment Interactions %%%%
disp('Compiling Interactions within Veins')
tic % Start timing.

% Venous Vessel Heat Flow
for Node = 1:size(obj.flow_obj.Vessel2VolumeVein,1)
    
    %%%%%% Venous Advection Interactions %%%%%%%%%%%%
    if Node~=1 % Node 1 does not contain a segment.
        if ~obj.flow_obj.SplitVein(Node,1) % If current segment is not split.
            NodeConn = obj.flow_obj.veins.tree(Node,7); % Establish connecting node.
            if obj.flow_obj.FdotVein(Node,1) > 0 % If flow is directed towards current node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+NodeConn,abs(obj.flow_obj.FdotVein(Node,1))*obj.bloodCp]; % Current node interaction for connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node,-abs(obj.flow_obj.FdotVein(Node,1))*obj.bloodCp]; % Current node interaction for current node.
            end
            if obj.flow_obj.FdotVein(Node,2) < 0 % If flow is directed towards connecting node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+Node,abs(obj.flow_obj.FdotVein(Node,2))*obj.bloodCp]; % Connecting node interaction for current node.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+NodeConn,-abs(obj.flow_obj.FdotVein(Node,2))*obj.bloodCp]; % Connecting node interaction for current node.
            end
        end
    end
    % Do advection flow for the downstream node of current segment.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Venous Conduction Interactions %%%%%%%%%%%
    if Node ~= 1 % Node 1 does not contain a segment.
        NodeConn = obj.flow_obj.veins.tree(Node,7); % Establish connecting node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+NodeConn, obj.flow_obj.geometry.vein.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.vein.L(Node))]; % Current node interaction for connecting node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node, -obj.flow_obj.geometry.vein.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.vein.L(Node))]; % Current node interaction for current node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+Node, obj.flow_obj.geometry.vein.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.vein.L(Node))]; % Connecting node interaction for current node.
        T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+NodeConn, -obj.flow_obj.geometry.vein.Aavg(Node)*(obj.bloodKc/obj.flow_obj.geometry.vein.L(Node))]; % Connecting node interaction for current node.
    end
    % Do conduction for the segment. Both current and connecting nodes
    % are updated here so that connections do not have to be searched
    % for in each node.
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Inter-Domain Convection %%%%%%%%%%%%%%%%%%
    if Node~=1 % Node 1 does not contain a segment.
        
        if obj.flow_obj.FdotVein(Node,1) > 0
            NodeDown = Node;
        else
            NodeDown = NodeConn;
        end
        % This determines which is the downstream node (NodeDown) of the
        % line segments.
        
        for NumInts = 1:size(obj.flow_obj.Vessel2VolumeVein{Node},1) % For all intersections of venous segment.
            I = obj.flow_obj.Vessel2VolumeVein{Node}(NumInts,1); % Get location I for ntersecting voxel.
            J = obj.flow_obj.Vessel2VolumeVein{Node}(NumInts,2); % Get location J for ntersecting voxel.
            K = obj.flow_obj.Vessel2VolumeVein{Node}(NumInts,3); % Get location K for ntersecting voxel.
            BRow = DomTotConvert(I,J,K) + RowAdjustment; % Find corresponding row value for intersecting voxel (blood domain).
            VoxMdot = abs(obj.flow_obj.Vessel2VolumeVein{Node}(NumInts,4)); % Get mass transfer for ntersecting voxel.
            
            if obj.flow_obj.SplitVein(Node,1) % If segment is split.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,BRow,(obj.flow_obj.SplitVein(Node,2))*VoxMdot*obj.bloodCp]; % Current node interaction for blood domain.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node,-(obj.flow_obj.SplitVein(Node,2))*VoxMdot*obj.bloodCp]; % Current node interaction for current node.
                
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,BRow,(1-obj.flow_obj.SplitVein(Node,2))*VoxMdot*obj.bloodCp]; % Connecting node interaction for blood domain.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeConn,NumDomRows+VesselRow+NodeConn,-(1-obj.flow_obj.SplitVein(Node,2))*VoxMdot*obj.bloodCp]; % Connecting node interaction for connecting node.
                % Deliver the corresponding amount of blood to both nodes
                % on the segment.
                
            else % If segment is not split.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,BRow,VoxMdot*obj.bloodCp]; % Downstream node interaction for blood domain.
                T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+NodeDown,NumDomRows+VesselRow+NodeDown,-VoxMdot*obj.bloodCp]; % Downstream node interaction for downstream node.
                % Deliver the blood to the downstream node only.
            end
        end
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Transient Heat Transfer %%%%%%%%%%%%%%%%%%
    %     if Option_TransientSolve
            T_DataTemp = [T_DataTemp;NumDomRows+VesselRow+Node,NumDomRows+VesselRow+Node,-2*obj.flow_obj.geometry.vein.V(Node)*obj.flow_obj.bloodDensity*obj.bloodCp/Timestep];
    %     end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
    %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
    if size(T_DataTemp,1) ~= 0
        T_Super(T_DataRowcount:T_DataRowcount-1+size(T_DataTemp,1),:) = T_DataTemp; % Add T_DataTemp to T_Super.
        T_DataRowcount = T_DataRowcount + size(T_DataTemp,1); % Update T_DataRowcount.
        T_DataTemp = []; % Clear T_DataTemp.
        % Compiling interactions into T_Super list.
    end
    %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    
    
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Creating Sparse Matrix %%%%%%%%%%%%%%%%%%%
T_Boolean = logical(T_Super(:,1));
T_Super = T_Super(T_Boolean,:); % Delete any excess rows from T_Super.
T_Super(:,3) = 0.5*T_Super(:,3);
T_Solve = sparse(T_Super(:,1),T_Super(:,2),T_Super(:,3)); % Convert P_Super into sparse matrix.
clearvars T_Super T_Boolean % Delete variables to free up memory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

obj.T_Solve = T_Solve;
%obj.D = D;

end

