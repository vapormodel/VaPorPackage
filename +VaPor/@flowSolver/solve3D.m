function [obj] = solve3D(obj)
%SOLVE3D Summary of this function goes here
%   Detailed explanation goes here
% Solving the flowrates in the 3D porous region.


disp('%%%%%%% Solving 3D Flowrates %%%%%%%%%%%%%%%%%%%%%%') % Display.


%%%%%% Setting up Reference Lists for obj.grey_white %%%%%
row_convert = find(obj.grey_white);  % goes from row number to I,J,K
grey_white_convert = zeros(size(obj.grey_white));
grey_white_convert(row_convert) = 1:numel(obj.grey_white(obj.grey_white)); % goes from I,J,K to row number
% The sparse matrix to solve only needs those voxels that are within obj.grey_white.
% These lists help convert sparse matrix row number to location in obj.grey_white
% (given by I,J,K) and from locations in obj.grey_white to row numbers in the sparse
% matrix.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Initialising Storage Matrices %%%%%%%%%%%%
P_DataTemp = [];
P_DataRowcount = 1;
P_Super = zeros(1e8,3);
D = zeros(numel(obj.grey_white(obj.grey_white)),1);
% P_super will be a list of all interactions that will be fed into the
% sparse matrix stored as [row, col, value]. To create the list, it is fed
% by P_DataTemp on every iteration through the domain which is filled up
% by interactions. P_DataRowcount keeps track of where to add P_DataTemp
% at the end of P_Super.
% D contains all the values that are not variable dependent (in this case
% all the obj.MdotVoxelsOverall values). This will be the RHS of the equation in the linear
% solver (Mx = D where M is the sparse matrix, x is the list of variables)
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Compiling Interactions %%%%%%%%%%%%%%%%%%%
disp('Compiling Interactions within Flow Domain') % Display.
tic % Start timing.

Check = 0; % Used as a condition to set the first voxel to P=0.
for I = 1:size(obj.grey_white,1)
    for J = 1:size(obj.grey_white,2)
        for K = 1:size(obj.grey_white,3)
            if obj.grey_white(I,J,K)
                % For all voxels, if they are within obj.grey_white.
                
                %%%%%% Set Pressure Condition %%%%%%%%%%%%%%%%%%%
                if Check == 0
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K),1]; % Set P(I,J,K)=0;
                    Check = 1; % Changes condition.
                end
                % Sets the first value to P = 0 so that the solution is
                % bounded.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Mass Source Term %%%%%%%%%%%%%%%%%%%%%%%%%
                D(grey_white_convert(I,J,K)) = D(grey_white_convert(I,J,K)) + obj.MdotVoxelsOverall(I,J,K);
                % Source term from inter-domain mass transfer.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Interactions with Neighbouring Voxels %%%%
                if I > 1 && obj.grey_white(I-1,J,K)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I-1,J,K)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I-1,J,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at I-1.
                
                if I < size(obj.grey_white,1) && obj.grey_white(I+1,J,K)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I+1,J,K)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I+1,J,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at I+1.
                
                if J > 1 && obj.grey_white(I,J-1,K)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I,J-1,K)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J-1,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at J-1.
                
                if J < size(obj.grey_white,2) && obj.grey_white(I,J+1,K)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I,J+1,K)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J+1,K),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at J+1.
                
                if K > 1 && obj.grey_white(I,J,K-1)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I,J,K-1)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K-1),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at K-1.
                
                if K < size(obj.grey_white,3) && obj.grey_white(I,J,K+1)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I,J,K+1)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K+1),-G]; % Interaction for neighbouring voxel.
                    P_DataTemp = [P_DataTemp;grey_white_convert(I,J,K),grey_white_convert(I,J,K),G]; % Interaction for current voxel.
                end
                % Exchange with voxel at K+1.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
                %%%%%% Updating Sparse Matrix List %%%%%%%%%%%%%%
                P_Super(P_DataRowcount:P_DataRowcount-1+size(P_DataTemp,1),:) = P_DataTemp; % Add P_DataTemp to P_Super.
                P_DataRowcount = P_DataRowcount + size(P_DataTemp,1); % Update P_DataRowcount.
                P_DataTemp = []; % Clear P_DataTemp.
                % Compiling interactions into P_Super list.
                %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
                
                
            end
        end
    end
end
toc % Finish timing.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Creating Sparse Matrix %%%%%%%%%%%%%%%%%%%
P_Boolean = logical(P_Super(:,1));
P_Super = P_Super(P_Boolean,:); % Delete any excess rows from P_Super.
P_Solve = sparse(P_Super(:,1),P_Super(:,2),P_Super(:,3)); % Convert P_Super into sparse matrix.
clearvars P_super P_Boolean % Delete variables to free up memory.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Solving Linear System %%%%%%%%%%%%%%%%%%%%
disp('Starting Linear Solver') % Display.
disp(['Solving Matrix Inversion size: ' num2str(size(P_Solve,1)) ' by ' num2str(size(P_Solve,2))]) % Display.
tic % Start timing.

P_N = P_Solve\D; % Solving the linear system

toc % Finish timing.
disp('Matrix Inversion Completed')
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Reshaping Resutls from Linear Solver %%%%%
P = zeros(size(obj.grey_white));
P(obj.grey_white) = P_N;
P(~obj.grey_white) = NaN;
P = P-min(min(min(P)));
% Converting the results back into the same shape as obj.grey_white domain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Deriving Velocities %%%%%%%%%%%%%%%%%%%%%%
obj.porous.U=zeros(size(obj.grey_white,1)+1,size(obj.grey_white,2),size(obj.grey_white,3));
obj.porous.V=zeros(size(obj.grey_white,1),size(obj.grey_white,2)+1,size(obj.grey_white,3));
obj.porous.W=zeros(size(obj.grey_white,1),size(obj.grey_white,2),size(obj.grey_white,3)+1);
% Initialise velocity matrices.

for I = 1:size(obj.grey_white,1)
    for J = 1:size(obj.grey_white,2)
        for K = 1:size(obj.grey_white,3)
            if obj.grey_white(I,J,K)
                % For all voxels, if they are within obj.grey_white.
                
                if I < size(obj.grey_white,1) && obj.grey_white(I+1,J,K)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I+1,J,K)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    obj.porous.U(I+1,J,K) = 1/obj.bloodDensity*1/obj.voxelSize^2*G*(P(I,J,K)-P(I+1,J,K)); % Calculate velocity.
                end
                % create U velocity at voxel boundary I-1.
                
                if J < size(obj.grey_white,2) && obj.grey_white(I,J+1,K)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I,J+1,K)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    obj.porous.V(I,J+1,K) = 1/obj.bloodDensity*1/obj.voxelSize^2*G*(P(I,J,K)-P(I,J+1,K)); % Calculate velocity.
                end
                % create V velocity at voxel boundary J-1.
                
                if K < size(obj.grey_white,3) && obj.grey_white(I,J,K+1)
                    AvgPorosity = 0.5*(obj.porosity(I,J,K) + obj.porosity(I,J,K+1)); % Average porosity of two voxels.
                    G = obj.bloodDensity*obj.voxelSize*AvgPorosity*pi*obj.capDia^2/(32*obj.bloodVisc*obj.tortuosity); % Calculate conductance.
                    obj.porous.W(I,J,K+1) = 1/obj.bloodDensity*1/obj.voxelSize^2*G*(P(I,J,K)-P(I,J,K+1)); % Calculate velocity.
                end
                % create W velocity at voxel boundary K-1.
                
            end
        end
    end
end
% Velocities are derived from the pressure differences of each voxel. Only
% the faces at I-1, J-1 & K-1 need to be done for each voxel. Velocities at
% boudries are zero by default.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Calculate Derived Perfusion %%%%%%%%%%%%%%
if ~obj.options.counterCurrentFlow 
    measurePerfusion
    obj.measuredPerfusion = MeasuredPerfusion;
    obj.massFlow = MassFlow;
end
% This calculates perfusion values from flows into every voxel. Allows
% comparison with predicted Perfusion values inputted into the model.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


%%%%%% Creating Colocated Velocities %%%%%%%%%%%%
obj.porous.UU = 0.5*(obj.porous.U(1:size(obj.grey_white,1),1:size(obj.grey_white,2),1:size(obj.grey_white,3))+obj.porous.U(2:size(obj.grey_white,1)+1,1:size(obj.grey_white,2),1:size(obj.grey_white,3)));
obj.porous.VV = 0.5*(obj.porous.V(1:size(obj.grey_white,1),1:size(obj.grey_white,2),1:size(obj.grey_white,3))+obj.porous.V(1:size(obj.grey_white,1),2:size(obj.grey_white,2)+1,1:size(obj.grey_white,3)));
obj.porous.WW = 0.5*(obj.porous.W(1:size(obj.grey_white,1),1:size(obj.grey_white,2),1:size(obj.grey_white,3))+obj.porous.W(1:size(obj.grey_white,1),1:size(obj.grey_white,2),2:size(obj.grey_white,3)+1));
obj.porous.Velocity = sqrt(obj.porous.UU.^2 + obj.porous.VV.^2 + obj.porous.WW.^2);
% The velocities used in the solver are stored at voxel boundaries which
% are difficult to visualise alongside data stored at voxel centres.
% Therefore combining the average values of two faces gives an
% approximation for values at the centre. These are only used for display
% purposes.

obj.porous.UU(~obj.grey_white) = NaN;
obj.porous.VV(~obj.grey_white) = NaN;
obj.porous.WW(~obj.grey_white) = NaN;
obj.porous.Velocity(~obj.grey_white) = NaN;
% Setting all values outside of the domain to NaN. This facilitates
% visualisation of data with 'planecut'.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


disp('%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%') % Display.

end

