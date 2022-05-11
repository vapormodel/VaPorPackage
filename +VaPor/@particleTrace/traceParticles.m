function [obj] = traceParticles(obj, numParticles, parallel, overrideInlet)
%TRACEPARTCILES Summary of this function goes here
%   Detailed explanation goes here
% initialise concentration measurements
if nargin < 3
    parallel = false;
end
Timestep = 0.1;
SaveTime = Timestep;
TotalTime = 120;
NoTimesteps = ceil(TotalTime/Timestep);
TotalTime = NoTimesteps*Timestep;

InjectionLength = 1;
T_TransientStore = zeros(size(obj.Tree,1),TotalTime/SaveTime);

disp('Simulating Particles')
tic

ParticlesGenerated = numParticles;


Option_Streamlines = true;
StreamlineFigureTrigger = true;

timelist = zeros(ParticlesGenerated,4);
if nargin == 4
    Option_OverrideInlet = true;
    newInlet = overrideInlet;
    obj.new_inlet = newInlet;
else
    Option_OverrideInlet = false;
    obj.new_inlet = 0;
end
%%% Check if the script should
if ~parallel
    T_TransientStore = particleTracerParInner(ParticlesGenerated,obj.flowobj.arteries.termPoints,obj.flowobj.arteries.termFlows,...
        obj.flowobj.veins.termPoints, obj.NumGM_WM, obj.NumArt, obj.Tree, TotalTime, Timestep,...
        SaveTime, InjectionLength, false, 1);
    
else
    numberOfBatches = ceil(ParticlesGenerated/1e6)*4;
    fprintf('Using %d batches. \n', numberOfBatches);
    particlesToInsert = ceil(ParticlesGenerated/numberOfBatches);
    for run=1:numberOfBatches
        if ~Option_OverrideInlet
            F(run) = parfeval(@particleTracerParInner, 1, particlesToInsert,obj.flowobj.arteries.termPoints,obj.flowobj.arteries.termFlows,...
                obj.flowobj.veins.termPoints, obj.NumGM_WM, obj.NumArt, obj.Tree, TotalTime, Timestep,...
                SaveTime, InjectionLength, false, 1);
        else
            F(run) = parfeval(@particleTracerParInner, 1, particlesToInsert,newInlet,-obj.flowobj.FdotArt(newInlet),...
                obj.flowobj.veins.termPoints, obj.NumGM_WM, obj.NumArt, obj.Tree, TotalTime, Timestep,...
                SaveTime, InjectionLength, false, 1);
        end
    end
    %Build a waitbar to track progress
    h = waitbar(0,'Simulating Particles in Parallel...');
    
    for run=1:numberOfBatches
        [completedRun,thisResult] = fetchNext(F);
        T_TransientStore = T_TransientStore + thisResult;
        waitbar(run/numberOfBatches,h, sprintf('Simulated %d Particles So Far', run*particlesToInsert));
        
        % Discard the read result
        F(completedRun) = [];
    end
    close(h);
end

% T_TransientStore = [zeros(size(Tree,1),1), [T_TransientStore(NumArt+1:NumArt+NumGM_WM,:);T_TransientStore(1:NumArt,:);T_TransientStore(NumArt+NumGM_WM+1:end,:)]];
%obj.T_TransientStore_Save = T_TransientStore;

T_TransientStore = T_TransientStore/ParticlesGenerated;
%
for N = 1:size(T_TransientStore,2)
    T_TransientStore(1:obj.NumGM_WM,N) = T_TransientStore(1:obj.NumGM_WM,N)./((obj.flowobj.voxelSize^3)*obj.flowobj.porosity(obj.flowobj.grey_white));
    T_TransientStore(obj.NumGM_WM+2:obj.NumGM_WM+obj.NumArt,N) = T_TransientStore(obj.NumGM_WM+2:obj.NumGM_WM+obj.NumArt,N)./obj.flowobj.geometry.art.V(2:end);
    T_TransientStore(obj.NumGM_WM+1,N) = T_TransientStore(obj.NumGM_WM+2,N);
    T_TransientStore(obj.NumGM_WM+obj.NumArt+2:end,N) = T_TransientStore(obj.NumGM_WM+obj.NumArt+2:end,N)./obj.flowobj.geometry.vein.V(2:end);
    T_TransientStore(obj.NumGM_WM+obj.NumArt+1,N) = T_TransientStore(obj.NumGM_WM+obj.NumArt+2,N);
end

T_TransientStore = [zeros(size(obj.Tree,1),1), T_TransientStore];
% T_TransientStore(NumGM_WM+InletPoints,:) = [ones(1,round(InjectionLength/Timestep)+1), zeros(1,round(TotalTime/Timestep-(InjectionLength/Timestep)));...
%     ones(1,round(InjectionLength/Timestep)+1), zeros(1,round(TotalTime/Timestep-(InjectionLength/Timestep)));...
%     ones(1,round(InjectionLength/Timestep)+1), zeros(1,round(TotalTime/Timestep-(InjectionLength/Timestep)))];

if ~Option_OverrideInlet
    AIF = mean([T_TransientStore(obj.NumDomTot+obj.flowobj.arteries.termPoints(1),2:end); T_TransientStore(obj.NumDomTot+obj.flowobj.arteries.termPoints(2),2:end); T_TransientStore(obj.NumDomTot+obj.flowobj.arteries.termPoints(3),2:end)]);
else
    AIF = T_TransientStore(obj.NumDomTot+newInlet+1,2:end);
end
obj.T_TransientStore = T_TransientStore/max(AIF);
end

