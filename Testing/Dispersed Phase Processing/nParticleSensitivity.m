%% Preamble

clear variables;
close all;
clc;
evalc('delete(gcp(''nocreate''));');

if exist('/mnt/Processing/Data', 'dir')
    saveLocation = '/mnt/Processing/Data';
else
    saveLocation = '~/Data';
end

nProc = maxNumCompThreads - 2; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%%

load('~/Downloads/Windsor_SB_wW.mat');
load('/mnt/Processing/Data/Numerical/MATLAB/LagData/Windsor_SB_wW_Upstream_SC/plane/T12500_T40000_F200_cumulative.mat');

LagData = LagData.X_P020225;
% LagData = LagData.X_P124625;

% removeNthParcel = 2;
keepNthParcel = 3;

parcelMass = 8e-11;

nParcels = 75e6;

cellSize.target = 8e-3;

normalise = false;

%%

% nParticleAll = cell2mat(LagData.nParticle);
% dAll = cell2mat(LagData.d);
% origIdAll = cell2mat(LagData.origId);
% 
% massAll = 1000 * (nParticleAll .* ((tau * dAll.^3) / 12));


%%

evalc('parpool(''threads'');');

%%%%

nTimes = height(LagData.time);

% Identify Empty Time Instances Not Related to Case Initialisation
emptyTimes = cellfun(@isempty, LagData.origId);
firstValidTime = find(emptyTimes == false, 1, 'first');

suspectTimes = find(emptyTimes(firstValidTime:end) == true) + (firstValidTime - 1);

if ~isempty(suspectTimes)
    
    for i = 1:height(suspectTimes)
        disp(['        WARNING: Time ''', num2str(LagData.time(suspectTimes(i))), ''' Is Unexpectedly Empty']);
    end
    clear i;
    
end

% Shift Data Origin
for i = 1:nTimes

    if ~isempty(LagData.positionCartesian{i})
        LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
    end

end
clear i;

% Specify Map Boundaries
mapPerim = [];

xLimsData = double(LagData.positionCartesian{end}(1,1));
yLimsData = [-0.522; 0.522];
zLimsData = [0; 0.6264];

% Collate Particles of Interest
index = cell(nTimes,1);
for i = 1:nTimes
            
    if ~isempty(LagData.positionCartesian{i})
        index{i} = find((LagData.positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                        (LagData.positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                        (LagData.positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                        (LagData.positionCartesian{i}(:,3) <= zLimsData(2)));
    end

end
clear i;

% Remove Unnecessary Data
LagFields = fieldnames(LagData);
reqFields = {'time'; 'd'; 'nParticle'; 'origId'; 'positionCartesian'};

LagData = rmfield(LagData, LagFields(~ismember(LagFields, reqFields)));
LagProps = LagProps(ismember(LagProps, reqFields));

for i = 1:nTimes
    
    for j = 1:height(LagProps)
        LagData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
    end
    clear j;

end
clear i;

% Artificially Reduce Parcel Count
% removeIDs = (1:removeNthParcel:nParcels)';
removeIDs = (keepNthParcel:keepNthParcel:nParcels)'; removeIDs = setdiff((1:nParcels)', removeIDs);

nParcelsRemoved = height(removeIDs);
massRemoved = nParcelsRemoved * parcelMass;
massToAdd = massRemoved / (nParcels - nParcelsRemoved);
parcelMass = parcelMass + massToAdd;

for i = 1:nTimes
    LagData.nParticle{i} = parcelMass ./ (1000 * ((tau * LagData.d{i}.^3) / 12));
end

if height(removeIDs) > 0
    
    tic;
    for i = 1:nTimes
        index = ismember(LagData.origId{i}, removeIDs);
        
        for j = 1:height(LagProps)
            LagData.(LagProps{j}){i}(index,:) = [];
        end
        clear j;
    
    end
    clear i;
    toc;
    
end
        
% Generate Presentation Grid
cellSize.x = cellSize.target;
cellSize.y = (yLimsData (2) - yLimsData (1)) / round(((yLimsData (2) - yLimsData (1)) / cellSize.target));
cellSize.z = (zLimsData (2) - zLimsData (1)) / round(((zLimsData (2) - zLimsData (1)) / cellSize.target));

cellSize.area = cellSize.y * cellSize.z;

[y, z] = ndgrid(yLimsData(1):cellSize.y:yLimsData(2), zLimsData(1):cellSize.z:zLimsData(2));

mapData.positionGrid = zeros([height(y(:)),3]);
mapData.positionGrid(:,1) = xLimsData;
mapData.positionGrid(:,(2:3)) = [y(:), z(:)];

% Assign Particles to Grid Cells
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Closest Mesh Node

positionGrid = mapData.positionGrid;
positionCartesian = LagData.positionCartesian;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index{i} = dsearchn(positionGrid, positionCartesian{i});
    end

    % Remove Unnecessary Data
    positionCartesian{i} = [];
end
clear i positionGrid positionCartesian;

% Generate Instantaneous Contaminant Maps
mapData.time = LagData.time;

% Calculate Instantaneous Field Variables
nParticles = cell(nTimes,1); nParticles(:) = {zeros([height(mapData.positionGrid),1], 'single')}; % Number of Particles in Cell
mass = nParticles; % Total Mass in Cell
d32 = nParticles; % Sauter Mean Diameter in Cell
dTemp = zeros([height(mapData.positionGrid),1], 'single');
nParticle = LagData.nParticle;
d = LagData.d;
parfor i = 1:nTimes
    d30Temp = dTemp;
    d20Temp = dTemp;
    
    if totalParticles(i) > 0
        
        for j = 1:totalParticles(i)
            nParticles{i}(index{i}(j)) = nParticles{i}(index{i}(j)) + ...
                                         nParticle{i}(j);
            
            mass{i}(index{i}(j)) = mass{i}(index{i}(j)) + ...
                                   (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
            
            d30Temp(index{i}(j)) = d30Temp(index{i}(j)) + ...
                                   (nParticle{i}(j) * (d{i}(j)^3));
            
            d20Temp(index{i}(j)) = d20Temp(index{i}(j)) + ...
                                   (nParticle{i}(j) * (d{i}(j)^2));
        end
        
    end
    
    % Remove Unnecessary Data
    index{i} = [];
    nParticle{i} = [];
    d{i} = [];

    % Calculate Derived Variables
    mass{i} = 1000 * mass{i};
    d32{i} = (d30Temp ./ d20Temp) * 1e6;
    
    % Set Empty Cells Back to Zero
    d32{i}(isnan(d32{i})) = 0;
end
clear i j dTemp nParticle d;

clear LagData;

mapData.inst.nParticles = nParticles; clear nParticles;
mapData.inst.mass = mass; clear mass;
mapData.inst.d32 = d32; clear d32;

% Calculate Instantaneous Centre of Mass
mapData.inst.CoM = cell(nTimes,1); mapData.inst.CoM(:) = {zeros([1,3], 'single')};

for i = 1:nTimes
    mapData.inst.CoM{i}(1) = mapData.positionGrid(1,1);
    mapData.inst.CoM{i}(2) = sum(mapData.inst.mass{i} .* mapData.positionGrid(:,2)) / sum(mapData.inst.mass{i});
    mapData.inst.CoM{i}(3) = sum(mapData.inst.mass{i} .* mapData.positionGrid(:,3)) / sum(mapData.inst.mass{i});
end
clear i;

% Generate Time-Averaged Contaminant Map
mapData.mean.nParticles = zeros([height(mapData.positionGrid),1], 'single');
mapData.mean.mass = mapData.mean.nParticles;
mapData.mean.d32 = mapData.mean.nParticles;

for i = 1:nTimes
    mapData.mean.nParticles = mapData.mean.nParticles + mapData.inst.nParticles{i};
    mapData.mean.mass = mapData.mean.mass + mapData.inst.mass{i};
    mapData.mean.d32 = mapData.mean.d32 + mapData.inst.d32{i};
end
clear i;

mapData.mean.nParticles = mapData.mean.nParticles / nTimes;
mapData.mean.mass = mapData.mean.mass / nTimes;
mapData.mean.d32 = mapData.mean.d32 / nTimes;

% Calculate RMS
mapData.inst.massPrime = mapData.inst.mass;
mapData.mean.massRMS = zeros([height(mapData.mean.mass),1]);

for i = 1:nTimes
    mapData.inst.massPrime{i} = mapData.inst.massPrime{i} - mapData.mean.mass;
    mapData.mean.massRMS = mapData.mean.massRMS + mapData.inst.massPrime{i}.^2;
end
clear i;

mapData.mean.massRMS = sqrt((1 / nTimes) * mapData.mean.massRMS);

% Calculate Time-Averaged Centre of Mass
mapData.mean.CoM = zeros([1,3], 'single');
    
mapData.mean.CoM(1) = mapData.positionGrid(1,1);
mapData.mean.CoM(2) = sum(mapData.mean.mass .* mapData.positionGrid(:,2)) / sum(mapData.mean.mass);
mapData.mean.CoM(3) = sum(mapData.mean.mass .* mapData.positionGrid(:,3)) / sum(mapData.mean.mass);

%%%%

evalc('delete(gcp(''nocreate''));');


%%

plotVars = fieldnames(mapData.mean); plotVars = plotVars(1:(end - 1));

orientation = 'YZ';
xLimsPlot = [0.31875; 2.73575];
yLimsPlot = [-0.522; 0.522];
zLimsPlot = [0; 0.6264];
    
positionData = mapData.positionGrid;
nPlanes = 1;
planeNo = 1;
cMap = flipud(viridis(32));
refPoint = [];
figTitle = '-'; % Leave Blank ('-') for Formatting Purposes
figSubtitle = ' ';

for i = 1:height(plotVars)
    scalarData = mapData.mean.(plotVars{i});
    
    figName = ['Time_Averaged_', plotVars{i}, '_Map_', num2str(nParcels - nParcelsRemoved)];
    
    contourlines = [];
    
    if strcmp(plotVars{i}, 'nParticles')
%         cLims = [0; single(7.4216729e+03)];
        cLims = [0; single(1.2171055e+04)];
    elseif strcmp(plotVars{i}, 'mass')
%         cLims = [0; single(6.7297306e-09)];
        cLims = [0; single(1.1430570e-08)];
    elseif strcmp(plotVars{i}, 'd32')
%         cLims = [0; 27.1185875];
        cLims = [0; 44.0660515];
    elseif strcmp(plotVars{i}, 'massRMS')
%         cLims = [0; single(4.2842041e-09)];
        cLims = [0; single(6.3777281e-09)];
    else
        continue;
    end

    [fig, planeNo] = plotPlanarScalarField(orientation, xLimsData, yLimsData, zLimsData, positionData, scalarData, ...
                                           mapPerim, nPlanes, planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           xDims, yDims, zDims, refPoint, figTitle, figSubtitle, cLims, ...
                                           xLimsPlot, yLimsPlot, zLimsPlot, normalise, false);
end
clear i;


%%

CoM = mapData.mean.CoM

totalParcels = sum(double(mapData.mean.nParticles))

totalMass = sum(double(mapData.mean.mass))

meanD = mean(double(mapData.mean.d32))

meanRMS = mean(double(mapData.mean.massRMS))