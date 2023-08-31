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

nProc = 4; % Number of Processors Used for Process-Based Parallelisation

fig = 0; % Initialise Figure Tracking
figHold = 0; % Enable Overwriting of Figures


%% Planar Lagrangian Contaminant Mapper v3.0

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

figSave = false; % Save .fig File(s)

normalise = true; % Normalisation of Dimensions

normDensity = 1;

disp('==============================');
disp('Planar Contaminant Mapper v3.0');
disp('==============================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit (Base Contamination and Far-Field Extraction Plane)
% v1.1 - Updated Name and Added Time-Averaging Functionality
% v2.0 - Rewrite, Accommodating New OpenFOAM Data Formats
% v2.1 - Improved Efficiency of Map Generation
% v2.2 - Added Support for Arithmetic and Sauter Mean Particle Diameters
% v2.3 - Changed to Thread-Based Parallelization to Reduce Memory Requirements
% v2.4 - Added Support for Full-Scale Windsor Model Simulations
% v2.5 - Changed Primary Output From Mass to Area Density
% v3.0 - Rewrite, Making Use of Sparse Arrays to Reduce Memory Requirements


%% Initialise Case

[caseFolder, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normalise, normLength] = initialiseCaseData(normalise);

disp(' ');
disp(' ');


%% Select Mapping Location

disp('Mapping Location');
disp('-----------------');

disp(' ');

disp('Possible Mapping Locations:');
disp('    A: Base Contamination');
disp('    B: Far-Field Spray Transport');

valid = false;
while ~valid
    disp(' ');
    selection = input('Select Mapping Location [A/B]: ', 's');

    if selection == 'a' | selection == 'A' %#ok<OR2>
        format = 'A';
        valid = true;
    elseif selection == 'b' | selection == 'B' %#ok<OR2>
        format = 'B';
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

disp(' ');
disp(' ');


%% Initialise Lagrangian Data

maxNumCompThreads(nProc);

switch format
    
    case 'A'
        [dataID, LagProps, ~, LagData, ~, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                              cloudName, false, true, ...
                                                                              false, timeDirs, deltaT, ...
                                                                              timePrecision, maxNumCompThreads);
                                                                                
    case 'B'
        [dataID, LagProps, LagData, ~, ~, sampleInterval] = initialiseLagData(saveLocation, caseFolder, caseID, ...
                                                                              cloudName, true, false, ...
                                                                              false, timeDirs, deltaT, ...
                                                                              timePrecision, maxNumCompThreads);

        % Select Plane of Interest
        planes = fieldnames(LagData);
        
        valid = false;
        while ~valid
            disp(' ');
            [index, valid] = listdlg('listSize', [300, 300], ...
                                     'selectionMode', 'single', ...
                                     'name', 'Select Plane of Interest', ...
                                     'listString', planes);
            
            if ~valid
                disp('WARNING: No Plane of Interest Selected');
            end
        
        end
        clear valid;

        LagData = LagData.(planes{index});
        planeID = erase(planes{index}, '.');
        
        clear planes;
        
        disp(['Plane of Interest: ', planeID]);
end

disp(' ');
disp(' ');


%% Select Mapping Options

disp('Mapping Options');
disp('----------------');

if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    dLimsDefault = [1; 147];
else
    dLimsDefault = [20; 400];
end

valid = false;
while ~valid
    disp(' ');
    selection = input('Filter Particle Diameters? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        dLims = dLimsDefault;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        dLims = zeros([2,1]);
        
        dLims(1) = inputD('Min');

        if dLims(1) == -1
            continue;
        end

        dLims(2) = inputD('Max');

        if dLims(2) == -1
            continue;
        end
        
        dLims = sort(dLims);
        dLims(1) = floor(dLims(1));
        dLims(2) = ceil(dLims(2));
        
        if dLims(2) < dLimsDefault(1) || dLims(1) > dLimsDefault(2)
            disp('        WARNING: No Lagrangian Data in Diameter Range');
            
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid dLimsDefault;

% Generate Dataset ID
if normalise
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2)), '_Norm'];
else
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2))];
end

disp(' ');
disp(' ');


%% Generate Contaminant Maps

disp('Contaminant Mapping');
disp('--------------------');

disp(' ');

disp('***********');
disp('  RUNNING ');

tic;

evalc('parpool(''threads'');');

%%%%

disp(' ');

disp('    Initialising...');

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
if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    
    for i = 1:nTimes
        
        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    clear i;
    
end

% Normalise Coordinate System
if normalise
    disp('        Normalising Coordinate System');
    
    % Initialise Progress Bar
    wB = waitbar(0, 'Normalising Coordinate System', 'name', 'Progress');
    wB.Children.Title.Interpreter = 'none';
    
    % Perform Normalisation
    for i = 1:nTimes
        
        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}  = round((LagData.positionCartesian{i} / normLength), ...
                                                  spacePrecision);
        end
        
        % Update Waitbar
        waitbar((i / nTimes), wB);
    end
    clear i;
    
    delete(wB);
end

% Specify Map Boundaries
switch format
    
    case 'A'

        % Identify Model Base
        parts = fieldnames(geometry);
        for i = 1:height(parts)
            
            if max(geometry.(parts{i}).vertices(:,1)) == xDims(2)
                parts = parts{i};
                break;
            end
            
            if i == height(parts)
                error('Mismatch Between ''xDims'' and Geometry Bounding Box')
            end
            
        end
        clear i;
    
        geoPoints = geometry.(parts).vertices;
        basePoints = geoPoints((geoPoints(:,1) == xDims(2)),:);
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 1);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
        basePoly = polybuffer(basePoly, -0.005, 'jointType', 'square');
        mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
        mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(mapPerim(1,:) == mapPerim(end,:))
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        xLimsData = xDims(2)
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case 'B'
        mapPerim = [];
    
        if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
            xLimsData = double(LagData.positionCartesian{end}(1,1));
            yLimsData = [-0.6264; 0.6264];
            zLimsData = [0; 0.6264];
        else
            xLimsData = double(LagData.positionCartesian{end}(1,1));
            yLimsData = [-2.5056; 2.5056];
            zLimsData = [0; 2.5056];
        end
        
        if normalise
            yLimsData = round((yLimsData / normLength), spacePrecision);
            zLimsData = round((zLimsData / normLength), spacePrecision);
        end
        
end

disp(' ');

% Collate Particles of Interest
disp('    Collating Particles of Interest...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Centre of Mass', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

index = cell(nTimes,1);

switch format

    case 'A'
        
        for i = 1:nTimes
            
            if ~isempty(LagData.positionCartesian{i})
                index{i} = find(((LagData.d{i} * 1e6) >= dLims(1)) & ...
                                ((LagData.d{i} * 1e6) <= dLims(2)) & ...
                                (round(LagData.positionCartesian{i}(:,1), spacePrecision) == single(round(xLimsData, spacePrecision))) & ...
                                (LagData.positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                                (LagData.positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,3) <= zLimsData(2)));
            end
            
            % Update Waitbar
            waitbar((i / nTimes), wB);
        end
        clear i;

    case 'B'
        
        for i = 1:nTimes
            
            if ~isempty(LagData.positionCartesian{i})
                index{i} = find(((LagData.d{i} * 1e6) >= dLims(1)) & ...
                                ((LagData.d{i} * 1e6) <= dLims(2)) & ...
                                (LagData.positionCartesian{i}(:,2) >= yLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,2) <= yLimsData(2)) & ...
                                (LagData.positionCartesian{i}(:,3) >= zLimsData(1)) & ...
                                (LagData.positionCartesian{i}(:,3) <= zLimsData(2)));
            end
            
            % Update Waitbar
            waitbar((i / nTimes), wB);
        end
        clear i;

end

delete(wB);

% Remove Unnecessary Data
disp('        Removing Unnecessary Data');

LagFields = fieldnames(LagData);
reqFields = {'time'; 'timeExact'; 'd'; 'nParticle'; 'positionCartesian'};

LagData = rmfield(LagData, LagFields(~ismember(LagFields, reqFields)));
LagProps = LagProps(ismember(LagProps, reqFields));

for i = 1:nTimes
    LagData.timeExact{i} = LagData.timeExact{i}(index{i});
    
    for j = 1:height(LagProps)
        LagData.(LagProps{j}){i} = LagData.(LagProps{j}){i}(index{i},:);
    end
    clear j;

end
clear i;

disp(' ');

% Generate Presentation Grid
disp('    Generating Presentation Grid...');

% Set Target Spatial Resolution
if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
    cellSize.target = 8e-3;
else
    cellSize.target = 32e-3;
end

% Adjust Uniform Cell Size to Fit Region of Interest
if normalise
    cellSize.target = round((cellSize.target / normLength), spacePrecision);
end

cellSize.x = cellSize.target;
cellSize.y = (yLimsData (2) - yLimsData (1)) / round(((yLimsData (2) - yLimsData (1)) / cellSize.target));
cellSize.z = (zLimsData (2) - zLimsData (1)) / round(((zLimsData (2) - zLimsData (1)) / cellSize.target));

cellSize.area = cellSize.y * cellSize.z;

[y, z] = ndgrid(yLimsData(1):cellSize.y:yLimsData(2), zLimsData(1):cellSize.z:zLimsData(2));

mapData.positionGrid = zeros([height(y(:)),3]);
mapData.positionGrid(:,1) = xLimsData;
mapData.positionGrid(:,(2:3)) = [y(:), z(:)]; clear y z;

nCells = height(mapData.positionGrid);

% Assign Particles to Grid Cells
disp('        Assigning Particles Grid Cells');

% Initialise Progress Bar
wB = waitbar(0, 'Assigning Particles to Grid Cells', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, nTimes);

% Perform Assignment
totalParticles = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Nearest Mesh Node

xVals = unique(mapData.positionGrid(:,1));
yVals = unique(mapData.positionGrid(:,2));
zVals = unique(mapData.positionGrid(:,3));
positionCartesian = LagData.positionCartesian;
positionGrid = mapData.positionGrid;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        index3D = zeros([totalParticles(i),3], 'uint32');
        positions = zeros([totalParticles(i),3]);
        
        for j = 1:totalParticles(i)
            [~, index3D(j,1)] = min(abs(positionCartesian{i}(j,1) - xVals));
            [~, index3D(j,2)] = min(abs(positionCartesian{i}(j,2) - yVals));
            [~, index3D(j,3)] = min(abs(positionCartesian{i}(j,3) - zVals));
            
            positions(j,:) = [xVals(index3D(j,1)), yVals(index3D(j,2)), zVals(index3D(j,3))];
        end
        
        [~, index{i}] = ismember(positions, positionGrid, 'rows');
    end
    
    % Remove Unnecessary Data
    positionCartesian{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear xVals yVals zVals positionCartesian positionGrid;

delete(wB);

disp(' ');

% Generate Instantaneous Contaminant Maps
disp('    Generating Instantaneous Contaminant Maps...');

mapData.time = LagData.time;

% Calculate Instantaneous Field Variables
disp('        Calculating Instantaneous Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, nTimes);

% Perform Calculation
nParticles = cell(nTimes,1); nParticles(:) = {zeros([nCells,1])}; % Number of Particles in Cell
density = nParticles; % Spray Density in Cell
d32 = nParticles; % Sauter Mean Particle Diameter in Cell
d10 = nParticles; % Arithmetic Mean Particle Diameter in Cell

d_tmp = zeros([nCells,1]);
nParticle = cellfun(@double, LagData.nParticle, 'uniformOutput', false);
d = cellfun(@double, LagData.d, 'uniformOutput', false);
cellArea = cellSize.area;
clear LagData;
parfor i = 1:nTimes
    
    if totalParticles(i) > 0
        d30 = d_tmp;
        d20 = d_tmp;
        
        for j = 1:totalParticles(i)
            nParticles{i}(index{i}(j)) = nParticles{i}(index{i}(j)) + ...
                                         nParticle{i}(j);
            
            density{i}(index{i}(j)) = density{i}(index{i}(j)) + ...
                                      (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
            
            d30(index{i}(j)) = d30(index{i}(j)) + ...
                               (nParticle{i}(j) * (d{i}(j)^3));
            
            d20(index{i}(j)) = d20(index{i}(j)) + ...
                               (nParticle{i}(j) * (d{i}(j)^2));
            
            d10{i}(index{i}(j)) = d10{i}(index{i}(j)) + ...
                                  (nParticle{i}(j) * d{i}(j));
        end
        
        % Calculate Derived Variables
        density{i} = ((1000 * density{i}) / cellArea) / normDensity;
        d32{i} = (d30 ./ d20) * 1e6;
        d10{i} = (d10{i} ./ nParticles{i}) * 1e6;
        
        % Set Empty Cells Back to Zero
        d32{i}(isnan(d32{i})) = 0;
        d10{i}(isnan(d10{i})) = 0;
    end
    
    % Make Arrays Sparse
    nParticles{i} = sparse(nParticles{i});
    density{i} = sparse(density{i});
    d32{i} = sparse(d32{i});
    d10{i} = sparse(d10{i});

    % Remove Unnecessary Data
    index{i} = [];
    nParticle{i} = [];
    d{i} = [];
    
    % Update Waitbar
    send(dQ, []);
end
clear index nParticle d d_tmp;

delete(wB);

mapData.inst.nParticles = nParticles; clear nParticles;
mapData.inst.density = density; clear density;
mapData.inst.d32 = d32; clear d32;
mapData.inst.d10 = d10; clear d10;

% Calculate Instantaneous Centre of Mass
disp('        Calculating Instantaneous Centre of Mass');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Centre of Mass', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.inst.CoM = cell(nTimes,1); mapData.inst.CoM(:) = {zeros([1,3], 'single')};

for i = 1:nTimes
    mapData.inst.CoM{i}(1) = mapData.positionGrid(1,1);
    mapData.inst.CoM{i}(2) = sum(mapData.inst.density{i} .* mapData.positionGrid(:,2)) / sum(mapData.inst.density{i});
    mapData.inst.CoM{i}(3) = sum(mapData.inst.density{i} .* mapData.positionGrid(:,3)) / sum(mapData.inst.density{i});
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

disp(' ');

% Generate Time-Averaged Contaminant Map
disp('    Generating Time-Averaged Contaminant Maps...');

% Calculate Time-Averaged Field Variables
disp('        Calculating Instantaneous Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.mean.nParticles = zeros([nCells,1]);
mapData.mean.density = mapData.mean.nParticles;
mapData.mean.d32 = mapData.mean.nParticles;
mapData.mean.d10 = mapData.mean.nParticles;

for i = 1:nTimes
    mapData.mean.nParticles = mapData.mean.nParticles + mapData.inst.nParticles{i};
    mapData.mean.density = mapData.mean.density + mapData.inst.density{i};
    mapData.mean.d32 = mapData.mean.d32 + mapData.inst.d32{i};
    mapData.mean.d10 = mapData.mean.d10 + mapData.inst.d10{i};
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

mapData.mean.nParticles = sparse(mapData.mean.nParticles / nTimes);
mapData.mean.density = sparse(mapData.mean.density / nTimes);
mapData.mean.d32 = sparse(mapData.mean.d32 / nTimes);
mapData.mean.d10 = sparse(mapData.mean.d10 / nTimes);

% Calculate Time-Averaged Centre of Mass
disp('        Calculating Time-Averaged Centre of Mass');

mapData.mean.CoM = zeros([1,3], 'single');
    
mapData.mean.CoM(1) = mapData.positionGrid(1,1);
mapData.mean.CoM(2) = sum(mapData.mean.density .* mapData.positionGrid(:,2)) / sum(mapData.mean.density);
mapData.mean.CoM(3) = sum(mapData.mean.density .* mapData.positionGrid(:,3)) / sum(mapData.mean.density);

%%%%

evalc('delete(gcp(''nocreate''));');

executionTime = toc;

disp(' ');

disp(['    Run Time: ', num2str(executionTime), 's']);

disp(' ');

disp('  SUCCESS  ');
disp('***********');

disp(' ');
disp(' ');


%% Select Presentation Options

disp('Presentation Options');
disp('---------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Time-Averaged Map(s)? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotMean = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotMean = true;
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

valid = false;
while ~valid
    disp(' ');
    selection = input('Plot Instantaneous Maps? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotInst = false;
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotInst = true;
        
        startFrame = inputFrames(nTimes, 'Start');
        
        if startFrame == -1
            continue;
        end
        
        endFrame = inputFrames(nTimes, 'End');
        
        if endFrame == -1
            continue;
        elseif endFrame < startFrame
            disp('        WARNING: Invalid Time Format (''endFrame'' Precedes ''startFrame'')');
            
            continue;
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;

if plotInst || plotMean
    
    % Select Variable(s) of Interest
    plotVars = fieldnames(mapData.mean);
    plotVars = plotVars(1:(end - 1));

    valid = false;
    while ~valid
        disp(' ');
        [index, valid] = listdlg('listSize', [300, 300], ...
                                 'selectionMode', 'multiple', ...
                                 'name', 'Select Variable(s) to Plot', ...
                                 'listString', plotVars);

        if ~valid
            disp('WARNING: No Mapping Variables Selected');
        end
        
    end
    clear valid;

    plotVars = plotVars(index);
    
    disp(['Plotting ', num2str(height(plotVars)), ' Variable(s) of Interest']);
end

disp(' ');
disp(' ');


%% Present Contaminant Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotInst || plotMean
    
    % Define Plot Limits
    switch format

        case 'A'
            orientation = 'YZ';

            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                xLimsPlot = [0.31875; 1.43075];
                yLimsPlot = [-0.2445; 0.2445];
                zLimsPlot = [0; 0.389];
            else
                xLimsPlot = [1.275; 1.43075];
                yLimsPlot = [-0.978; 0.978];
                zLimsPlot = [0; 1.538];
            end

        case 'B'
            orientation = 'YZ';

            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                xLimsPlot = [0.31875; 2.73575];
                yLimsPlot = [-0.522; 0.522];
                zLimsPlot = [0; 0.522];
            else
                xLimsPlot = [1.275; 19.295];
                yLimsPlot = [-2.088; 2.088];
                zLimsPlot = [0; 2.088];
            end

    end
    
    if normalise
        xLimsPlot = round((xLimsPlot / normLength), spacePrecision);
        yLimsPlot = round((yLimsPlot / normLength), spacePrecision);
        zLimsPlot = round((zLimsPlot / normLength), spacePrecision);
    end
    
    spatialRes = cellSize.target / 8;
    positionData = mapData.positionGrid;
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    figTitle = '-'; % Leave Blank ('-') for Formatting Purposes

end

if plotMean
    
    for i = 1:height(plotVars)
        disp(['    Presenting Time-Averaged ''', plotVars{i}, ''' Data...']);
        
        scalarData = full(mapData.mean.(plotVars{i}));
        
        switch format
            
            case 'A'
                figName = ['Time_Averaged_Base_', plotVars{i}, '_Map'];
                
            case 'B'
                figName = ['Time_Averaged_', planeID, '_', plotVars{i}, '_Map'];
                
        end
        
        if strcmp(plotVars{i}, 'density')
            
            if normDensity == 1
                contourlines = [0.02; 0.02] * max(scalarData);
            else
                contourLines = [0.02; 0.02];
            end
            
        else
            contourlines = [];
        end
        
        figSubtitle = ' ';
        
        if any(strcmp(plotVars{i}, {'d10', 'd32'}))
            
            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                cLims = [0; 150];
            else
                cLims = [0; 400];
            end
            
        elseif strcmp(plotVars{i}, 'density') && normDensity == 1
            cLims = [0; max(scalarData)];
        else
            cLims = [0; 1.2];
        end
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
                                               yLimsPlot, zLimsPlot, normalise, figSave);
    end
    clear i;
    
    disp(' ');
end

if plotInst
    
    for i = 1:height(plotVars)
        disp(['    Presenting Instantaneous ''', plotVars{i}, ''' Data...']);
        
        if strcmp(plotVars{i}, 'density')
            
            if normDensity == 1
                contourlines = [0.02; 0.02] * max(cellfun(@max, mapData.inst.(plotVars{i})));
            else
                contourLines = [0.02; 0.02] * 1;
            end
            
        else
            contourlines = [];
        end
        
        if any(strcmp(plotVars{i}, {'d10', 'd32'}))
            
            if contains(caseID, 'Run_Test') || (contains(caseID, 'Windsor') && contains(caseID, 'Upstream'))
                cLims = [0; 150];
            else
                cLims = [0; 400];
            end
            
        elseif strcmp(plotVars{i}, 'density') && normDensity == 1
            cLims = [0; max(cellfun(@max, mapData.inst.(plotVars{i})))];
        else
            cLims = [0; 1];
        end
        
        figHold = fig;
        
        for j = startFrame:endFrame
            
            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            scalarData = full(mapData.inst.(plotVars{i}){j});
            figTime = num2str(mapData.time(j), ['%.', num2str(timePrecision), 'f']);
            
            switch format
                
                case 'A'
                    figName = ['Instantaneous_Base_', plotVars{i}, '_Map_T', erase(figTime, '.')];
                
                case 'B'
                    figName = ['Instantaneous_', planeID, '_', plotVars{i}, '_Map_T', erase(figTime, '.')];
            end
            
            figSubtitle = [figTime, ' \it{s}'];
        
            [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                                   xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                                   planeNo, fig, figName, cMap, geometry, contourlines, ...
                                                   refPoint, figTitle, figSubtitle, cLims, xLimsPlot, ...
                                                   yLimsPlot, zLimsPlot, normalise, figSave);
        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if ~plotMean && ~plotInst
    disp('    Skipping Map Presentation');

    disp(' ');
end

disp(' ');


%% Save Map Data

disp('Data Save Options');
disp('------------------');

valid = false;
while ~valid
    disp(' ');
    selection = input('Save Data for Future Use? [y/n]: ', 's');
    
    if selection == 'n' | selection == 'N' %#ok<OR2>
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        
        switch format
            
            case 'A'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base'], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base']);
                end
                
            case 'B'
                
                if ~exist([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID], 'dir')
                    mkdir([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/base/', dataID, '.mat'], ...
                      'caseID', 'dataID', 'mapData', 'cellSize', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID, '/', dataID, '.mat']);
                save([saveLocation, '/Numerical/MATLAB/planarContaminantMap/', caseID, '/', planeID, '/', dataID, '.mat'], ...
                      'caseID', 'planeID', 'dataID', 'mapData', 'cellSize', 'sampleInterval', 'timePrecision', 'dLims', 'normalise', '-v7.3', '-noCompression');
                disp('        Success');
        
        end
        
        valid = true;
    else
        disp('    WARNING: Invalid Entry');
    end

end
clear valid;


%% Local Functions

function D = inputD(type)

    D = str2double(input(['    ', type, ' Diameter of Interest [', char(956), 'm]: '], 's'));
    
    if isnan(D) || length(D) > 1 || D < 1
        disp('        WARNING: Invalid Entry');
        
        D = -1;
    end
    
end


function frameNo = inputFrames(Nt, type)

    frameNo = str2double(input(['    Input Desired ', type, ' Frame [1-', num2str(Nt), ']: '], 's'));
    
    if isnan(frameNo) || frameNo < 1 || frameNo > Nt
        disp('        WARNING: Invalid Entry');
        
        frameNo = -1;
    end

end