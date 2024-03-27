%% Planar Lagrangian Parcel Counts v1.0
% ----
% Lorem ipsum


%% Preamble

run preamble;

%#ok<*UNRCH>

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name
    
normDims = true; % Normalise Spatial Dimensions in Plots

figSave = false; % Save .fig File(s)

disp('=========================');
disp('Planar Parcel Counts v1.0');
disp('=========================');

disp(' ');
disp(' ');


%% Changelog

% v1.0 - Initial Commit


%% Initialise Case

[caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normLength] = initialiseCaseData(geoLoc);

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
        [dataID, LagProps, LagData, ~, ...
         ~, sampleInt, dataFormat] = initialiseLagData(saveLoc, caseFolder, campaignID, ...
                                                       caseID, cloudName, true, false, ...
                                                       false, timeDirs, deltaT, ...
                                                       timePrecision, maxNumCompThreads);
        
    case 'B'
        [dataID, LagProps, ~, LagData, ...
         ~, sampleInt, dataFormat] = initialiseLagData(saveLoc, caseFolder, campaignID, ...
                                                       caseID, cloudName, false, true, ...
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

if strcmp(campaignID, 'Windsor_fullScale')
    dLimsDefault = [20; 400]; % um
else
    dLimsDefault = [1; 147]; % um
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
dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2)), '_', dataFormat];

disp(' ');
disp(' ');


%% Generate Spray Maps

disp('Spray Mapping');
disp('--------------');

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

% Adjust Data Origin
if strcmp(campaignID, 'Windsor_Upstream_2023')
    
    for i = 1:nTimes
        
        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    clear i;
    
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
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);
        
        if strcmp(campaignID, 'Windsor_fullScale')
            basePoly = polybuffer(basePoly, -16e-3, 'jointType', 'square');
        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
            basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
        else
            basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
        end
        
        mapPerim = ones(height(basePoly.Vertices),3) * mapPerim(1,1);
        mapPerim(:,[2,3]) = basePoly.Vertices(:,[1,2]);

        if ~all(mapPerim(1,:) == mapPerim(end,:))
            mapPerim = [mapPerim; mapPerim(1,:)]; % Close Boundary
        end
        
        xLimsData = xDims(2);
        yLimsData = [min(mapPerim(:,2)); max(mapPerim(:,2))];
        zLimsData = [min(mapPerim(:,3)); max(mapPerim(:,3))];
        
    case 'B'
        mapPerim = [];
        
        xLimsData = double(LagData.positionCartesian{end}(1,1));
        yLimsData = [-0.5; 0.5] * normLength;
        zLimsData = [0; 0.5] * normLength;
        
end

disp(' ');

% Collate Particles of Interest
disp('    Collating Particles of Interest...');

% Initialise Progress Bar
wB = waitbar(0, 'Collating Particles of Interest', 'name', 'Progress');
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
reqFields = {'time'; 'positionCartesian'};

LagData = rmfield(LagData, LagFields(~ismember(LagFields, reqFields)));
LagProps = LagProps(ismember(LagProps, reqFields));

for i = 1:nTimes
    
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
if strcmp(campaignID, 'Windsor_fullScale')
    cellSize.target = 32e-3;
elseif strcmp(campaignID, 'Windsor_Upstream_2023')
    cellSize.target = 8e-3;
else
    cellSize.target = 8e-3;
end

% Adjust Uniform Cell Size to Fit Region of Interest
nPy = (diff(yLimsData) / cellSize.target) + 1;
nPz = (diff(zLimsData) / cellSize.target) + 1;

sizeY = diff(linspace(yLimsData(1), yLimsData(2), nPy));
sizeZ = diff(linspace(zLimsData(1), zLimsData(2), nPz));

cellSize.y = sizeY(1); clear sizeY;
cellSize.z = sizeZ(1); clear sizeZ;
cellSize.area = cellSize.y * cellSize.z;

% Generate Grid
[y, z] = ndgrid(linspace(yLimsData(1), yLimsData(2), nPy), ...
                linspace(zLimsData(1), zLimsData(2), nPz));

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
totalParcels = cellfun(@height, LagData.positionCartesian);
index = cell(nTimes,1); % Array Position of Nearest Mesh Node

xVals = unique(mapData.positionGrid(:,1));
yVals = unique(mapData.positionGrid(:,2));
zVals = unique(mapData.positionGrid(:,3));
positionCartesian = LagData.positionCartesian;
positionGrid = mapData.positionGrid;
parfor i = 1:nTimes
    
    if totalParcels(i) > 0
        index3D = zeros([totalParcels(i),3], 'uint32');
        positions = zeros([totalParcels(i),3]);
        
        for j = 1:totalParcels(i)
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

mapData.time = LagData.time;

% Generate Instantaneous Parcel Counts
disp('    Calculating Instantaneous Parcel Counts...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Parcel Counts', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';
dQ = parallel.pool.DataQueue;
afterEach(dQ, @parforWaitBar);

parforWaitBar(wB, nTimes);

% Perform Calculation
nParcels = cell(nTimes,1); nParcels(:) = {zeros([nCells,1])}; % Number of Parcels in Cell

parfor i = 1:nTimes
    
    if totalParcels(i) > 0
        
        for j = 1:totalParcels(i)
            nParcels{i}(index{i}(j)) = nParcels{i}(index{i}(j)) + 1;
        end
        
    end
    
    % Make Arrays Sparse
    nParcels{i} = sparse(nParcels{i});

    % Remove Unnecessary Data
    index{i} = -1;
    
    % Update Waitbar
    send(dQ, []);
end

delete(wB);

clear LagData index;

mapData.nParcels.inst = nParcels; clear nParcels;

disp(' ');

% Generate Instantaneous Parcel Counts
disp('    Calculating Time-Averaged Parcel Counts...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Parcel Counts', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.nParcels.mean = sparse(nCells,1);

for i = 1:nTimes
    mapData.nParcels.mean = mapData.nParcels.mean + mapData.nParcels.inst{i};
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

mapData.nParcels.mean = mapData.nParcels.mean / nTimes;

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
    
    selection = input('Plot Time-Averaged Map? [y/n]: ', 's');

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

if plotMean || plotInst
    
    % Normalise Coordinate System
    if normDims
        disp(' ');

        disp('Normalising Spatial Dimensions...');

        parts = fieldnames(geometry);
        for i = 1:height(parts)
            geometry.(parts{i}).vertices = geometry.(parts{i}).vertices / normLength;
        end
        clear i parts;

        xDims = xDims / normLength;
        yDims = yDims / normLength;
        zDims = zDims / normLength;

        mapPerim = mapPerim / normLength;

        xLimsData = xLimsData / normLength;
        yLimsData = yLimsData / normLength;
        zLimsData = zLimsData / normLength;

        cellSize.target = cellSize.target / normLength;
        cellSize.y = cellSize.y / normLength;
        cellSize.z = cellSize.z / normLength;
        cellSize.area = cellSize.area / (normLength^2);

        mapData.positionGrid = mapData.positionGrid / normLength;
    end
    
end

disp(' ');
disp(' ');


%% Present Spray Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotInst
    orientation = 'YZ';
    positionData = mapData.positionGrid;
    
    if strcmp(campaignID, 'Windsor_fullScale')
        spatialRes = 2e-3;
    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
        spatialRes = 0.5e-3;
    else
        spatialRes = 0.5e-3;
    end
    
    if normDims
        spatialRes = spatialRes / normLength;
    end
    
    % Offset Particles From Surface to Improve Visibility
    switch format
        
        case 'A'
            xLimsData = xLimsData + spatialRes;
            positionData(:,1) = xLimsData;
            
    end
    
    nPlanes = 1;
    planeNo = 1;
    cMap = flipud(viridis(32));
    refPoint = [];
    
    switch format

        case 'A'
            xLimsPlot = [0.3; 4.6257662];
            yLimsPlot = [-0.25; 0.25];
            zLimsPlot = [0; 0.4];
            
        case 'B'
            xLimsPlot = [0.3; 4.6257662];
            yLimsPlot = [-0.5; 0.5];
            zLimsPlot = [0; 0.5];
            
    end
    
    if ~normDims
        xLimsPlot = xLimsPlot * normLength;
        yLimsPlot = yLimsPlot * normLength;
        zLimsPlot = zLimsPlot * normLength;
    end
    
end

if plotMean
    
    disp('Presenting Time-Averaged Parcel Counts...');

    scalarData = full(mapData.nParcels.mean);

    switch format

        case 'A'
            figName = ['Average_Base_nParcels_', caseID];

        case 'B'
            figName = ['Average_', planeID, '_nParcels_', caseID];

    end
    
    contourlines = [];
    figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
    cLims = [0; max(scalarData)];

    [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                           xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                           planeNo, fig, figName, cMap, geometry, contourlines, ...
                                           refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                           zLimsPlot, normDims, figSave);
    
    disp(' ');
end

if plotInst
    disp('Presenting Instantaneous Parcel Counts...');

    contourlines = [];

    instMax = zeros([nTimes,1]);

    for j = 1:nTimes
        instMax(j) = full(max(mapData.nParcels.inst{j}));
    end

    cLims = [0; max(instMax)];

    figHold = fig;

    for j = startFrame:endFrame

        if j ~= startFrame
            clf(fig);
            fig = figHold;
        end

        scalarData = full(mapData.nParcels.inst{j});
        figTime = num2str(mapData.time(j), ['%.', num2str(timePrecision), 'f']);

        switch format

            case 'A'
                figName = ['Inst_Base_nParcels_T', erase(figTime, '.'), '_', caseID];

            case 'B'
                figName = ['Inst_', planeID, '_nParcels_T', erase(figTime, '.'), '_', caseID];
        end

        figTitle = ['{', figTime, ' \it{s}}'];

        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end
    clear j;
    
    disp(' ');
end

if ~plotMean && ~plotInst
    disp('Skipping Map Presentation...');
    
    disp(' ');
end


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