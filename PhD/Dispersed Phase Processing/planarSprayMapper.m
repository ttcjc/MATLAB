%% Planar Lagrangian Spray Mapper v4.0
% ----
% Load, Process and Present Planar Lagrangian Data Acquired Using OpenFOAM v7


%% Preamble

run preamble;

cloudName = 'kinematicCloud'; % OpenFOAM Cloud Name

normDims = true; % Normalise Spatial Dimensions

figSave = false; % Save .fig File(s)

disp('========================');
disp('Planar Spray Mapper v4.0');
disp('========================');

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
% v3.1 - Minor Update to Shift Preamble Into Separate Script
% v4.0 - Update To Include Calculation of Fluctuating Field Variables


%% Initialise Case

[caseFolder, campaignID, caseID, timeDirs, deltaT, timePrecision, geometry, ...
 xDims, yDims, zDims, spacePrecision, normDims, normLength] = initialiseCaseData(normDims);

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
    dLimsDefault = [1; 147];
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
if normDims
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2)), '_', dataFormat, '_normDims'];
else
    dataID = [dataID, '_D', num2str(dLims(1)), '_D', num2str(dLims(2)), '_', dataFormat];
end

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

% Shift Data Origin
if contains(caseID, 'Run_Test') || strcmp(campaignID, 'Windsor_Upstream_2023')
    
    for i = 1:nTimes
        
        if ~isempty(LagData.positionCartesian{i})
            LagData.positionCartesian{i}(:,1) = LagData.positionCartesian{i}(:,1) + 1.325;
        end
        
    end
    clear i;
    
end

% Normalise Coordinate System
if normDims
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
        
        mapPerim = boundary(basePoints(:,2), basePoints(:,3), 0.95);
        mapPerim = basePoints(mapPerim,:);
        basePoly = polyshape(mapPerim(:,2), mapPerim(:,3), 'keepCollinearPoints', true);

        if normDims
            basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
        else
            
            if strcmp(campaignID, 'Windsor_fullScale')
                basePoly = polybuffer(basePoly, -16e-3, 'jointType', 'square');
            else
                basePoly = polybuffer(basePoly, -4e-3, 'jointType', 'square');
            end
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
        yLimsData = [-0.5; 0.5];
        zLimsData = [0; 0.5];
        
        if ~normDims

            if strcmp(campaignID, 'Windsor_fullScale')
                yLimsData = round((yLimsData * 4.176), spacePrecision);
                zLimsData = round((zLimsData * 4.176), spacePrecision);
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                yLimsData = round((yLimsData * 1.044), spacePrecision);
                zLimsData = round((zLimsData * 1.044), spacePrecision);
            end

        end
        
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
if normDims
    cellSize.target = 8e-3;
else
    
    if strcmp(campaignID, 'Windsor_fullScale')
        cellSize.target = 32e-3;
    else
        cellSize.target = 8e-3;
    end
    
end

% Adjust Uniform Cell Size to Fit Region of Interest
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

% Generate Instantaneous Spray Maps
disp('    Generating Instantaneous Spray Maps...');

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
areaDensity = nParticles; % Spray Density in Cell
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
            
            areaDensity{i}(index{i}(j)) = areaDensity{i}(index{i}(j)) + ...
                                          (nParticle{i}(j) * ((1 / 12) * tau * (d{i}(j)^3)));
            
            d30(index{i}(j)) = d30(index{i}(j)) + ...
                               (nParticle{i}(j) * (d{i}(j)^3));
            
            d20(index{i}(j)) = d20(index{i}(j)) + ...
                               (nParticle{i}(j) * (d{i}(j)^2));
            
            d10{i}(index{i}(j)) = d10{i}(index{i}(j)) + ...
                                  (nParticle{i}(j) * d{i}(j));
        end
        
        % Calculate Derived Variables
        areaDensity{i} = ((1000 * areaDensity{i}) / cellArea);
        d32{i} = (d30 ./ d20) * 1e6;
        d10{i} = (d10{i} ./ nParticles{i}) * 1e6;
        
        % Set Empty Cells Back to Zero
        d32{i}(isnan(d32{i})) = 0;
        d10{i}(isnan(d10{i})) = 0;
    end
    
    % Make Arrays Sparse
    nParticles{i} = sparse(nParticles{i});
    areaDensity{i} = sparse(areaDensity{i});
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

mapData.nParticles.inst = nParticles; clear nParticles;
mapData.areaDensity.inst = areaDensity; clear areaDensity;
mapData.d32.inst = d32; clear d32;
mapData.d10.inst = d10; clear d10;

% Calculate Instantaneous Centre of Mass
disp('        Calculating Instantaneous Centre of Mass');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Centre of Mass', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.CoM.inst = cell(nTimes,1); mapData.CoM.inst(:) = {zeros([1,3], 'single')};

for i = 1:nTimes
    mapData.CoM.inst{i}(1) = mapData.positionGrid(1,1);
    mapData.CoM.inst{i}(2) = full(sum(mapData.areaDensity.inst{i} .* mapData.positionGrid(:,2)) / ...
                             sum(mapData.areaDensity.inst{i}));
    mapData.CoM.inst{i}(3) = full(sum(mapData.areaDensity.inst{i} .* mapData.positionGrid(:,3)) / ...
                             sum(mapData.areaDensity.inst{i}));
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

disp(' ');

% Generate Time-Averaged Spray Map
disp('    Generating Time-Averaged Spray Maps...');

% Calculate Time-Averaged Field Variables
disp('        Calculating Time-Averaged Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Time-Averaged Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.nParticles.mean = sparse(nCells,1);
mapData.areaDensity.mean = mapData.nParticles.mean;
mapData.d32.mean = mapData.nParticles.mean;
mapData.d10.mean = mapData.nParticles.mean;

for i = 1:nTimes
    mapData.nParticles.mean = mapData.nParticles.mean + mapData.nParticles.inst{i};
    mapData.areaDensity.mean = mapData.areaDensity.mean + mapData.areaDensity.inst{i};
    mapData.d32.mean = mapData.d32.mean + mapData.d32.inst{i};
    mapData.d10.mean = mapData.d10.mean + mapData.d10.inst{i};
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

mapData.nParticles.mean = mapData.nParticles.mean / nTimes;
mapData.areaDensity.mean = mapData.areaDensity.mean / nTimes;
mapData.d32.mean = mapData.d32.mean / nTimes;
mapData.d10.mean = mapData.d10.mean / nTimes;

% Calculate Time-Averaged Centre of Mass
disp('        Calculating Time-Averaged Centre of Mass');

mapData.CoM.mean = zeros([1,3], 'single');
    
mapData.CoM.mean(1) = mapData.positionGrid(1,1);
mapData.CoM.mean(2) = full(sum(mapData.areaDensity.mean .* mapData.positionGrid(:,2)) / ...
                      sum(mapData.areaDensity.mean));
mapData.CoM.mean(3) = full(sum(mapData.areaDensity.mean .* mapData.positionGrid(:,3)) / ...
                      sum(mapData.areaDensity.mean));

disp(' ');

% Generate Time-Averaged Spray Map
disp('    Generating Fluctuating Spray Maps...');

% Calculate Instantaneous Field Fluctuations
disp('        Calculating Instantaneous Field Fluctuations...');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating Instantaneous Field Fluctuations', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.nParticles.prime = mapData.nParticles.inst;
mapData.areaDensity.prime = mapData.areaDensity.inst;
mapData.d32.prime = mapData.d32.inst;
mapData.d10.prime = mapData.d10.inst;

for i = 1:nTimes
    mapData.nParticles.prime{i} = mapData.nParticles.prime{i} - mapData.nParticles.mean;
    mapData.areaDensity.prime{i} = mapData.areaDensity.prime{i} - mapData.areaDensity.mean;
    mapData.d32.prime{i} = mapData.d32.prime{i} - mapData.d32.mean;
    mapData.d10.prime{i} = mapData.d10.prime{i} - mapData.d10.mean;
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

% Calculate RMS of Field Variables
disp('        Calculating RMS of Field Variables');

% Initialise Progress Bar
wB = waitbar(0, 'Calculating RMS of Field Variables', 'name', 'Progress');
wB.Children.Title.Interpreter = 'none';

% Perform Calculation
mapData.nParticles.RMS = sparse(nCells,1);
mapData.areaDensity.RMS = mapData.nParticles.RMS;
mapData.d32.RMS = mapData.nParticles.RMS;
mapData.d10.RMS = mapData.nParticles.RMS;

for i = 1:nTimes
    mapData.nParticles.RMS = mapData.nParticles.RMS + mapData.nParticles.prime{i}.^2;
    mapData.areaDensity.RMS = mapData.areaDensity.RMS + mapData.areaDensity.prime{i}.^2;
    mapData.d32.RMS = mapData.d32.RMS + mapData.d32.prime{i}.^2;
    mapData.d10.RMS = mapData.d10.RMS + mapData.d10.prime{i}.^2;
    
    % Update Waitbar
    waitbar((i / nTimes), wB);
end
clear i;

delete(wB);

mapData.nParticles.RMS = sqrt((1 / nTimes) * mapData.nParticles.RMS);
mapData.areaDensity.RMS = sqrt((1 / nTimes) * mapData.areaDensity.RMS);
mapData.d32.RMS = sqrt((1 / nTimes) * mapData.d32.RMS);
mapData.d10.RMS = sqrt((1 / nTimes) * mapData.d10.RMS);

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
    
    selection = input('Plot RMS Map? [y/n]: ', 's');

    if selection == 'n' | selection == 'N' %#ok<OR2>
        plotRMS = false;
        
        valid = true;
    elseif selection == 'y' | selection == 'Y' %#ok<OR2>
        plotRMS = true;
        
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

if plotMean || plotRMS || plotInst
    disp(' ');
    
    % Select Variable(s) of Interest
    disp('Select Variable(s) to Plot...');

    mapDataVars = fieldnames(mapData);
    nonFieldVars = {'positionGrid'; 'time'; 'CoM'};
    fieldVars = setdiff(mapDataVars, nonFieldVars);
    clear mapDataVars nonFieldVars;

    valid = false;
    while ~valid
        [index, valid] = listdlg('listSize', [300, 300], ...
                                 'selectionMode', 'multiple', ...
                                 'name', 'Select Variable(s) to Plot', ...
                                 'listString', fieldVars);

        if ~valid
            disp('    WARNING: No Mapping Variable Selected');
        end

    end
    clear valid;

    plotVars = fieldVars(index); clear fieldVars;
    
    disp(['Plotting ', num2str(height(plotVars)), ' Variable(s) of Interest']);
end

disp(' ');
disp(' ');


%% Present Spray Maps

disp('Map Presentation');
disp('-----------------');

disp(' ');

if plotMean || plotRMS || plotInst
    orientation = 'YZ';
    positionData = mapData.positionGrid;
    
    if normDims
        spatialRes = 0.5e-3;
    else

        if strcmp(campaignID, 'Windsor_fullScale')
            spatialRes = 2e-3;
        else
            spatialRes = 0.5e-3;
        end

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

        if strcmp(campaignID, 'Windsor_fullScale')
            xLimsPlot = xLimsPlot * 4.176;
            yLimsPlot = yLimsPlot * 4.176;
            zLimsPlot = zLimsPlot * 4.176;
        elseif strcmp(campaignID, 'Windsor_Upstream_2023')
            xLimsPlot = xLimsPlot * 1.044;
            yLimsPlot = yLimsPlot * 1.044;
            zLimsPlot = zLimsPlot * 1.044;
        end

    end
    
end

if plotMean
    
    for i = 1:height(plotVars)
        disp(['    Presenting Time-Averaged ''', plotVars{i}, ''' Data...']);
        
        scalarData = full(mapData.(plotVars{i}).mean);
        
        switch format
            
            case 'A'
                figName = ['Average_Base_', plotVars{i}, '_', caseID];
                
            case 'B'
                figName = ['Average_', planeID, '_', plotVars{i}, '_', caseID];
                
        end
        
        switch format
            
            case 'A'
                contourlines = [];
                
            case 'B'
                
                if strcmp(plotVars{i}, 'areaDensity')
                    
                    if strcmp(campaignID, 'Windsor_fullScale')
%                         contourlines = [0.02; 0.02] * 15.5e-3; % Coupled
                        contourlines = [0.02; 0.02] * 2.6; % Uncoupled
                    elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                        contourlines = [0.02; 0.02] * 0.1e-3;
                    end
                    
                end
                
                if ~exist('contourlines', 'var')
                    contourlines = [];
                end
                
        end
        
        figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
        
        if any(strcmp(plotVars{i}, {'d10', 'd32'}))
            
            if strcmp(campaignID, 'Windsor_fullScale')
                cLims = [0; 400];
            else
                cLims = [0; 150];
            end
            
        elseif strcmp(plotVars{i}, 'areaDensity')
            
            if strcmp(campaignID, 'Windsor_fullScale')
                
                switch format

                    case 'A'
%                         cLims = [0; 7.5e-3]; % Coupled
                        cLims = [0; 55e-3]; % Uncoupled
                        

                    case 'B'
%                         cLims = [0; 15.5]; % Coupled
                        cLims = [0; 2.6]; % Uncoupled

                end
                
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
            
                switch format

                    case 'A'
                        cLims = [0; 5.6e-6];

                    case 'B'
                        cLims = [0; 0.1e-3];

                end
                
            end
            
        end
            
        if ~exist('cLims', 'var')
            cLims = [0; max(scalarData)];
        end
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end
    clear i;
    
    disp(' ');
end

if plotRMS
    
    for i = 1:height(plotVars)
        disp(['    Presenting RMS of ''', plotVars{i}, ''' Data...']);
        
        scalarData = full(mapData.(plotVars{i}).RMS) / mean(mapData.(plotVars{i}).mean(mapData.(plotVars{i}).mean > 0));
        
        switch format
            
            case 'A'
                figName = ['RMS_Base_', plotVars{i}, '_', caseID];
                
            case 'B'
                figName = ['RMS_', planeID, '_', plotVars{i}, '_', caseID];
                
        end
        
        contourlines = [];
        figTitle = '{ }'; % Leave Blank ('{ }') for Formatting Purposes
        
        if strcmp(plotVars{i}, 'areaDensity')
            
            if strcmp(campaignID, 'Windsor_fullScale')
                
                switch format

                    case 'A'
%                         cLims = [0; 16]; % Coupled
                        cLims = [0; 6];

                    case 'B'
%                         cLims = [0; 12.3]; % Coupled
                        cLims = [0; 4.4]; % Uncoupled

                end
                
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                
                switch format

                    case 'A'
                        cLims = [0; 4.1];

                    case 'B'
                        cLims = [0; 6.8];

                end
                
            end
            
        end
            
        if ~exist('cLims', 'var')
            cLims = [0; max(scalarData)];
        end
        
        [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                               xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                               planeNo, fig, figName, cMap, geometry, contourlines, ...
                                               refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                               zLimsPlot, normDims, figSave);
    end
    clear i;
    
    disp(' ');
end

if plotInst
    
    for i = 1:height(plotVars)
        disp(['    Presenting Instantaneous ''', plotVars{i}, ''' Data...']);
        
        contourlines = [];
        
        if any(strcmp(plotVars{i}, {'d10', 'd32'}))
            
            if strcmp(campaignID, 'Windsor_fullScale')
                cLims = [0; 400];
            else
                cLims = [0; 150];
            end
            
        elseif strcmp(plotVars{i}, 'areaDensity')
            
            if strcmp(campaignID, 'Windsor_fullScale')
                
                switch format

                    case 'A'
%                         cLims = [0; 0.5]; % Coupled
                        cLims = [0; 1.4]; % Uncoupled

                    case 'B'
%                         cLims = [0; 75]; % Coupled
                        cLims = [0; 16]; % Uncoupled

                end
                
            elseif strcmp(campaignID, 'Windsor_Upstream_2023')
                
                switch format

                    case 'A'
                        cLims = [0; 0.18-3];

                    case 'B'
                        cLims = [0; 0.7e-3];

                end
                
            end
            
        end
            
        if ~exist('cLims', 'var')
            instMax = 0;
            
            for j = 1:nTimes
                instMax = max(instMax, max(mapData.(plotVars{i}).inst{j}));
            end
                
            cLims = full([0; instMax]);
        end
        
        figHold = fig;
        
        for j = startFrame:endFrame
            
            if j ~= startFrame
                clf(fig);
                fig = figHold;
            end
            
            scalarData = full(mapData.(plotVars{i}).inst{j});
            figTime = num2str(mapData.time(j), ['%.', num2str(timePrecision), 'f']);
            
            switch format
                
                case 'A'
                    figName = ['Inst_Base_', plotVars{i}, '_T', erase(figTime, '.'), '_', caseID];
                
                case 'B'
                    figName = ['Inst_', planeID, '_', plotVars{i}, '_T', erase(figTime, '.'), '_', caseID];
            end
            
            figTitle = ['{', figTime, ' \it{s}}'];
        
            [fig, planeNo] = plotPlanarScalarField(orientation, positionData, scalarData, spatialRes, ...
                                                   xLimsData, yLimsData, zLimsData, mapPerim, nPlanes, ...
                                                   planeNo, fig, figName, cMap, geometry, contourlines, ...
                                                   refPoint, figTitle, cLims, xLimsPlot, yLimsPlot, ...
                                                   zLimsPlot, normDims, figSave);
        end
        clear j;
        
    end
    clear i;
    
    disp(' ');
end

if ~plotMean && ~ plotRMS && ~plotInst
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
                
                if ~exist([saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/base'], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/base']);
                end
                
            case 'B'
                
                if ~exist([saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/', planeID], 'dir')
                    mkdir([saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/', planeID]);
                end
                
        end
        
        switch format
            
            case 'A'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/base/', dataID, '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/base/', dataID, '.mat'], ...
                     'campaignID', 'caseID', 'dataID', 'mapData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', 'dataFormat', 'normDims', '-v7.3', '-noCompression');
                disp('        Success');
                 
            case 'B'
                disp(['    Saving to: ', saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/', planeID, '/', dataID, '.mat']);
                save([saveLoc, '/Numerical/MATLAB/planarSprayMap/', campaignID, '/', caseID, '/', planeID, '/', dataID, '.mat'], ...
                     'campaignID', 'caseID', 'planeID', 'dataID', 'mapData', 'cellSize', 'sampleInt', 'timePrecision', 'dLims', 'dataFormat', 'normDims', '-v7.3', '-noCompression');
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
